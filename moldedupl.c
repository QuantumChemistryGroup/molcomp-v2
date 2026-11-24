#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <omp.h>

#include "mc_types.h"
#include "mc_atominfo.h"
#include "mc_getconnect.h"
#include "mc_io.h"
#include "mc_core.h"

enum {
    MC_NAME_SZ = 256,
    MC_MODE_SZ = 32,
    MC_LINE_SZ = 512,
    MC_PATH_SZ = 512,
    MC_COPY_SZ = 8192,
    MC_DESC_SZ = 4,
    MC_OUT_SZ  = 3
};

typedef struct {
    char name[MC_NAME_SZ];
    long int orig_idx;
    long file_offset;
    double *xyz;
    double desc[MC_DESC_SZ];
    int ndesc;
    bool removed;
} ConfItem;

typedef struct {
    int idx1;
    int idx2;
    double out[MC_OUT_SZ];
    double dE_kcal;
} DupRecord;

typedef struct {
    char mode[MC_MODE_SZ];
    McOptions mc_opts;
    double value;
    double rmsd_thld;
    double maxdist_thld;
    bool readtopo;
} DedupConfig;

static size_t get_cgroup_limit_bytes(void) {
    size_t lim = 0;
    FILE *f = fopen("/sys/fs/cgroup/memory.max", "r");
    if (f) {
        char buf[64] = "";
        if (fgets(buf, sizeof(buf), f) && strncmp(buf, "max", 3) != 0) {
            unsigned long long v = strtoull(buf, NULL, 10);
            if (v > 0 && v < (1ULL << 62)) lim = (size_t)v;
        }
        fclose(f);
    }
    if (lim == 0) {
        f = fopen("/sys/fs/cgroup/memory/memory.limit_in_bytes", "r");
        if (f) {
            unsigned long long v = 0;
            if (fscanf(f, "%llu", &v) == 1 && v > 0 && v < (1ULL << 62)) lim = (size_t)v;
            fclose(f);
        }
    }
    return lim;
}

static size_t get_mem_available_bytes(void) {
    FILE *f = fopen("/proc/meminfo", "r");
    if (!f) return 0;
    char buf[256];
    unsigned long long kb = 0;
    while (fgets(buf, sizeof(buf), f)) {
        if (sscanf(buf, "MemAvailable: %llu kB", &kb) == 1) break;
    }
    fclose(f);
    size_t avail_ = (size_t)(kb * 1024ULL);
size_t lim_ = get_cgroup_limit_bytes();
if (lim_ > 0 && (avail_ == 0 || lim_ < avail_)) avail_ = lim_;
return avail_;
}

typedef struct {
    uint64_t *w;
    size_t nwords;
    long int n;
    int active;
} PairSet;

static size_t pairset_words(long int n) {
    if (n < 2) return 0;
    uint64_t np = ((uint64_t)n * (uint64_t)(n - 1)) / 2ULL;
    return (size_t)((np + 63ULL) / 64ULL);
}

static void pairset_init(PairSet *ps, long int n, size_t max_bytes) {
    ps->w = NULL;
    ps->nwords = 0;
    ps->n = n;
    ps->active = 0;
    size_t nw = pairset_words(n);
    if (nw == 0) return;
    size_t need = nw * sizeof(uint64_t);
    if (need > max_bytes) {
        printf("[auto-mem] Pair cache disabled (needs %.2f GB > %.2f GB budget)\n",
               (double)need / 1073741824.0, (double)max_bytes / 1073741824.0);
        fflush(stdout);
        return;
    }
    ps->w = calloc(nw, sizeof(uint64_t));
    if (!ps->w) {
        printf("[auto-mem] Pair cache alloc failed -> disabled\n");
        fflush(stdout);
        return;
    }
    ps->nwords = nw;
    ps->active = 1;
    printf("[auto-mem] Pair cache: %.1f MB (lock-free bitset)\n", (double)need / 1048576.0);
    fflush(stdout);
}

static bool pairset_test_and_insert(PairSet *ps, long int a, long int b) {
    if (!ps->active) return false;
    if (a == b) return true;
    if (a > b) { long int t = a; a = b; b = t; }
    if (a < 0 || b >= ps->n) return true;
    uint64_t ua = (uint64_t)a, ub = (uint64_t)b, un = (uint64_t)ps->n;
    uint64_t idx = ua * un - (ua * (ua + 1ULL)) / 2ULL + (ub - ua - 1ULL);
    size_t w = (size_t)(idx >> 6);
    uint64_t bit = 1ULL << (idx & 63ULL);
    uint64_t old = __atomic_fetch_or(&ps->w[w], bit, __ATOMIC_RELAXED);
    return (old & bit) != 0;
}

static void pairset_free(PairSet *ps) {
    free(ps->w);
    ps->w = NULL;
    ps->nwords = 0;
    ps->active = 0;
}

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int natural_cmp(const char *a, const char *b) {
    while (*a && *b) {
        if (isdigit((unsigned char)*a) && isdigit((unsigned char)*b)) {
            char *ea, *eb;
            unsigned long va = strtoul(a, &ea, 10);
            unsigned long vb = strtoul(b, &eb, 10);
            if (va != vb) return (va < vb) ? -1 : 1;
            a = ea; b = eb;
        } else {
            if (*a != *b) return ((unsigned char)*a < (unsigned char)*b) ? -1 : 1;
            a++; b++;
        }
    }
    return (*a) ? 1 : ((*b) ? -1 : 0);
}

static int cmp_str_natural(const void *a, const void *b) {
    const char *sa = *(const char * const *)a;
    const char *sb = *(const char * const *)b;
    return natural_cmp(sa, sb);
}

static bool find_assign_value(const char *line, char var, double *out) {
    char lower = (var >= 'A' && var <= 'Z') ? (char)(var + 32) : var;
    for (const char *p = line; *p; p++) {
        if (*p != var && *p != lower) continue;
        if (p > line && isalpha((unsigned char)p[-1])) continue;
        const char *q = p + 1;
        while (*q == ' ' || *q == '\t') q++;
        if (*q != '=') continue;
        q++;
        while (*q == ' ' || *q == '\t') q++;
        char *endp = NULL;
        double val = strtod(q, &endp);
        if (endp != q && isfinite(val)) {
            *out = val;
            return true;
        }
    }
    return false;
}

static bool parse_energy_line(const char *line, double *out_e, bool *is_kcal) {
    *is_kcal = false;
    if (!line) return false;
    bool has_kcal = (strstr(line, "kcal") != NULL);
    if (find_assign_value(line, 'E', out_e)) {
        *is_kcal = has_kcal;
        return true;
    }
    double dval;
    if (find_assign_value(line, 'D', &dval)) return false;

    char buf[MC_LINE_SZ];
    strncpy(buf, line, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    char *saveptr = NULL;
    char *tok = strtok_r(buf, " \t\r\n", &saveptr);
    double last_val = 0.0;
    bool found = false;
    while (tok) {
        char *endp = NULL;
        double val = strtod(tok, &endp);
        if (endp != tok && *endp == '\0' && isfinite(val)) {
            last_val = val;
            found = true;
        }
        tok = strtok_r(NULL, " \t\r\n", &saveptr);
    }
    if (found) {
        *out_e = last_val;
        *is_kcal = has_kcal;
        return true;
    }
    return false;
}

static void init_dedup_config(DedupConfig *cfg, const char *mode, bool readtopo) {
    memset(cfg, 0, sizeof(*cfg));
    snprintf(cfg->mode, sizeof(cfg->mode), "%s", mode);
    cfg->readtopo = readtopo;
    cfg->mc_opts = MC_DEFAULT_OPTIONS();
    cfg->mc_opts.STEREO = false;

    if (strcmp(mode, "rmsd-kabsch") == 0) {
        cfg->mc_opts.CONNECT = false;
        strcpy(cfg->mc_opts.ALG, "kabsch");
        cfg->value = 0.125;
    } else if (strcmp(mode, "rmsd-hungarian") == 0) {
        cfg->mc_opts.CONNECT = false;
        strcpy(cfg->mc_opts.ALG, "std");
        cfg->value = 0.125;
    } else if (strcmp(mode, "rmsd-graph") == 0) {
        cfg->mc_opts.CONNECT = true;
        strcpy(cfg->mc_opts.ALG, "std");
        cfg->value = 0.125;
    } else if (strcmp(mode, "angle") == 0) {
        cfg->mc_opts.CONNECT = true;
        strcpy(cfg->mc_opts.ALG, "std");
        cfg->value = 30.0;
        cfg->rmsd_thld = 0.5;
        cfg->maxdist_thld = 1.5;
    } else if (strcmp(mode, "dist") == 0) {
        printf("Dist is not implemented yet. Exit.\n");
        exit(0);
    } else {
        printf("Unknown mode '%s'. Use: rmsd-kabsch, rmsd-hungarian, rmsd-graph, angle\n", mode);
        exit(0);
    }
}

static bool is_duplicate(const DedupConfig *cfg, const McResult *res) {
    if (strcmp(cfg->mode, "rmsd-graph") == 0 ||
        strcmp(cfg->mode, "rmsd-hungarian") == 0 ||
        strcmp(cfg->mode, "rmsd-kabsch") == 0) {
        return (res->rmsd < cfg->value);
    } else if (strcmp(cfg->mode, "angle") == 0) {
        return (res->rmsd < cfg->rmsd_thld &&
                res->dmax < cfg->maxdist_thld &&
                res->amax < cfg->value);
    }
    return false;
}

static double calc_nucrep_flat(const double *xyz, const double *Z_arr, long int natoms) {
    double E = 0.0;
    const double ang_to_bohr = 1.8897261259077824;
    for (long int i = 0; i < natoms - 1; i++) {
        double q1 = Z_arr[i];
        double xi = xyz[3 * i + 0];
        double yi = xyz[3 * i + 1];
        double zi = xyz[3 * i + 2];
        for (long int j = i + 1; j < natoms; j++) {
            double dx = xi - xyz[3 * j + 0];
            double dy = yi - xyz[3 * j + 1];
            double dz = zi - xyz[3 * j + 2];
            double d = sqrt(dx * dx + dy * dy + dz * dz) * ang_to_bohr;
            if (d > 1e-12) E += (q1 * Z_arr[j]) / d;
        }
    }
    return E;
}

static bool load_one_xyz(
    FILE *multi_fp,
    bool single_file,
    const ConfItem *item,
    long int natoms,
    double *xyz_out
) {
    char line[MC_LINE_SZ];
    if (single_file) {
        if (fseek(multi_fp, item->file_offset, SEEK_SET) != 0) return false;
        if (!fgets(line, sizeof(line), multi_fp)) return false;
        if (!fgets(line, sizeof(line), multi_fp)) return false;
        for (long int a = 0; a < natoms; a++) {
            if (!fgets(line, sizeof(line), multi_fp)) return false;
            char el_buf[8];
            if (sscanf(line, "%7s %lf %lf %lf",
                       el_buf, &xyz_out[3*a+0], &xyz_out[3*a+1], &xyz_out[3*a+2]) != 4) {
                return false;
            }
        }
        return true;
    } else {
        FILE *fp = fopen(item->name, "r");
        if (!fp) return false;
        if (!fgets(line, sizeof(line), fp) || !fgets(line, sizeof(line), fp)) {
            fclose(fp);
            return false;
        }
        for (long int a = 0; a < natoms; a++) {
            if (!fgets(line, sizeof(line), fp)) {
                fclose(fp);
                return false;
            }
            char el_buf[8];
            if (sscanf(line, "%7s %lf %lf %lf",
                       el_buf, &xyz_out[3*a+0], &xyz_out[3*a+1], &xyz_out[3*a+2]) != 4) {
                fclose(fp);
                return false;
            }
        }
        fclose(fp);
        return true;
    }
}

static void load_xyz_or_die(FILE *fp, bool single_file, const ConfItem *item,
                            long int natoms, double *out) {
    if (!out) { fprintf(stderr, "FATAL: xyz buffer is NULL\n"); exit(1); }
    if (single_file && !fp) { fprintf(stderr, "FATAL: input stream is NULL\n"); exit(1); }
    if (!load_one_xyz(fp, single_file, item, natoms, out)) {
        fprintf(stderr, "FATAL: cannot read conformer '%s' (offset %ld)\n",
                item->name, item->file_offset);
        exit(1);
    }
}

static ConfItem *read_multixyz(const char *filename, long int *out_nconfs, long int *out_natoms,
                               char (**out_elements)[4], bool *out_mix, bool *out_kcal,
                               bool *out_disk_mode) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { printf("Cannot open file %s\n", filename); exit(0); }
    long int cap = 256, count = 0;
    ConfItem *list = calloc((size_t)cap, sizeof(ConfItem));
    long int expected_natoms = -1;
    char (*elements)[4] = NULL;
    bool any_missing_e = false, detected_kcal = false;
    char line[MC_LINE_SZ];
    
    while (1) {
        long frame_off = ftell(fp);
        if (!fgets(line, sizeof(line), fp)) break;
        char *p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0') continue;
        long int nat = strtol(p, NULL, 10);
        if (nat <= 0) break;
        
        if (expected_natoms < 0) {
            expected_natoms = nat;
            elements = calloc((size_t)nat, sizeof(*elements));
        } else if (nat != expected_natoms) {
            printf("Error: inconsistent NATOMS in %s (%ld vs %ld)\n", filename, nat, expected_natoms);
            exit(0);
        }
        
        if (!fgets(line, sizeof(line), fp)) break;
        double e_val = 0.0;
        bool is_kcal = false;
        bool has_e = parse_energy_line(line, &e_val, &is_kcal);
        if (!has_e) any_missing_e = true;
        else if (is_kcal) detected_kcal = true;
        
        for (long int a = 0; a < nat; a++) {
            if (!fgets(line, sizeof(line), fp)) {
                printf("Unexpected EOF in %s at frame %ld\n", filename, count + 1);
                exit(0);
            }
            char el_buf[8] = "";
            double tx, ty, tz;
            if (sscanf(line, "%7s %lf %lf %lf", el_buf, &tx, &ty, &tz) != 4) {
                printf("Bad coordinate line in %s at frame %ld: %s\n", filename, count + 1, line);
                exit(0);
            }
            if (count == 0) {
                snprintf(elements[a], 4, "%s", el_buf);
            } else if (strncmp(elements[a], el_buf, 3) != 0) {
                printf("Error: atom order mismatch in %s at frame %ld atom %ld (%s vs %s)\n",
                       filename, count + 1, a + 1, el_buf, elements[a]);
                exit(1);
            }
        }
        
        if (count >= cap) { cap *= 2; list = realloc(list, (size_t)cap * sizeof(ConfItem)); }
        memset(&list[count], 0, sizeof(ConfItem));
        snprintf(list[count].name, sizeof(list[count].name), "conf_%ld.xyz", count + 1);
        list[count].orig_idx = count;
        list[count].file_offset = frame_off;
        list[count].xyz = NULL;
        *(list[count].desc + 0) = has_e ? e_val : 0.0;
        list[count].ndesc = has_e ? 1 : 0;
        count++;
    }
    
    size_t avail = get_mem_available_bytes();
    size_t coords_bytes = (count > 0 && expected_natoms > 0)
        ? ((size_t)count * (size_t)expected_natoms * 3ULL * sizeof(double))
        : 0;
    bool use_disk = (avail > 0 && coords_bytes > (avail * 2) / 5);
    const char *env_disk = getenv("MC_DISK_MODE");
    if (env_disk && atoi(env_disk) == 1) use_disk = true;
    
    if (!use_disk && count > 0) {
        rewind(fp);
        for (long int c = 0; c < count; c++) {
            list[c].xyz = malloc((size_t)(3 * expected_natoms) * sizeof(double));
            if (!list[c].xyz || !load_one_xyz(fp, true, &list[c], expected_natoms, list[c].xyz)) {
                for (long int k = 0; k <= c; k++) { free(list[k].xyz); list[k].xyz = NULL; }
                use_disk = true;
                break;
            }
        }
    }
    fclose(fp);
    
    *out_nconfs = count;
    *out_natoms = (expected_natoms > 0) ? expected_natoms : 0;
    *out_elements = elements;
    *out_mix = any_missing_e;
    *out_kcal = detected_kcal;
    *out_disk_mode = use_disk;
    return list;
}

static ConfItem *read_dir_xyz(long int *out_nconfs, long int *out_natoms,
                              char (**out_elements)[4], bool *out_mix, bool *out_kcal,
                              bool *out_disk_mode) {
    DIR *d = opendir(".");
    if (!d) { printf("Cannot open current directory.\n"); exit(0); }
    char **fnames = NULL;
    long int nf = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        size_t len = strlen(ent->d_name);
        if (len > 4 && strcasecmp(ent->d_name + len - 4, ".xyz") == 0) {
            fnames = realloc(fnames, (size_t)(nf + 1) * sizeof(char *));
            fnames[nf] = strdup(ent->d_name);
            nf++;
        }
    }
    closedir(d);
    if (nf == 0) { printf("No .xyz files found in current directory.\n"); exit(0); }
    qsort(fnames, (size_t)nf, sizeof(char *), cmp_str_natural);
    
    ConfItem *list = calloc((size_t)nf, sizeof(ConfItem));
    long int expected_natoms = -1;
    char (*elements)[4] = NULL;
    bool any_missing_e = false, detected_kcal = false;
    
    for (long int i = 0; i < nf; i++) {
        FILE *fp = fopen(fnames[i], "r");
        if (!fp) { printf("Cannot open file %s\n", fnames[i]); exit(1); }
        char line[MC_LINE_SZ];
        long int nat = 0;
        if (fgets(line, sizeof(line), fp)) nat = strtol(line, NULL, 10);
        double e_val = 0.0;
        bool is_kcal = false, has_e = false;
        if (fgets(line, sizeof(line), fp)) {
            has_e = parse_energy_line(line, &e_val, &is_kcal);
        }
        if (!has_e) any_missing_e = true;
        else if (is_kcal) detected_kcal = true;
        
        if (expected_natoms < 0) {
            expected_natoms = nat;
            elements = calloc((size_t)nat, sizeof(*elements));
        } else if (nat != expected_natoms) {
            printf("Error: NAT1 != NAT2 in %s (%ld vs %ld)\n", fnames[i], nat, expected_natoms);
            exit(0);
        }
        
        for (long int a = 0; a < nat; a++) {
            if (!fgets(line, sizeof(line), fp)) {
                printf("Unexpected EOF in %s at atom %ld\n", fnames[i], a + 1);
                exit(0);
            }
            char el_buf[8] = "";
            double tx, ty, tz;
            if (sscanf(line, "%7s %lf %lf %lf", el_buf, &tx, &ty, &tz) != 4) {
                printf("Bad coordinate line in %s at atom %ld\n", fnames[i], a + 1);
                exit(0);
            }
            if (i == 0) {
                snprintf(elements[a], 4, "%s", el_buf);
            } else if (strncmp(elements[a], el_buf, 3) != 0) {
                printf("Error: atom order mismatch in %s at atom %ld (%s vs %s)\n",
                       fnames[i], a + 1, el_buf, elements[a]);
                exit(1);
            }
        }
        fclose(fp);
        strncpy(list[i].name, fnames[i], sizeof(list[i].name) - 1);
        list[i].orig_idx = i;
        list[i].file_offset = 0;
        list[i].xyz = NULL;
        *(list[i].desc + 0) = has_e ? e_val : 0.0;
        list[i].ndesc = has_e ? 1 : 0;
        free(fnames[i]);
    }
    free(fnames);
    
    size_t avail = get_mem_available_bytes();
    size_t coords_bytes = (nf > 0 && expected_natoms > 0)
        ? ((size_t)nf * (size_t)expected_natoms * 3ULL * sizeof(double))
        : 0;
    bool use_disk = (avail > 0 && coords_bytes > (avail * 2) / 5);
    const char *env_disk = getenv("MC_DISK_MODE");
    if (env_disk && atoi(env_disk) == 1) use_disk = true;
    
    if (!use_disk) {
        for (long int i = 0; i < nf; i++) {
            list[i].xyz = malloc((size_t)(3 * expected_natoms) * sizeof(double));
            if (!list[i].xyz || !load_one_xyz(NULL, false, &list[i], expected_natoms, list[i].xyz)) {
                for (long int k = 0; k <= i; k++) { free(list[k].xyz); list[k].xyz = NULL; }
                use_disk = true;
                break;
            }
        }
    }
    
    *out_nconfs = nf;
    *out_natoms = expected_natoms;
    *out_elements = elements;
    *out_mix = any_missing_e;
    *out_kcal = detected_kcal;
    *out_disk_mode = use_disk;
    return list;
}

static int G_SORT_KEY = 0;
static bool G_SORT_REV = false;

static int cmp_conf_ptr(const void *a, const void *b) {
    const ConfItem *ca = *(const ConfItem * const *)a;
    const ConfItem *cb = *(const ConfItem * const *)b;
    double va = ca->desc[G_SORT_KEY];
    double vb = cb->desc[G_SORT_KEY];
    if (va < vb) return G_SORT_REV ? 1 : -1;
    if (va > vb) return G_SORT_REV ? -1 : 1;
    return natural_cmp(ca->name, cb->name);
}

static int cmp_duprec_canon(const void *a, const void *b) {
    const DupRecord *da = (const DupRecord *)a;
    const DupRecord *db = (const DupRecord *)b;
    if (da->idx1 != db->idx1) return (da->idx1 < db->idx1) ? -1 : 1;
    if (da->idx2 != db->idx2) return (da->idx2 < db->idx2) ? -1 : 1;
    return 0;
}

static int cmp_duprec(const void *a, const void *b) {
    const DupRecord *da = (const DupRecord *)a;
    const DupRecord *db = (const DupRecord *)b;
    if (da->dE_kcal < db->dE_kcal) return -1;
    if (da->dE_kcal > db->dE_kcal) return 1;
    return cmp_duprec_canon(a, b);
}

static bool path_exists(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0);
}

static void copy_file_raw(const char *src, const char *dst) {
    FILE *fi = fopen(src, "rb");
    if (!fi) return;
    FILE *fo = fopen(dst, "wb");
    if (!fo) { fclose(fi); return; }
    char buf[MC_COPY_SZ];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), fi)) > 0) fwrite(buf, 1, n, fo);
    fclose(fi);
    fclose(fo);
}

int main(int argc, char **argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 2) {
        printf("Usage:\n  %s <conformers.xyz> <mode> [readtopo]\n  %s <mode> [readtopo]\n", argv[0], argv[0]);
        printf("Modes: rmsd-kabsch, rmsd-hungarian, rmsd-graph, angle\n");
        return 0;
    }
    bool single_file = false;
    char input_file[MC_PATH_SZ] = "";
    char mode[MC_MODE_SZ] = "";
    bool readtopo = false;
    
    if (strstr(argv[1], ".xyz") != NULL || strstr(argv[1], ".XYZ") != NULL) {
        if (argc < 3) { printf("Please specify mode after xyz file.\n"); return 0; }
        single_file = true;
        strncpy(input_file, argv[1], sizeof(input_file) - 1);
        strncpy(mode, argv[2], sizeof(mode) - 1);
    } else {
        single_file = false;
        strncpy(mode, argv[1], sizeof(mode) - 1);
    }
    
    for (int i = 0; mode[i]; i++) mode[i] = (char)tolower((unsigned char)mode[i]);
    if (strcasecmp(argv[argc - 1], "readtopo") == 0) readtopo = true;
    if (readtopo && !path_exists("topo.txt")) {
        printf("Your topo.txt file is missing. Stop.\n");
        return 0;
    }
    
    char base_name[MC_PATH_SZ] = "";
    char out_unique_dir[MC_PATH_SZ] = "";
    char out_unique_xyz[MC_PATH_SZ] = "";
    char out_dup_txt[MC_PATH_SZ] = "";
    char out_dupxyz_txt[MC_PATH_SZ] = "";
    
    if (!single_file) {
        snprintf(out_unique_dir, sizeof(out_unique_dir), "unique-%s", mode);
        snprintf(out_dup_txt, sizeof(out_dup_txt), "dup-%s.txt", mode);
        snprintf(out_dupxyz_txt, sizeof(out_dupxyz_txt), "DUPXYZ-%s.txt", mode);
        if (path_exists(out_unique_dir)) { printf("Folder %s exists. Remove it first.\n", out_unique_dir); return 0; }
        if (path_exists(out_dup_txt)) { printf("File %s exists. Remove it first.\n", out_dup_txt); return 0; }
    } else {
        size_t len = strlen(input_file);
        size_t blen = (len > 4) ? (len - 4) : len;
        if (blen >= sizeof(base_name)) blen = sizeof(base_name) - 1;
        memcpy(base_name, input_file, blen);
        base_name[blen] = '\0';
        snprintf(out_unique_xyz, sizeof(out_unique_xyz), "%s-unique-%s.xyz", base_name, mode);
        snprintf(out_dup_txt, sizeof(out_dup_txt), "%s-dup-%s.txt", base_name, mode);
        snprintf(out_dupxyz_txt, sizeof(out_dupxyz_txt), "%s-DUPXYZ-%s.txt", base_name, mode);
        if (path_exists(out_unique_xyz)) { printf("File %s exists. Remove it first.\n", out_unique_xyz); return 0; }
        if (path_exists(out_dup_txt)) { printf("File %s exists. Remove it first.\n", out_dup_txt); return 0; }
    }
    
    DedupConfig cfg;
    init_dedup_config(&cfg, mode, readtopo);
    long int nconfs = 0, natoms = 0;
    char (*elements)[4] = NULL;
    bool mix_mode = false, is_kcal_units = false, disk_mode = false;
    
    ConfItem *confs = single_file
        ? read_multixyz(input_file, &nconfs, &natoms, &elements, &mix_mode, &is_kcal_units, &disk_mode)
        : read_dir_xyz(&nconfs, &natoms, &elements, &mix_mode, &is_kcal_units, &disk_mode);
        
    if (nconfs <= 0) { printf("No conformers loaded.\n"); return 0; }
    
    size_t avail_ram = get_mem_available_bytes();
    printf("[setup] nconfs=%ld  natoms=%ld  storage=%s  avail_ram=%.2f GB\n",
           nconfs, natoms,
           disk_mode ? "DISK-STREAM (low-RAM)" : "RAM",
           (double)avail_ram / 1073741824.0);
    fflush(stdout);
    
    double *Z_arr = malloc((size_t)natoms * sizeof(double));
    for (long int a = 0; a < natoms; a++) Z_arr[a] = (double)MC_ATOMINFO(elements[a]).Z;
    
    double *xyz_ref0 = malloc((size_t)(3 * natoms) * sizeof(double));
    double *xyz_refN = malloc((size_t)(3 * natoms) * sizeof(double));
    if (!disk_mode) {
        memcpy(xyz_ref0, confs[0].xyz, (size_t)(3 * natoms) * sizeof(double));
        memcpy(xyz_refN, confs[nconfs - 1].xyz, (size_t)(3 * natoms) * sizeof(double));
    } else {
        FILE *ftmp = single_file ? fopen(input_file, "r") : NULL;
        load_xyz_or_die(ftmp, single_file, &confs[0], natoms, xyz_ref0);
        load_xyz_or_die(ftmp, single_file, &confs[nconfs - 1], natoms, xyz_refN);
        if (ftmp) fclose(ftmp);
    }
    
    long int *nn_ref = NULL;
    long int (*nbrs_ref)[23] = NULL;
    if (cfg.mc_opts.CONNECT) {
        McAtom *ref_atoms = calloc((size_t)natoms, sizeof(McAtom));
        for (long int a = 0; a < natoms; a++) {
            strcpy(ref_atoms[a].El, elements[a]);
            ref_atoms[a].xyz[0] = xyz_ref0[3*a+0];
            ref_atoms[a].xyz[1] = xyz_ref0[3*a+1];
            ref_atoms[a].xyz[2] = xyz_ref0[3*a+2];
        }
        if (readtopo) {
            FILE *ft = fopen("topo.txt", "r");
            MC_TXYZREAD(ft, &ref_atoms, natoms);
            fclose(ft);
        } else {
            MC_GETCONNECT(&ref_atoms, natoms);
        }
        MC_SYMMETRYCHECK(&ref_atoms, natoms);
        nn_ref = calloc((size_t)natoms, sizeof(long int));
        nbrs_ref = calloc((size_t)natoms, sizeof(*nbrs_ref));
        for (long int a = 0; a < natoms; a++) {
            nn_ref[a] = ref_atoms[a].nn;
            for (long int k = 0; k < ref_atoms[a].nn; k++) nbrs_ref[a][k] = ref_atoms[a].n[k];
        }
        free(ref_atoms);
    }
    
    McPrecomp *pc = MC_PRECOMP_CREATE(natoms, elements, xyz_ref0, nn_ref, nbrs_ref, &cfg.mc_opts);
    if (mix_mode) {
        for (long int c = 0; c < nconfs; c++) confs[c].ndesc = 0;
    }
    
    McOptions rmsd_opts = MC_DEFAULT_OPTIONS();
    rmsd_opts.CONNECT = false; rmsd_opts.STEREO = false;
    strcpy(rmsd_opts.ALG, "std");
    
    #pragma omp parallel
    {
        FILE *t_fp = (disk_mode && single_file) ? fopen(input_file, "r") : NULL;
        double *t_xyz = disk_mode ? malloc((size_t)(3 * natoms) * sizeof(double)) : NULL;
        #pragma omp for schedule(dynamic, 16)
        for (long int c = 0; c < nconfs; c++) {
            const double *c_xyz = confs[c].xyz;
            if (disk_mode) {
                load_xyz_or_die(t_fp, single_file, &confs[c], natoms, t_xyz);
                c_xyz = t_xyz;
            }
            int idx = confs[c].ndesc;
            confs[c].desc[idx] = calc_nucrep_flat(c_xyz, Z_arr, natoms);
            McResult r1 = MC_COMPARE_FLAT(natoms, elements, NULL, NULL, c_xyz, xyz_ref0, &rmsd_opts);
            confs[c].desc[idx + 1] = (r1.status == 0) ? r1.rmsd : 0.0;
            if (mix_mode) {
                McResult r2 = MC_COMPARE_FLAT(natoms, elements, NULL, NULL, c_xyz, xyz_refN, &rmsd_opts);
                confs[c].desc[idx + 2] = (r2.status == 0) ? r2.rmsd : 0.0;
            }
        }
        if (t_fp) fclose(t_fp);
        free(t_xyz);
    }
    
    for (long int c = 0; c < nconfs; c++) confs[c].ndesc += mix_mode ? 3 : 2;
    int t_desc = confs[0].ndesc;
    free(xyz_ref0);
    free(xyz_refN);
    
    ConfItem **active = malloc((size_t)nconfs * sizeof(ConfItem *));
    long int n_active = nconfs;
    for (long int c = 0; c < nconfs; c++) active[c] = &confs[c];
    
    G_SORT_KEY = 0; G_SORT_REV = false;
    qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);
    
    avail_ram = get_mem_available_bytes();
    size_t cache_budget = (avail_ram > 0) ? (avail_ram / 2) : (2ULL * 1024 * 1024 * 1024);
    size_t max_cap = 8ULL * 1024 * 1024 * 1024;
    if (cache_budget > max_cap) cache_budget = max_cap;
    const char *ev = getenv("MC_CACHE_GB");
    if (ev && atof(ev) > 0) cache_budget = (size_t)(atof(ev) * 1073741824.0);
    
    PairSet tested_pairs;
    pairset_init(&tested_pairs, nconfs, cache_budget);
    
    DupRecord *dups = NULL;
    long int ndups = 0, dups_cap = 0, total_tested = 0;
    int chunk_base = 24, chunk_sz = chunk_base, iter_max = 3000, perm_max = 8, pm = 0;
    unsigned int rng_seed = 42;
    long int pair_budget = 200000000L;
    double time_budget = 6.0 * 3600.0;
    
    ev = getenv("MC_PAIR_BUDGET"); if (ev && atol(ev) > 0) pair_budget = atol(ev);
    ev = getenv("MC_MAX_HOURS");   if (ev && atof(ev) > 0) time_budget = atof(ev) * 3600.0;
    ev = getenv("MC_MAX_ITER");    if (ev && atoi(ev) > 0) iter_max = atoi(ev);
    
    double t_start = now_sec();
    printf("Deduplication started (%ld conformers)...\n", nconfs);
    fflush(stdout);
    
    for (int it = 0; it < iter_max; it++) {
        long int prev_dups = ndups;
        long int n_chunks = (n_active + chunk_sz - 1) / chunk_sz;
        DupRecord *iter_dups = NULL;
        long int n_iter_dups = 0, iter_cap = 0;
        
        #pragma omp parallel
        {
            FILE *t_fp = (disk_mode && single_file) ? fopen(input_file, "r") : NULL;
            double *chunk_xyz = NULL;
            long int chunk_xyz_cap = 0;
            DupRecord *loc_dups = NULL;
            long int loc_ndups = 0, loc_cap = 0, loc_tested = 0;
            
            #pragma omp for schedule(dynamic, 1)
            for (long int ch = 0; ch < n_chunks; ch++) {
                long int start = ch * chunk_sz;
                long int end = start + chunk_sz;
                if (end > n_active) end = n_active;
                long int len = end - start;
                if (len <= 1) continue;
                
                if (disk_mode) {
                    if (len > chunk_xyz_cap) {
                        chunk_xyz_cap = len;
                        chunk_xyz = realloc(chunk_xyz, (size_t)(chunk_xyz_cap * 3 * natoms) * sizeof(double));
                    }
                    for (long int k = 0; k < len; k++) {
                        load_xyz_or_die(t_fp, single_file, active[start + k], natoms,
                                     chunk_xyz + (size_t)(k * 3 * natoms));
                    }
                }
                
                bool *loc_removed = calloc((size_t)len, sizeof(bool));
                for (long int i = 0; i < len - 1; i++) {
                    if (loc_removed[i]) continue;
                    ConfItem *ci = active[start + i];
                    const double *xyz_i = disk_mode ? (chunk_xyz + (size_t)(i * 3 * natoms)) : ci->xyz;
                    for (long int m = i + 1; m < len; m++) {
                        if (loc_removed[m]) continue;
                        ConfItem *cm = active[start + m];
                        const double *xyz_m = disk_mode ? (chunk_xyz + (size_t)(m * 3 * natoms)) : cm->xyz;
                        ConfItem *keep = ci, *cand = cm;
                        const double *xyz_keep = xyz_i, *xyz_cand = xyz_m;
                        
                        if (*(cm->desc + 0) < *(ci->desc + 0)) {
                            keep = cm; cand = ci; xyz_keep = xyz_m; xyz_cand = xyz_i;
                        } else if (*(cm->desc + 0) == *(ci->desc + 0)) {
                            if (natural_cmp(cm->name, ci->name) < 0) {
                                keep = cm; cand = ci; xyz_keep = xyz_m; xyz_cand = xyz_i;
                            }
                        }
                        
                        long int id_a = (keep->orig_idx < cand->orig_idx) ? keep->orig_idx : cand->orig_idx;
                        long int id_b = (keep->orig_idx < cand->orig_idx) ? cand->orig_idx : keep->orig_idx;
                        
                        if (pairset_test_and_insert(&tested_pairs, id_a, id_b)) continue;
                        loc_tested++;
                        
                        McResult res = MC_COMPARE_PRECOMP(pc, xyz_keep, xyz_cand, &cfg.mc_opts);
                        if (res.status == 0 && is_duplicate(&cfg, &res)) {
                            double dE = fabs(*(cand->desc + 0) - *(keep->desc + 0));
                            double dE_kcal = is_kcal_units ? dE : (dE * 627.51);
                            if (loc_ndups >= loc_cap) {
                                loc_cap = (loc_cap == 0) ? 16 : loc_cap * 2;
                                loc_dups = realloc(loc_dups, (size_t)loc_cap * sizeof(DupRecord));
                            }
                            loc_dups[loc_ndups].idx1 = (int)keep->orig_idx;
                            loc_dups[loc_ndups].idx2 = (int)cand->orig_idx;
                            loc_dups[loc_ndups].out[0] = res.rmsd;
                            loc_dups[loc_ndups].out[1] = res.dmax;
                            loc_dups[loc_ndups].out[2] = res.amax;
                            loc_dups[loc_ndups].dE_kcal = dE_kcal;
                            loc_ndups++;
                            
                            if (cand == cm) loc_removed[m] = true;
                            else { loc_removed[i] = true; break; }
                        }
                    }
                }
                free(loc_removed);
            }
            
            #pragma omp critical(merge_dups)
            {
                total_tested += loc_tested;
                for (long int k = 0; k < loc_ndups; k++) {
                    if (n_iter_dups >= iter_cap) {
                        iter_cap = (iter_cap == 0) ? 64 : iter_cap * 2;
                        iter_dups = realloc(iter_dups, (size_t)iter_cap * sizeof(DupRecord));
                    }
                    iter_dups[n_iter_dups++] = loc_dups[k];
                }
            }
            free(loc_dups);
            free(chunk_xyz);
            if (t_fp) fclose(t_fp);
        }
        
        if (n_iter_dups > 1) qsort(iter_dups, (size_t)n_iter_dups, sizeof(DupRecord), cmp_duprec_canon);
        for (long int k = 0; k < n_iter_dups; k++) {
            int keep_orig = iter_dups[k].idx1;
            int dup_orig  = iter_dups[k].idx2;
            if (!confs[keep_orig].removed && !confs[dup_orig].removed) {
                confs[dup_orig].removed = true;
                if (ndups >= dups_cap) {
                    dups_cap = (dups_cap == 0) ? 64 : dups_cap * 2;
                    dups = realloc(dups, (size_t)dups_cap * sizeof(DupRecord));
                }
                dups[ndups++] = iter_dups[k];
            }
        }
        free(iter_dups);
        
        double elapsed = now_sec() - t_start;
        if (ndups > prev_dups || it == 0 || (it % 50 == 0)) {
            printf("Iter %4d (perm %d, chunk %3d): tested=%ld | DUP=%ld | active=%ld | %.1f min\n",
                   it, pm, chunk_sz, total_tested, ndups, n_active, elapsed / 60.0);
            fflush(stdout);
        }
        
        if (total_tested >= pair_budget || elapsed > time_budget) {
            printf("[stop] Budget reached (tested=%ld, %.1f min)\n", total_tested, elapsed / 60.0);
            fflush(stdout);
            break;
        }
        
        if (ndups == prev_dups && pm >= perm_max) {
            break;
        } else if (ndups == prev_dups && pm < perm_max) {
            int coin = 0;
            if (t_desc > 1) coin = (int)(rand_r(&rng_seed) % (unsigned int)t_desc);
            int step = (128 - chunk_base) / perm_max + coin;
            if (step < 1) step = 1;
            chunk_sz += step;
            G_SORT_KEY = coin;
            G_SORT_REV = (chunk_sz % 2 == 0);
            qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);
            pm++;
        } else {
            long int w = 0;
            for (long int r = 0; r < n_active; r++) {
                if (!active[r]->removed) active[w++] = active[r];
            }
            n_active = w;
            pm = 0;
            chunk_sz = chunk_base;
        }
    }
    
    G_SORT_KEY = 0; G_SORT_REV = false;
    qsort(active, (size_t)n_active, sizeof(ConfItem *), cmp_conf_ptr);
    
    bool do_exact = false;
    if (disk_mode) {
        printf("[exact] SKIPPED: low-RAM disk streaming mode active\n");
    } else if (!mix_mode && n_active > 1) {
        double e_span_kcal = fabs(*(active[n_active - 1]->desc + 0) - *(active[0]->desc + 0))
                           * (is_kcal_units ? 1.0 : 627.51);
        if (e_span_kcal > 2.0 && n_active <= 1000) {
            do_exact = true;
            printf("[exact] Running exact pass on %ld structures (span=%.2f kcal/mol)...\n",
                   n_active, e_span_kcal);
        } else if (n_active > 1000) {
            printf("[exact] SKIPPED: n_active=%ld > 1000 (chunk phase covers large ensembles; MC_FORCE_EXACT=1 to override)\n", n_active);
        } else {
            printf("[exact] SKIPPED: energy span %.2f kcal/mol <= 2.0\n", e_span_kcal);
        }
    } else if (mix_mode) {
        printf("[exact] SKIPPED: no energies parsed (mix_mode)\n");
    }
    
    const char *force_exact = getenv("MC_FORCE_EXACT");
    if (!disk_mode && force_exact && atoi(force_exact) == 1) {
        do_exact = true;
        printf("[exact] FORCED via MC_FORCE_EXACT=1\n");
    }
    fflush(stdout);
    
    if (n_active > 1 && do_exact) {
        DupRecord *exact_dups = NULL;
        long int n_exact = 0, exact_cap = 0, exact_tested = 0;
        long int done_i = 0;
        int exact_timed_out = 0;
        long int step_i = (n_active > 10) ? ((n_active - 1) / 10) : 1;
        
        #pragma omp parallel
        {
            DupRecord *loc_dups = NULL;
            long int loc_ndups = 0, loc_cap = 0, loc_tested = 0;
            
            #pragma omp for schedule(dynamic, 4)
            for (long int i = 0; i < n_active - 1; i++) {
                long int cur_done = __atomic_add_fetch(&done_i, 1, __ATOMIC_RELAXED);
                if (cur_done % step_i == 0) {
                    printf("  [exact] %ld / %ld (%.0f%%) | %.1f min\n",
                           cur_done, n_active - 1,
                           100.0 * (double)cur_done / (double)(n_active - 1),
                           (now_sec() - t_start) / 60.0);
                    fflush(stdout);
                }
                if ((now_sec() - t_start) > time_budget) { __atomic_store_n(&exact_timed_out, 1, __ATOMIC_RELAXED); continue; }
                
                ConfItem *ci = active[i];
                double d0_i = *(ci->desc + 0);
                for (long int m = i + 1; m < n_active; m++) {
                    ConfItem *cm = active[m];
                    double d0_m = *(cm->desc + 0);
                    if (!mix_mode) {
                        double dE_kcal = fabs(d0_m - d0_i) * (is_kcal_units ? 1.0 : 627.51);
                        if (dE_kcal > 2.0) break;
                    }
                    long int id_a = (ci->orig_idx < cm->orig_idx) ? ci->orig_idx : cm->orig_idx;
                    long int id_b = (ci->orig_idx < cm->orig_idx) ? cm->orig_idx : ci->orig_idx;
                    
                    if (pairset_test_and_insert(&tested_pairs, id_a, id_b)) continue;
                    loc_tested++;
                    
                    McResult res = MC_COMPARE_PRECOMP(pc, ci->xyz, cm->xyz, &cfg.mc_opts);
                    if (res.status == 0 && is_duplicate(&cfg, &res)) {
                        double dE = fabs(d0_m - d0_i);
                        double dE_kcal = mix_mode ? 0.0 : (is_kcal_units ? dE : (dE * 627.51));
                        if (loc_ndups >= loc_cap) {
                            loc_cap = (loc_cap == 0) ? 16 : loc_cap * 2;
                            loc_dups = realloc(loc_dups, (size_t)loc_cap * sizeof(DupRecord));
                        }
                        loc_dups[loc_ndups].idx1 = (int)ci->orig_idx;
                        loc_dups[loc_ndups].idx2 = (int)cm->orig_idx;
                        loc_dups[loc_ndups].out[0] = res.rmsd;
                        loc_dups[loc_ndups].out[1] = res.dmax;
                        loc_dups[loc_ndups].out[2] = res.amax;
                        loc_dups[loc_ndups].dE_kcal = dE_kcal;
                        loc_ndups++;
                    }
                }
            }
            
            #pragma omp critical(merge_exact)
            {
                exact_tested += loc_tested;
                for (long int k = 0; k < loc_ndups; k++) {
                    if (n_exact >= exact_cap) {
                        exact_cap = (exact_cap == 0) ? 64 : exact_cap * 2;
                        exact_dups = realloc(exact_dups, (size_t)exact_cap * sizeof(DupRecord));
                    }
                    exact_dups[n_exact++] = loc_dups[k];
                }
            }
            free(loc_dups);
        }
        
        if (n_exact > 1) qsort(exact_dups, (size_t)n_exact, sizeof(DupRecord), cmp_duprec_canon);
        for (long int k = 0; k < n_exact; k++) {
            int keep_orig = exact_dups[k].idx1;
            int dup_orig  = exact_dups[k].idx2;
            if (!confs[keep_orig].removed && !confs[dup_orig].removed) {
                confs[dup_orig].removed = true;
                if (ndups >= dups_cap) {
                    dups_cap = (dups_cap == 0) ? 64 : dups_cap * 2;
                    dups = realloc(dups, (size_t)dups_cap * sizeof(DupRecord));
                }
                dups[ndups++] = exact_dups[k];
            }
        }
        free(exact_dups);
        
        long int w = 0;
        for (long int r = 0; r < n_active; r++) {
            if (!active[r]->removed) active[w++] = active[r];
        }
        n_active = w;
        total_tested += exact_tested;
        printf("[exact] %sTested %ld extra pairs, total DUP = %ld\n", exact_timed_out ? "PARTIAL (time budget hit): " : "", exact_tested, ndups);
        fflush(stdout);
    }
    
    printf("UNIQUE = %ld  (DUP = %ld, tested = %ld, elapsed = %.1f min)\n",
           n_active, ndups, total_tested, (now_sec() - t_start) / 60.0);
    fflush(stdout);
    
    if (ndups > 1) qsort(dups, (size_t)ndups, sizeof(DupRecord), cmp_duprec);
    
    FILE *fd = fopen(out_dup_txt, "w");
    FILE *fdx = fopen(out_dupxyz_txt, "w");
    if (fd && fdx) {
        fprintf(fd, "DUP = %ld\n", ndups);
        fprintf(fdx, "DUP = %ld (2nd is the DUP)\n", ndups);
        for (long int i = 0; i < ndups; i++) {
            const ConfItem *c1 = &confs[dups[i].idx1];
            const ConfItem *c2 = &confs[dups[i].idx2];
            fprintf(fdx, "%s %s\n", c1->name, c2->name);
            fprintf(fd, "%s %s %.3f %.3f %.3f %.2f\n", c1->name, c2->name,
                    dups[i].out[0], dups[i].out[1], dups[i].out[2], dups[i].dE_kcal);
        }
    }
    if (fd) fclose(fd);
    if (fdx) fclose(fdx);
    
    if (!single_file) {
        mkdir(out_unique_dir, 0755);
        for (long int i = 0; i < n_active; i++) {
            char dst[MC_PATH_SZ];
            snprintf(dst, sizeof(dst), "%s/%s", out_unique_dir, active[i]->name);
            copy_file_raw(active[i]->name, dst);
        }
    } else {
        FILE *fo = fopen(out_unique_xyz, "w");
        FILE *fi = disk_mode ? fopen(input_file, "r") : NULL;
        double *tmp_xyz = disk_mode ? malloc((size_t)(3 * natoms) * sizeof(double)) : NULL;
        if (fo) {
            for (long int i = 0; i < n_active; i++) {
                const ConfItem *ci = active[i];
                const double *c_xyz = ci->xyz;
                if (disk_mode) {
                    load_xyz_or_die(fi, true, ci, natoms, tmp_xyz);
                    c_xyz = tmp_xyz;
                }
                fprintf(fo, "%ld\n", natoms);
                if (mix_mode) fprintf(fo, "%s\n", ci->name);
                else if (is_kcal_units) fprintf(fo, "Conf %ld  E = %.6f kcal/mol\n", i + 1, *(ci->desc + 0));
                else fprintf(fo, "%.10f\n", *(ci->desc + 0));
                for (long int a = 0; a < natoms; a++) {
                    fprintf(fo, "%3s %14.8f %14.8f %14.8f\n",
                            elements[a], c_xyz[3*a+0], c_xyz[3*a+1], c_xyz[3*a+2]);
                }
            }
            fclose(fo);
        }
        if (fi) fclose(fi);
        free(tmp_xyz);
    }
    
    MC_PRECOMP_FREE(pc);
    for (long int c = 0; c < nconfs; c++) free(confs[c].xyz);
    free(confs);
    free(elements);
    free(Z_arr);
    free(nn_ref);
    free(nbrs_ref);
    free(active);
    pairset_free(&tested_pairs);
    free(dups);
    
    return 0;
}
