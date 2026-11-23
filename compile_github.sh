#!/bin/bash
set -e

CC=${CC:-gcc}
# -march=x86-64 гарантирует запуск на любом 64-битном проце с 2003 года
CFLAGS="-std=c11 -O3 -march=x86-64 -mtune=generic -fopenmp -Wall -Wextra"
LDFLAGS="-static -O3 -fopenmp"

OUT_DIR="release"
mkdir -p "$OUT_DIR"

$CC -c mc_atominfo.c    $CFLAGS
$CC -c mc_inertia.c     $CFLAGS
$CC -c mc_kabshangle.c  $CFLAGS
$CC -c mc_getconnect.c  $CFLAGS
$CC -c mc_stdfunc.c     $CFLAGS
$CC -c mc_hungarian.c   $CFLAGS
$CC -c mc_renumb_hung.c $CFLAGS
$CC -c mc_core.c        $CFLAGS
$CC -c mc_io.c          $CFLAGS
$CC -c main.c           $CFLAGS
$CC -c moldedupl.c      $CFLAGS

OBJS="mc_atominfo.o mc_inertia.o mc_kabshangle.o mc_getconnect.o mc_stdfunc.o mc_hungarian.o mc_renumb_hung.o mc_core.o mc_io.o"

$CC -o "$OUT_DIR/molcomp"   $OBJS main.o      $LDFLAGS -lm
$CC -o "$OUT_DIR/moldedupl" $OBJS moldedupl.o $LDFLAGS -lm

rm -f *.o

# Уменьшаем размер бинарников (убираем отладочные символы)
strip "$OUT_DIR/molcomp" "$OUT_DIR/moldedupl"

echo "Built portable static binaries in $OUT_DIR/:"
ls -lh "$OUT_DIR/molcomp" "$OUT_DIR/moldedupl"
ldd "$OUT_DIR/moldedupl" || true
