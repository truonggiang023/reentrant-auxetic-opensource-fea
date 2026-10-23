#!/bin/bash
# Debian amd64 cloud build, with signed apt metadata and artifact verification.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
root=${FEA_BUILD_ROOT:-/workspace}
mkdir -p "$root/apt/sources" "$root/apt/lists/partial" "$root/apt/cache/archives/partial" "$root/compiler" "$root/ccx-runtime" "$root/ccx-user" "$root/bin"
cat > "$root/apt/sources/debian.list" <<'EOF'
deb [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian bookworm main
deb-src [signed-by=/usr/share/keyrings/debian-archive-keyring.gpg] https://deb.debian.org/debian bookworm main
EOF
opts=(-o "Dir::Etc::sourcelist=$root/apt/sources/debian.list" -o Dir::Etc::sourceparts=- -o "Dir::State::lists=$root/apt/lists" -o "Dir::Cache=$root/apt/cache" -o "Acquire::https::Proxy=${HTTPS_PROXY:?}")
cd "$root/apt"
apt-get "${opts[@]}" update
apt-get "${opts[@]}" download gfortran-12 libgfortran-12-dev libspooles2.2 libspooles-dev
for p in gfortran-12 libgfortran-12-dev; do
    for deb in "$p"_*.deb; do dpkg-deb -x "$deb" "$root/compiler"; done
done
for p in libspooles2.2 libspooles-dev; do
    for deb in "$p"_*.deb; do dpkg-deb -x "$deb" "$root/ccx-runtime"; done
done
apt-get "${opts[@]}" source --download-only calculix-ccx=2.20-1
src="$root/apt/calculix-ccx-2.20"
if [ ! -d "$src" ]; then dpkg-source -x calculix-ccx_2.20-1.dsc "$src"; fi
python "$repo/native/patch_hourglass.py" "$src/ccx_2.20/src"
python "$repo/native/patch_timepoints.py" "$src/ccx_2.20/src"
if ! cmp -s "$repo/native/timepoints.h" "$src/ccx_2.20/src/reentrant_timepoints.h"; then
    cp "$repo/native/timepoints.h" "$src/ccx_2.20/src/reentrant_timepoints.h"
    touch "$src/ccx_2.20/src/checkconvergence.c"
fi
cd "$src/ccx_2.20/src"
# GCC 14 diagnoses legacy C declarations as errors. Preserve the upstream
# source and use warning compatibility flags; Fortran uses its matching driver.
make -j"${FEA_BUILD_JOBS:-4}" CC=gcc FC="$root/compiler/usr/bin/gfortran-12" \
    "CFLAGS=-O2 -Wno-error=implicit-function-declaration -Wno-error=incompatible-pointer-types -Wno-error=int-conversion -I$root/ccx-runtime/usr/include/spooles -DARCH=Linux -DSPOOLES -DARPACK -DMATRIXSTORAGE -DCALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT" \
    'FFLAGS=-O2 -fallow-argument-mismatch' \
    "LIBS=-L$root/ccx-runtime/usr/lib/x86_64-linux-gnu -lspooles -lpthread /lib/x86_64-linux-gnu/libarpack.so.2 /lib/x86_64-linux-gnu/liblapack.so.3 /lib/x86_64-linux-gnu/libblas.so.3 -lm -lc -ldl" \
    'LDFLAGS=-B/usr/lib/gcc/x86_64-linux-gnu/14'
# Build then rename: never truncate a library mapped by an active solver.
gcc -O2 -Wall -Wextra -fPIC -shared "$repo/native/ductile.c" \
    /lib/x86_64-linux-gnu/liblapack.so.3 -lm -o "$root/ccx-user/libREENTRANT.next.so"
mv "$root/ccx-user/libREENTRANT.next.so" "$root/ccx-user/libREENTRANT.so"
cat > "$root/bin/ccx-damage" <<EOF
#!/bin/sh
export LD_LIBRARY_PATH="$root/ccx-user:$root/ccx-runtime/usr/lib/x86_64-linux-gnu\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
exec "$src/ccx_2.20/src/ccx_2.20" "\$@"
EOF
chmod +x "$root/bin/ccx-damage"
sha256sum "$src/ccx_2.20/src/ccx_2.20" "$root/ccx-user/libREENTRANT.so" "$repo/native/ductile.c" > "$root/ccx-user/build-sha256.txt"
OMP_NUM_THREADS=2 CCX_EXECUTABLE="$root/bin/ccx-damage" CCX_DAMAGE_EXECUTABLE="$root/bin/ccx-damage" \
    "$root/fea-venv/bin/python" -m pytest "$repo/tests/test_native_damage.py" -q
