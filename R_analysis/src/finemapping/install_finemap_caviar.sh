#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Install FINEMAP + CAVIAR for use with the finemapr R package.
#
# Adapted from finemapr's own instructions:
#   https://github.com/variani/finemapr/blob/master/misc/install-finemaping-tools.md
# with two deliberate deviations:
#
#   1. Everything is installed under <project>/tools/ (gitignored) instead of
#      ~/.local/apps/, so the third-party binaries live with the project and
#      never enter version control.
#   2. FINEMAP is pinned to v1.1. This is NOT arbitrary: finemapr's
#      run_finemap() writes a v1.1-style master file
#          z;ld;snp;config;log;n-ind
#      and a 2-column "<snp> <zscore>" .z file. FINEMAP >=1.2 renamed the
#      master fields (n-ind -> n_samples, added a `cred` column) and requires a
#      7-column .z file with a header (rsid chromosome position allele1
#      allele2 maf beta se). Installing v1.4.x would break the wrapper.
#
# This script is designed to run inside WSL (Ubuntu). FINEMAP ships Linux and
# macOS builds only, and CAVIAR is a C++ program that expects a POSIX
# toolchain, so neither can be built natively on Windows.
#
# Build prerequisites (already present on this machine; install with apt if
# you are setting up elsewhere):
#   build-essential wget git libgsl-dev liblapack-dev libblas-dev
#
# Usage (from Windows):
#   wsl -d Ubuntu -- bash /mnt/h/.../src/finemapping/install_finemap_caviar.sh
# ---------------------------------------------------------------------------
set -euo pipefail

# Resolve the project root from this script's location (src/finemapping/..).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

TOOLS_DIR="${PROJECT_ROOT}/tools"
BIN_DIR="${TOOLS_DIR}/bin"
DL_DIR="${TOOLS_DIR}/downloads"

FINEMAP_VER="finemap_v1.1_x86_64"
FINEMAP_URL="http://www.christianbenner.com/${FINEMAP_VER}.tgz"
CAVIAR_REPO="https://github.com/fhormoz/caviar.git"

mkdir -p "${BIN_DIR}" "${DL_DIR}"

echo "=============================================================="
echo "Project root : ${PROJECT_ROOT}"
echo "Install dir  : ${TOOLS_DIR}"
echo "=============================================================="

# ---------------------------------------------------------------------------
# 1. FINEMAP v1.1 (precompiled binary)
# ---------------------------------------------------------------------------
echo
echo "--- [1/2] FINEMAP ${FINEMAP_VER} ---"

if [[ -x "${BIN_DIR}/finemap" ]]; then
  echo "already installed at ${BIN_DIR}/finemap - skipping"
else
  cd "${DL_DIR}"
  if [[ ! -f "${FINEMAP_VER}.tgz" ]]; then
    echo "downloading ${FINEMAP_URL}"
    wget --no-verbose --timeout=60 --tries=3 "${FINEMAP_URL}"
  else
    echo "tarball already downloaded - skipping fetch"
  fi

  tar -xzf "${FINEMAP_VER}.tgz"
  cp "${FINEMAP_VER}/${FINEMAP_VER}" "${BIN_DIR}/finemap"
  chmod +x "${BIN_DIR}/finemap"

  # The bundled example/ (region1.z, region1.ld, ...) is what the smoke test
  # and finemapr's own examples use; keep it next to the binary.
  if [[ -d "${FINEMAP_VER}/example" ]]; then
    rm -rf "${TOOLS_DIR}/finemap_example"
    cp -r "${FINEMAP_VER}/example" "${TOOLS_DIR}/finemap_example"
  fi
  echo "installed -> ${BIN_DIR}/finemap"
fi

# ---------------------------------------------------------------------------
# 1b. GSL (dependency of CAVIAR), built from source into tools/
# ---------------------------------------------------------------------------
# This machine has only the *runtime* libs (libblas3, liblapack3) and no GSL at
# all, and `sudo` is password-protected - so we cannot `apt install libgsl-dev`
# non-interactively. GSL is therefore built from source into tools/gsl, which
# needs no root and keeps the dependency inside the gitignored folder.
#
# If you are on a machine where you *do* have root, this whole step is
# replaceable by:
#     sudo apt-get install -y libgsl-dev liblapack-dev libblas-dev

GSL_VER="2.7.1"
GSL_PREFIX="${TOOLS_DIR}/gsl"
LIB_DIR="${TOOLS_DIR}/lib"

echo
echo "--- [1b/2] GSL ${GSL_VER} (CAVIAR dependency) ---"

if [[ -f "${GSL_PREFIX}/include/gsl/gsl_linalg.h" ]]; then
  echo "already built at ${GSL_PREFIX} - skipping"
else
  cd "${DL_DIR}"
  if [[ ! -f "gsl-${GSL_VER}.tar.gz" ]]; then
    echo "downloading GSL ${GSL_VER}"
    wget --no-verbose --timeout=60 --tries=3 \
      "https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VER}.tar.gz"
  fi
  tar -xzf "gsl-${GSL_VER}.tar.gz"
  cd "gsl-${GSL_VER}"
  echo "configuring GSL (this takes a few minutes)..."
  ./configure --prefix="${GSL_PREFIX}" --enable-static --disable-shared \
    > "${DL_DIR}/gsl_configure.log" 2>&1
  echo "building GSL..."
  make -j"$(nproc)" > "${DL_DIR}/gsl_make.log" 2>&1
  make install > "${DL_DIR}/gsl_install.log" 2>&1
  echo "installed -> ${GSL_PREFIX}"
fi

# CAVIAR's Makefile links -llapack -lblas, which needs the unversioned .so
# names supplied by the -dev packages. Only libblas.so.3 / liblapack.so.3 exist
# here, so provide the missing symlinks locally and link against those.
mkdir -p "${LIB_DIR}"
for lib in blas lapack; do
  if [[ ! -e "${LIB_DIR}/lib${lib}.so" ]]; then
    src="$(ls /usr/lib/x86_64-linux-gnu/lib${lib}.so.3* 2>/dev/null | head -1)"
    if [[ -n "${src}" ]]; then
      ln -sf "${src}" "${LIB_DIR}/lib${lib}.so"
      echo "symlinked lib${lib}.so -> ${src}"
    else
      echo "WARNING: no lib${lib}.so.3 found; CAVIAR link step may fail"
    fi
  fi
done

# ---------------------------------------------------------------------------
# 2. CAVIAR (compiled from source)
# ---------------------------------------------------------------------------
echo
echo "--- [2/2] CAVIAR ---"

if [[ -x "${BIN_DIR}/CAVIAR" ]]; then
  echo "already installed at ${BIN_DIR}/CAVIAR - skipping"
else
  CAVIAR_SRC="${TOOLS_DIR}/caviar-src"
  if [[ ! -d "${CAVIAR_SRC}/.git" ]]; then
    echo "cloning ${CAVIAR_REPO}"
    git clone --depth 1 "${CAVIAR_REPO}" "${CAVIAR_SRC}"
  else
    echo "source already cloned - skipping"
  fi

  cd "${CAVIAR_SRC}/CAVIAR-C++"
  echo "compiling (make)..."

  # Invoke g++ directly rather than `make`: the shipped Makefile hardcodes its
  # include/link flags, and we need to inject the locally-built GSL and the
  # local lib*.so symlinks. Source list and flags mirror Makefile:17, with the
  # bundled header-only armadillo used via -DARMA_DONT_USE_WRAPPER.
  g++ caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp \
    -I "$(pwd)/armadillo/include/" \
    -I "${GSL_PREFIX}/include" \
    -DARMA_DONT_USE_WRAPPER \
    -L "${GSL_PREFIX}/lib" -L "${LIB_DIR}" \
    -llapack -lblas -lgsl -lgslcblas \
    -o CAVIAR

  cp CAVIAR "${BIN_DIR}/CAVIAR"
  chmod +x "${BIN_DIR}/CAVIAR"
  echo "installed -> ${BIN_DIR}/CAVIAR"
fi

# ---------------------------------------------------------------------------
# 2c. finemapr R package source (cloned + patched, installed separately)
# ---------------------------------------------------------------------------
# The R package is installed from this local clone by
# src/finemapping/finemapr_setup.R. One patch is applied here.
#
# PATCH (R 4.0 compatibility): run_caviar() guards its inputs with
#     stopifnot(class(ld) == "matrix")
#     stopifnot(class(dir_run) == "character")
# Since R 4.0.0, class(<matrix>) returns c("matrix", "array"), so the first
# comparison yields c(TRUE, FALSE) and stopifnot() aborts - run_caviar() is
# simply unusable on modern R. Rewriting both to is.matrix()/is.character()
# preserves the intended check. run_finemap() already uses is.matrix() and
# needs no patch. finemapr was last updated in 2018, predating this change.

FINEMAPR_SRC="${TOOLS_DIR}/finemapr-src"

echo
echo "--- [2c] finemapr R source ---"

if [[ ! -d "${FINEMAPR_SRC}/.git" ]]; then
  echo "cloning finemapr"
  git clone --depth 1 https://github.com/variani/finemapr.git "${FINEMAPR_SRC}"
else
  echo "source already cloned - skipping"
fi

if grep -q 'stopifnot(class(ld) == "matrix")' "${FINEMAPR_SRC}/R/caviar.R"; then
  sed -i \
    -e 's/stopifnot(class(ld) == "matrix")/stopifnot(is.matrix(ld))/' \
    -e 's/stopifnot(class(dir_run) == "character")/stopifnot(is.character(dir_run))/' \
    "${FINEMAPR_SRC}/R/caviar.R"
  echo "patched R/caviar.R for R >= 4.0 class() semantics"
else
  echo "R/caviar.R already patched - skipping"
fi

# ---------------------------------------------------------------------------
# 2b. Windows -> WSL shims
# ---------------------------------------------------------------------------
# R runs natively on Windows here (R-4.5.1) but finemap/CAVIAR are Linux ELF
# binaries, and finemapr calls them through system(). These .cmd shims let
# Windows R invoke them transparently.
#
# This works because finemapr::run_finemap()/run_caviar() setwd() into the run
# directory and then pass only *relative* filenames (region.z, region.ld,
# region.master). So the shim only has to translate the working directory - the
# arguments need no path rewriting at all.
#
# %~dp0 keeps the shims relocatable: they resolve the sibling Linux binary
# rather than hardcoding an absolute path.

echo
echo "--- [2b] Windows -> WSL shims ---"

write_shim() {
  local name="$1"
  # CRLF line endings: cmd.exe is happier with them.
  sed 's/$/\r/' > "${BIN_DIR}/${name}.cmd" <<EOF
@echo off
rem Auto-generated by src/finemapping/install_finemap_caviar.sh - do not edit.
rem Forwards to the Linux ${name} binary inside WSL, translating only the cwd.
setlocal enabledelayedexpansion
for /f "usebackq delims=" %%i in (\`wsl -d Ubuntu -e wslpath -a "%CD%"\`) do set "WSLCWD=%%i"
for /f "usebackq delims=" %%i in (\`wsl -d Ubuntu -e wslpath -a "%~dp0${name}"\`) do set "WSLBIN=%%i"
wsl -d Ubuntu -e bash -c "cd '!WSLCWD!' && '!WSLBIN!' %*"
exit /b %ERRORLEVEL%
EOF
  echo "wrote ${BIN_DIR}/${name}.cmd"
}

write_shim finemap
write_shim CAVIAR

# ---------------------------------------------------------------------------
# 3. Smoke tests
# ---------------------------------------------------------------------------
echo
echo "--- smoke tests ---"
TMP_TEST="$(mktemp -d)"
trap 'rm -rf "${TMP_TEST}"' EXIT

if [[ -d "${TOOLS_DIR}/finemap_example" ]]; then
  cp "${TOOLS_DIR}"/finemap_example/region1.z "${TMP_TEST}/" 2>/dev/null || true
  cp "${TOOLS_DIR}"/finemap_example/region1.ld "${TMP_TEST}/" 2>/dev/null || true
fi

cd "${TMP_TEST}"

# FINEMAP: build a v1.1 master file exactly like finemapr::run_finemap() does.
if [[ -f region1.z && -f region1.ld ]]; then
  printf 'z;ld;snp;config;log;n-ind\nregion1.z;region1.ld;region1.snp;region1.config;region1.log;5363\n' \
    > region1.master
  if "${BIN_DIR}/finemap" --sss --log --in-files region1.master > finemap_stdout.txt 2>&1; then
    echo "FINEMAP smoke test: PASS ($(wc -l < region1.snp) lines in region1.snp)"
  else
    echo "FINEMAP smoke test: FAIL - see output below"
    tail -20 finemap_stdout.txt
  fi

  # CAVIAR: -z takes the z-scores, -l the LD matrix.
  # NOTE: example 1 in finemapr's install doc has -z/-l swapped, which is why
  # it segfaults there. The correct order is used here.
  if "${BIN_DIR}/CAVIAR" -c 2 -z region1.z -l region1.ld -o caviar_log > caviar_stdout.txt 2>&1; then
    echo "CAVIAR smoke test:  PASS (causal set: $(tr '\n' ' ' < caviar_log_set))"
  else
    echo "CAVIAR smoke test:  FAIL - see output below"
    tail -20 caviar_stdout.txt
  fi
else
  echo "no example data found - skipping smoke tests"
fi

echo
echo "=============================================================="
echo "Done. Binaries:"
ls -la "${BIN_DIR}"
echo
echo "Point the finemapr R package at them with:"
echo "  options(finemapr_finemap = '<project>/tools/bin/finemap.cmd')"
echo "  options(finemapr_caviar  = '<project>/tools/bin/CAVIAR.cmd')"
echo "(.cmd shims are for Windows R; on Linux use the bare binaries.)"
echo "=============================================================="
