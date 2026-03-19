#!/usr/bin/env bash
set -euo pipefail

CODON_HOME_DIR="${CODON_HOME:-${HOME}/.codon}"
SEQ_PLUGIN_DIR="${CODON_HOME_DIR}/lib/codon/plugins/seq"
SEQ_PLUGIN_TOML="${SEQ_PLUGIN_DIR}/plugin.toml"

echo "Checking Codon/Seq host prerequisites"
echo "CODON_HOME=${CODON_HOME_DIR}"

if ! command -v codon >/dev/null 2>&1; then
  echo "ERROR: 'codon' is not on PATH." >&2
  echo "Install Codon with Exaloop's installer script and ensure the binary is on PATH." >&2
  exit 1
fi

CODON_BIN="$(command -v codon)"
CODON_VERSION="$(codon --version 2>/dev/null | head -n 1 || true)"

echo "codon_path=${CODON_BIN}"
echo "codon_version=${CODON_VERSION:-unknown}"

if [[ ! -f "${SEQ_PLUGIN_TOML}" ]]; then
  echo "ERROR: Seq plugin metadata not found at ${SEQ_PLUGIN_TOML}" >&2
  echo "Install the Seq plugin tarball into ${CODON_HOME_DIR}/lib/codon/plugins/seq" >&2
  exit 1
fi

SEQ_VERSION="$(awk -F'"' '/^version = /{print $2; exit}' "${SEQ_PLUGIN_TOML}")"
SEQ_SUPPORTED_CODON="$(awk -F'"' '/^supported = /{print $2; exit}' "${SEQ_PLUGIN_TOML}")"

echo "seq_plugin_dir=${SEQ_PLUGIN_DIR}"
echo "seq_version=${SEQ_VERSION:-unknown}"
echo "seq_supported_codon=${SEQ_SUPPORTED_CODON:-unknown}"

if [[ ! -f "${SEQ_PLUGIN_DIR}/build/libseq.dylib" && ! -f "${SEQ_PLUGIN_DIR}/build/libseq.so" ]]; then
  echo "WARNING: Seq plugin shared library not found under ${SEQ_PLUGIN_DIR}/build" >&2
fi

echo "Host preflight passed."
