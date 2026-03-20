#!/usr/bin/env bash
set -euo pipefail

REQUIRED_CODON_VERSION="0.16.3"
REQUIRED_SEQ_VERSION="0.11.3"
CODON_HOME_DIR="${CODON_HOME:-${HOME}/.codon}"
CODON_BIN_CONFIGURED="${CODON_BIN:-}"
SEQ_PLUGIN_DIR="${CODON_HOME_DIR}/lib/codon/plugins/seq"
SEQ_PLUGIN_TOML="${SEQ_PLUGIN_DIR}/plugin.toml"

extract_semver() {
  local value="${1:-}"
  printf '%s\n' "${value}" | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true
}

echo "Checking pinned Codon/Seq host prerequisites for every pipeline run"
echo "CODON_HOME=${CODON_HOME_DIR}"
echo "required_codon_version=${REQUIRED_CODON_VERSION}"
echo "required_seq_version=${REQUIRED_SEQ_VERSION}"

if [[ -n "${CODON_BIN_CONFIGURED}" ]]; then
  if [[ ! -x "${CODON_BIN_CONFIGURED}" ]]; then
    echo "ERROR: configured CODON_BIN is missing or not executable: ${CODON_BIN_CONFIGURED}" >&2
    exit 1
  fi
  CODON_BIN="${CODON_BIN_CONFIGURED}"
  CODON_PATH_SOURCE="configured"
else
  if ! command -v codon >/dev/null 2>&1; then
    echo "ERROR: 'codon' is not on PATH and CODON_BIN is not configured." >&2
    echo "Install Codon ${REQUIRED_CODON_VERSION} with scripts/install_codon_0.16.3.sh and ensure the binary is on PATH, or configure CODON_BIN explicitly." >&2
    exit 1
  fi
  CODON_BIN="$(command -v codon)"
  CODON_PATH_SOURCE="PATH"
fi

CODON_VERSION_RAW="$("${CODON_BIN}" --version 2>/dev/null | head -n 1 || true)"
CODON_VERSION="$(extract_semver "${CODON_VERSION_RAW}")"

echo "codon_path=${CODON_BIN}"
echo "codon_path_source=${CODON_PATH_SOURCE}"
echo "codon_version_raw=${CODON_VERSION_RAW:-unknown}"
echo "codon_version=${CODON_VERSION:-unknown}"

if [[ "${CODON_VERSION}" != "${REQUIRED_CODON_VERSION}" ]]; then
  echo "ERROR: Codon ${REQUIRED_CODON_VERSION} is required for every pipeline run, found '${CODON_VERSION_RAW:-unknown}'." >&2
  exit 1
fi

if [[ ! -f "${SEQ_PLUGIN_TOML}" ]]; then
  echo "ERROR: Seq plugin metadata not found at ${SEQ_PLUGIN_TOML}" >&2
  echo "Install Seq ${REQUIRED_SEQ_VERSION} into ${CODON_HOME_DIR}/lib/codon/plugins/seq" >&2
  exit 1
fi

SEQ_VERSION_RAW="$(awk -F'"' '/^version = /{print $2; exit}' "${SEQ_PLUGIN_TOML}")"
SEQ_VERSION="$(extract_semver "${SEQ_VERSION_RAW}")"
SEQ_SUPPORTED_CODON="$(awk -F'"' '/^supported = /{print $2; exit}' "${SEQ_PLUGIN_TOML}")"

echo "seq_plugin_dir=${SEQ_PLUGIN_DIR}"
echo "seq_version_raw=${SEQ_VERSION_RAW:-unknown}"
echo "seq_version=${SEQ_VERSION:-unknown}"
echo "seq_supported_codon=${SEQ_SUPPORTED_CODON:-unknown}"

if [[ "${SEQ_VERSION}" != "${REQUIRED_SEQ_VERSION}" ]]; then
  echo "ERROR: Seq ${REQUIRED_SEQ_VERSION} is required for every pipeline run, found '${SEQ_VERSION_RAW:-unknown}'." >&2
  exit 1
fi

if [[ ! -f "${SEQ_PLUGIN_DIR}/build/libseq.dylib" && ! -f "${SEQ_PLUGIN_DIR}/build/libseq.so" ]]; then
  echo "WARNING: Seq plugin shared library not found under ${SEQ_PLUGIN_DIR}/build" >&2
fi

echo "Host preflight passed for globally pinned pipeline toolchain."
