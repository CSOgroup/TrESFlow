#!/usr/bin/env bash
set -euo pipefail

CODON_VERSION="0.16.3"
SEQ_VERSION="0.11.3"
CODON_REPO="https://github.com/exaloop/codon/releases/download/v${CODON_VERSION}"
SEQ_REPO="https://github.com/exaloop/seq/releases/download/v${SEQ_VERSION}"

usage() {
  cat <<'EOF'
Install Codon 0.16.3 and Seq 0.11.3 into an environment prefix.

Usage:
  scripts/install_codon_0.16.3.sh [--prefix <path>]
  scripts/install_codon_0.16.3.sh <path>

Resolution order:
  1. Explicit --prefix / positional argument
  2. $CONDA_PREFIX

The installer writes:
  <prefix>/bin/codon
  <prefix>/lib/codon/plugins/seq/
EOF
}

resolve_prefix() {
  local explicit_prefix=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --prefix)
        [[ $# -ge 2 ]] || { echo "error: --prefix requires a value" >&2; exit 1; }
        explicit_prefix="$2"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        if [[ -n "${explicit_prefix}" ]]; then
          echo "error: only one prefix may be specified" >&2
          exit 1
        fi
        explicit_prefix="$1"
        shift
        ;;
    esac
  done

  if [[ -n "${explicit_prefix}" ]]; then
    printf '%s\n' "${explicit_prefix}"
    return 0
  fi

  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    printf '%s\n' "${CONDA_PREFIX}"
    return 0
  fi

  echo "error: no install prefix provided and CONDA_PREFIX is unset" >&2
  usage >&2
  exit 1
}

resolve_platform() {
  local os_raw arch_raw os arch
  os_raw="$(uname -s)"
  arch_raw="$(uname -m)"

  case "${os_raw}" in
    Linux)  os="linux" ;;
    Darwin) os="darwin" ;;
    *)
      echo "error: unsupported operating system '${os_raw}'. Supported: Linux, Darwin." >&2
      exit 1
      ;;
  esac

  case "${arch_raw}" in
    x86_64|amd64) arch="x86_64" ;;
    arm64|aarch64) arch="arm64" ;;
    *)
      echo "error: unsupported architecture '${arch_raw}'. Supported: x86_64, arm64." >&2
      exit 1
      ;;
  esac

  case "${os}-${arch}" in
    linux-x86_64|darwin-x86_64|darwin-arm64)
      printf '%s %s\n' "${os}" "${arch}"
      ;;
    *)
      echo "error: no prebuilt Codon/Seq release archives are available for ${os}-${arch}." >&2
      exit 1
      ;;
  esac
}

extract_semver() {
  local value="${1:-}"
  printf '%s\n' "${value}" | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -n 1 || true
}

validate_install() {
  local prefix="$1"
  local codon_bin="${prefix}/bin/codon"
  local seq_plugin_toml="${prefix}/lib/codon/plugins/seq/plugin.toml"
  local smoke_file="$2"
  local codon_version_raw codon_version seq_version_raw seq_version

  [[ -x "${codon_bin}" ]] || {
    echo "error: installed Codon binary is missing or not executable: ${codon_bin}" >&2
    exit 1
  }

  [[ -f "${seq_plugin_toml}" ]] || {
    echo "error: installed Seq plugin metadata is missing: ${seq_plugin_toml}" >&2
    exit 1
  }

  codon_version_raw="$("${codon_bin}" --version 2>/dev/null | head -n 1 || true)"
  codon_version="$(extract_semver "${codon_version_raw}")"
  if [[ "${codon_version}" != "${CODON_VERSION}" ]]; then
    echo "error: expected Codon ${CODON_VERSION}, found '${codon_version_raw:-unknown}'" >&2
    exit 1
  fi

  seq_version_raw="$(awk -F'"' '/^version = /{print $2; exit}' "${seq_plugin_toml}")"
  seq_version="$(extract_semver "${seq_version_raw}")"
  if [[ "${seq_version}" != "${SEQ_VERSION}" ]]; then
    echo "error: expected Seq ${SEQ_VERSION}, found '${seq_version_raw:-unknown}'" >&2
    exit 1
  fi

  "${codon_bin}" run -plugin seq -release "${smoke_file}" >/dev/null
}

PREFIX="$(resolve_prefix "$@")"
read -r OS ARCH < <(resolve_platform)

CODON_ARCHIVE="codon-${OS}-${ARCH}.tar.gz"
SEQ_ARCHIVE="seq-${OS}-${ARCH}.tar.gz"
CODON_URL="${CODON_REPO}/${CODON_ARCHIVE}"
SEQ_URL="${SEQ_REPO}/${SEQ_ARCHIVE}"

mkdir -p "${PREFIX}"
PREFIX="$(cd "${PREFIX}" && pwd -P)"
SEQ_PLUGIN_DIR="${PREFIX}/lib/codon/plugins/seq"

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/codon-install.XXXXXX")"
trap 'rm -rf "${TMPDIR}"' EXIT

echo "Installing Codon ${CODON_VERSION} into ${PREFIX}"
echo "Installing Seq ${SEQ_VERSION} into ${SEQ_PLUGIN_DIR}"
echo "Using release archives:"
echo "  ${CODON_URL}"
echo "  ${SEQ_URL}"

curl -fsSL "${CODON_URL}" -o "${TMPDIR}/${CODON_ARCHIVE}"
curl -fsSL "${SEQ_URL}" -o "${TMPDIR}/${SEQ_ARCHIVE}"

tar -xzf "${TMPDIR}/${CODON_ARCHIVE}" -C "${PREFIX}" --strip-components=1
mkdir -p "${SEQ_PLUGIN_DIR}"
tar -xzf "${TMPDIR}/${SEQ_ARCHIVE}" -C "${SEQ_PLUGIN_DIR}" --strip-components=1

cat > "${TMPDIR}/seq_smoke.codon" <<'EOF'
import bio
print("seq-ok")
EOF

validate_install "${PREFIX}" "${TMPDIR}/seq_smoke.codon"

echo "Install validation passed:"
echo "  codon_path=${PREFIX}/bin/codon"
echo "  codon_version=${CODON_VERSION}"
echo "  seq_plugin_dir=${SEQ_PLUGIN_DIR}"
echo "  seq_version=${SEQ_VERSION}"
