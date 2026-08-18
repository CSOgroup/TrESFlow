#!/usr/bin/env bash

set -euo pipefail

repository="${1:?repository directory is required}"
manifest_fallback="${2:-}"

# A tagged checkout is a released build. Development checkouts retain the
# nearest release as their baseline and add an offline, reproducible commit
# distance/hash suffix. Deliberately ignore working-tree dirtiness: it is not a
# release-version source and would make otherwise identical runs report
# different pipeline versions.
repository="$(cd "${repository}" && pwd -P)"
git_root="$(git -C "${repository}" rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -n "${git_root}" ]] && [[ "$(cd "${git_root}" && pwd -P)" == "${repository}" ]]; then
    exact_tag="$(git -C "${repository}" describe --tags --exact-match --match 'v[0-9]*' HEAD 2>/dev/null || true)"
    if [[ "${exact_tag}" =~ ^v[0-9]+\.[0-9]+\.[0-9]+([.-][0-9A-Za-z]+)*$ ]]; then
        printf '%s\n' "${exact_tag}"
        exit 0
    fi

    described="$(git -C "${repository}" describe --tags --long --abbrev=7 --match 'v[0-9]*' HEAD 2>/dev/null || true)"
    if [[ "${described}" == *-g* ]]; then
        commit_hash="${described##*-g}"
        prefix="${described%-g*}"
        commit_count="${prefix##*-}"
        release_tag="${prefix%-*}"
        if [[ "${release_tag}" =~ ^v[0-9]+\.[0-9]+\.[0-9]+([.-][0-9A-Za-z]+)*$ ]] \
            && [[ "${commit_count}" =~ ^[0-9]+$ ]] \
            && [[ "${commit_hash}" =~ ^[0-9a-f]{7,}$ ]]; then
            if [[ "${commit_count}" == "0" ]]; then
                printf '%s\n' "${release_tag}"
            else
                printf '%s+%s.g%s\n' "${release_tag}" "${commit_count}" "${commit_hash}"
            fi
            exit 0
        fi
    fi
fi

# Source archives may not carry .git. In that case use the repository's
# manifest metadata without converting a development version into a release.
manifest_fallback="${manifest_fallback//[[:space:]]/}"
if [[ "${manifest_fallback}" =~ ^v[0-9]+\.[0-9]+\.[0-9]+([.-]?[0-9A-Za-z]+)*$ ]]; then
    printf '%s\n' "${manifest_fallback}"
elif [[ "${manifest_fallback}" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-]?[0-9A-Za-z]+)*$ ]]; then
    printf 'v%s\n' "${manifest_fallback}"
else
    printf '%s\n' 'unreleased'
fi
