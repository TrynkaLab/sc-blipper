#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_CACHE_DIR="${TEST_CACHE_DIR:-${PROJECT_ROOT}/tests/.cache}"
LOG_DIR="${TEST_LOG_DIR:-${TEST_CACHE_DIR}/logs}"
RUN_ID="$(date '+%Y%m%d-%H%M%S')"
LOG_FILE="${TEST_LOG_FILE:-${LOG_DIR}/test-all-${RUN_ID}.log}"
NF_TEST_LOG_FILE="${NF_TEST_LOG_FILE:-${LOG_DIR}/nf-test-${RUN_ID}.log}"

PYTHON_BIN="${PYTHON:-python}"
NF_TEST_BIN="${NF_TEST:-nf-test}"

mkdir -p "${LOG_DIR}" "${TEST_CACHE_DIR}/pytest"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

require_command() {
    local label="$1"
    local cmd="$2"

    if ! command -v "${cmd}" >/dev/null 2>&1; then
        log "ERROR: ${label} command not found: ${cmd}"
        log "Set ${label}=<path> before running this script if it is not on PATH."
        exit 127
    fi
}

run_stage() {
    local name="$1"
    shift

    log "START ${name}"
    log "COMMAND $*"
    local start
    start="$(date '+%s')"

    set +e
    "$@"
    local status=$?
    set -e

    local end
    end="$(date '+%s')"
    log "END ${name} status=${status} elapsed=$((end - start))s"
    return "${status}"
}

main() {
    cd "${PROJECT_ROOT}"

    log "Test run started"
    log "Project root: ${PROJECT_ROOT}"
    log "Log file: ${LOG_FILE}"
    log "nf-test log file: ${NF_TEST_LOG_FILE}"

    require_command PYTHON "${PYTHON_BIN}"

    # Specific/unit-level tests first: fast feedback on standalone Python logic.
    run_stage \
        "python specific tests" \
        "${PYTHON_BIN}" -B -m pytest -vv -ra tests/python

    # Process-level tests second: slower Nextflow/nf-test coverage after unit checks pass.
    require_command NF_TEST "${NF_TEST_BIN}"
    run_stage \
        "nextflow process tests" \
        "${NF_TEST_BIN}" test tests/nextflow/cnmf_smoke.nf.test --verbose --log "${NF_TEST_LOG_FILE}"

    log "All test stages passed"
}

set +e
main "$@" 2>&1 | tee -a "${LOG_FILE}"
status="${PIPESTATUS[0]}"
set -e
exit "${status}"
