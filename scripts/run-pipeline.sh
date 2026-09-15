#!/usr/bin/env bash

# Use strict error handling: stop on failures, unset variables, and failed
set -Eeuo pipefail

# Find the repository root from this script's location so the pipeline works
# even when it is started from a different directory.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT_DIR/scripts"
OUTPUT_DIR="$ROOT_DIR/outputs"
LOG_DIR="$ROOT_DIR/logs"
RUNTIME_LOG="$OUTPUT_DIR/pipeline-runtime.csv"

# Create the output directories if they do not already exist.
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Keep a history of pipeline runs. The header is written only for a new log.
if [[ ! -f "$RUNTIME_LOG" ]]; then
    printf 'timestamp,script,runtime,status\n' > "$RUNTIME_LOG"
fi

# Exit with a helpful message when a required program is unavailable.
require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        printf 'Error: required command not found: %s\n' "$1" >&2
        exit 1
    fi
}

# Print a label, run the command, and record its elapsed time in the CSV log.
run_step() {
    local script="$1"
    shift
    local start_time
    local end_time
    local start_seconds
    local end_seconds
    local runtime_seconds
    local runtime_minutes
    local runtime_remainder
    local runtime
    local status

    start_time="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    start_seconds="$(date +%s)"

    printf '\n==> %s\n' "$script"
    if "$@"; then
        status="success"
    else
        status="failed"
    fi

    end_time="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    end_seconds="$(date +%s)"
    runtime_seconds=$((end_seconds - start_seconds))
    runtime_minutes=$((runtime_seconds / 60))
    runtime_remainder=$((runtime_seconds % 60))
    runtime="$(printf '%02d:%02d' "$runtime_minutes" "$runtime_remainder")"

    # Quote CSV fields so paths remain valid if they contain spaces or commas.
    printf '"%s","%s","%s","%s"\n' \
        "$start_time" "$script" "$runtime" "$status" >> "$RUNTIME_LOG"

    if [[ "$status" == "failed" ]]; then
        return 1
    fi
}

# Quarto is currently the only external program used by this pipeline.
require_command quarto

# Add data preparation, analysis, and reporting scripts here in dependency order.
# A trailing backslash continues one command across multiple lines for readability.

# TODO: Replace the log files with the task scripts once they are created. log-questions.qmd and log-issues.qmd are placeholders.
run_step "$LOG_DIR/log-questions.qmd" \
    quarto render "$LOG_DIR/log-questions.qmd" \
    --to html \
    --output-dir "$LOG_DIR"

run_step "$LOG_DIR/log-issues.qmd" \
    quarto render "$LOG_DIR/log-issues.qmd" \
    --to html \
    --output-dir "$LOG_DIR"

printf '\n====================\n'
printf 'Pipeline finished at: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
printf '\nPipeline completed successfully.\n'
printf '====================\n'