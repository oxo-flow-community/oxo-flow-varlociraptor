#!/usr/bin/env bash
# Acceptance test for the oxo-flow-varlociraptor port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
# oxo-flow v0.11.0 prints the plan to stderr; capture both streams
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards} ({log} stays literal)"
"$OXO" debug main.oxoflow 2>&1 | grep -E '\{(config\.|sample\}|group\}|input\[|output\[|threads\}|memory\})' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "PASS"
