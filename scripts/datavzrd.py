"""Render a datavzrd report from a yte template.

Port of the snakemake-wrapper-utils v9.14.0 utils/datavzrd wrapper (MIT,
Copyright 2017, Johannes Koester), adapted to a plain argv interface:

  datavzrd.py TEMPLATE.yaml VARIABLES.json --output OUTDIR --log LOG

VARIABLES.json carries the wrapper's `variables` dict (params, wildcards,
input). Pandas DataFrames are encoded as {"columns": [...], "data": [...]}
(plus optional "index") and rebuilt before yte processing, matching the
upstream wrapper semantics (process_yaml with require_use_yte=True, then
`datavzrd {processed} --output {output}`).
"""

import argparse
import json
import subprocess
from types import SimpleNamespace
import sys
import tempfile

import pandas as pd
from yte import process_yaml

DATA_FRAME_KEYS = {"samples", "group_annotations", "labels"}


def rebuild_frame(encoded):
    frame = pd.DataFrame(encoded["data"], columns=encoded.get("columns", []))
    if encoded.get("index"):
        frame.set_index(encoded["index"], inplace=True)
    return frame


def to_namespace(value):
    # The yte templates use attribute access (?input.csv, ?wildcards.group) —
    # plain JSON dicts do not support it (live: yte YteError "dict object
    # has no attribute csv"). Recursively wrap dicts in SimpleNamespace.
    # Frame-encoded dicts ({"columns", "data"}) become DataFrames at any
    # nesting level — the templates use pandas .loc on them (live:
    # SimpleNamespace has no attribute loc on params.samples).
    if isinstance(value, dict):
        if set(value.keys()) >= {"columns", "data"}:
            return rebuild_frame(value)
        return SimpleNamespace(**{k: to_namespace(v) for k, v in value.items()})
    if isinstance(value, list):
        return [to_namespace(v) for v in value]
    return value


def rebuild_variables(variables):
    for key, value in variables.items():
        if isinstance(value, dict) and set(value.keys()) >= {"columns", "data"}:
            variables[key] = rebuild_frame(value)
        else:
            variables[key] = to_namespace(value)
    return variables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", help="yte template (datavzrd config)")
    parser.add_argument("variables", help="variables JSON (params/wildcards/input)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--log", required=True, help="Log file (stdout/stderr redirect)")
    args = parser.parse_args()

    with open(args.variables) as f:
        variables = rebuild_variables(json.load(f))

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml") as processed, open(
        args.template
    ) as f:
        with open(args.log, "w") as log_file:
            old_stdout, old_stderr = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = log_file
            try:
                # support templating in the config file
                process_yaml(
                    f,
                    outfile=processed,
                    variables=variables,
                    require_use_yte=True,
                )
            finally:
                sys.stdout, sys.stderr = old_stdout, old_stderr
        processed.flush()

        subprocess.run(
            ["datavzrd", processed.name, "--output", args.output],
            check=True,
            stdout=open(args.log, "a"),
            stderr=subprocess.STDOUT,
        )


if __name__ == "__main__":
    main()
