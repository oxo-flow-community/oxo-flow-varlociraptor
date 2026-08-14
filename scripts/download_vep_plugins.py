"""Download the VEP plugins zip for an Ensembl release and unpack it.

Port of the snakemake-wrapper-utils v8.0.0 bio/vep/plugins wrapper (MIT,
Copyright 2020, Johannes Koester), adapted to a plain argv interface:
  download_vep_plugins.py --release 111 --output resources/vep/plugins

The archive root directory is skipped; all members are written directly into
the output directory.
"""

import argparse
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlretrieve
from zipfile import ZipFile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="Ensembl release, e.g. 111")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--log", required=True, help="Log file (stderr redirect)")
    args = parser.parse_args()

    if args.log:
        sys.stderr = open(args.log, "w")

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    with NamedTemporaryFile() as tmp:
        urlretrieve(
            "https://github.com/Ensembl/VEP_plugins/archive/release/{release}.zip".format(
                release=args.release
            ),
            tmp.name,
        )

        with ZipFile(tmp.name) as f:
            for member in f.infolist():
                memberpath = Path(member.filename)
                if len(memberpath.parts) == 1:
                    # skip root dir
                    continue
                targetpath = outdir / memberpath.relative_to(memberpath.parts[0])
                if member.is_dir():
                    targetpath.mkdir()
                else:
                    with open(targetpath, "wb") as out:
                        out.write(f.read(member.filename))


if __name__ == "__main__":
    main()
