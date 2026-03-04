#!/usr/bin/env python3
"""Fetch OmniPath-backed decoupler networks and save them as TSV files."""

import argparse
import decoupler as dc


def fetch_progeny(organism: str, top: int):
    try:
        return dc.op.progeny(organism=organism, top=top)
    except TypeError:
        # Compatibility with decoupler versions that do not expose `top`.
        return dc.op.progeny(organism=organism)


def fetch_collectri(organism: str):
    return dc.op.collectri(organism=organism)


def main():
    parser = argparse.ArgumentParser(
        description="Fetch PROGENy and CollectTRI networks from decoupler/OmniPath and save as TSV."
    )
    parser.add_argument(
        "--organism",
        default="human",
        help="Organism passed to decoupler (default: human).",
    )
    parser.add_argument(
        "--progeny-top",
        type=int,
        default=500,
        help="Top genes per pathway for PROGENy when supported (default: 500).",
    )
    parser.add_argument(
        "--progeny-out",
        default="progeny_network.tsv",
        help="Output TSV path for PROGENy network.",
    )
    parser.add_argument(
        "--collectri-out",
        default="collectri_network.tsv",
        help="Output TSV path for CollectTRI network.",
    )

    args = parser.parse_args()

    print(f"Fetching PROGENy for organism={args.organism}...")
    progeny = fetch_progeny(args.organism, args.progeny_top)
    progeny.to_csv(args.progeny_out, sep="\t", index=False)
    print(f"Wrote PROGENy TSV: {args.progeny_out} ({len(progeny)} rows)")

    print(f"Fetching CollectTRI for organism={args.organism}...")
    collectri = fetch_collectri(args.organism)
    collectri.to_csv(args.collectri_out, sep="\t", index=False)
    print(f"Wrote CollectTRI TSV: {args.collectri_out} ({len(collectri)} rows)")


if __name__ == "__main__":
    main()
