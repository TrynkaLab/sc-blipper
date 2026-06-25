#!/usr/bin/env python3
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def read_gsa(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(line.split())
    if not rows:
        raise ValueError(f"No data found in {path}")
    return pd.DataFrame(rows[1:], columns=rows[0])


def main():
    parser = argparse.ArgumentParser(
        description="Plot permutation null distributions for MAGMA gene-covar results. "
                    "Produces one page per condition in a multi-page PDF."
    )
    parser.add_argument("--stats", required=True, help="Permutation stats TSV (from calc_permutation_stats.py)")
    parser.add_argument("--permuted-gsa", required=True, help="Permuted MAGMA gsa.out file")
    parser.add_argument("--output", required=True, help="Output PDF")
    parser.add_argument("--trait", required=True, help="Trait name (for plot titles)")
    parser.add_argument("--database", required=True, help="Database/matrix name (for plot titles)")
    args = parser.parse_args()

    stats = pd.read_csv(args.stats, sep="\t")
    perm = read_gsa(args.permuted_gsa)
    perm["BETA"] = perm["BETA"].astype(float)

    with PdfPages(args.output) as pdf:
        for _, row in stats.iterrows():
            condition = str(row["condition"])
            real_beta = float(row["real_beta"])
            emp_p = float(row["empirical_pval"])
            norm_p = float(row["normal_pval"])
            n_perm = int(row["n_permutations"])

            pattern = re.compile(rf"^{re.escape(condition)}_perm_\d+$")
            perm_betas = perm.loc[perm["VARIABLE"].str.match(pattern), "BETA"].values

            if len(perm_betas) == 0:
                continue

            fig, ax = plt.subplots(figsize=(7, 4))
            ax.hist(perm_betas, bins=50, color="steelblue", alpha=0.75, label="Null distribution")
            ax.axvline(real_beta, color="crimson", linewidth=2,
                       label=f"Real effect (β = {real_beta:.4f})")
            ax.set_xlabel("BETA (effect size)")
            ax.set_ylabel("Count")
            ax.set_title(f"{args.trait}  |  {condition}")
            ax.legend(fontsize=8)

            ann = (
                f"Empirical p = {emp_p:.4g}\n"
                f"Normal p    = {norm_p:.4g}\n"
                f"n perm      = {n_perm}"
            )
            ax.text(0.98, 0.97, ann, transform=ax.transAxes, fontsize=8,
                    verticalalignment="top", horizontalalignment="right",
                    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Written {len(stats)} plots to {args.output}")


if __name__ == "__main__":
    main()
