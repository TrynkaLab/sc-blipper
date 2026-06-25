#!/usr/bin/env python3
import argparse
import math
import re
import numpy as np
import pandas as pd


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


def normal_pval(z):
    """Two-tailed p-value from z-score using the complementary error function (no scipy needed)."""
    if math.isnan(z):
        return float("nan")
    return math.erfc(abs(z) / math.sqrt(2))


def main():
    parser = argparse.ArgumentParser(
        description="Compute permutation statistics for MAGMA gene-covar results. "
                    "Expects a real gsa.out and a permuted gsa.out where permuted VARIABLE names "
                    "follow the pattern {condition}_perm_{i}."
    )
    parser.add_argument("--real-gsa", required=True, help="Real MAGMA gsa.out file")
    parser.add_argument("--permuted-gsa", required=True, help="Permuted MAGMA gsa.out file")
    parser.add_argument("--output", required=True, help="Output TSV with permutation statistics")
    parser.add_argument("--trait", required=True, help="Trait name (for output labelling)")
    parser.add_argument("--database", required=True, help="Database/matrix name (for output labelling)")
    args = parser.parse_args()

    real = read_gsa(args.real_gsa)
    perm = read_gsa(args.permuted_gsa)

    real["BETA"] = real["BETA"].astype(float)
    perm["BETA"] = perm["BETA"].astype(float)

    records = []
    for _, row in real.iterrows():
        condition = str(row["VARIABLE"])
        real_beta = float(row["BETA"])

        pattern = re.compile(rf"^{re.escape(condition)}_perm_\d+$")
        perm_betas = perm.loc[perm["VARIABLE"].str.match(pattern), "BETA"].values

        if len(perm_betas) == 0:
            print(f"WARNING: no permuted results found for condition '{condition}', skipping")
            continue

        perm_mean = float(np.mean(perm_betas))
        perm_sd = float(np.std(perm_betas, ddof=1))
        z_score = (real_beta - perm_mean) / perm_sd if perm_sd > 0 else float("nan")
        n = len(perm_betas)
        r = int(np.sum(np.abs(perm_betas) >= np.abs(real_beta)))
        empirical_pval = (r + 1) / (n + 1)
        norm_pval = normal_pval(z_score)

        records.append({
            "trait": args.trait,
            "database": args.database,
            "condition": condition,
            "real_beta": real_beta,
            "perm_mean": perm_mean,
            "perm_sd": perm_sd,
            "z_score": z_score,
            "empirical_pval": empirical_pval,
            "normal_pval": norm_pval,
            "n_permutations": len(perm_betas),
        })

    out = pd.DataFrame(records)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Written permutation stats for {len(out)} conditions to {args.output}")


if __name__ == "__main__":
    main()
