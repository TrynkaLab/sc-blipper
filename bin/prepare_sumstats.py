#!/usr/bin/env python3
"""
Standardize GWAS summary statistics to a consistent column format.

Reads a (optionally gzipped) summary statistics file, maps user-specified
column names to standardised names, and computes Z = beta / SE when a direct
Z-score column is not available.

Odds ratios are converted to log(OR) to produce BETA before Z is computed.

Output column order (columns absent from the input are omitted):
    SNP  CHR  POS  A2  A1  BETA  SE  Z  P  N  N_CAS  N_CON
"""
import argparse
import numpy as np
import pandas as pd

OUTPUT_COLS = ["SNP", "CHR", "POS", "A2", "A1", "BETA", "SE", "Z", "P", "N", "N_CAS", "N_CON"]


def read_sumstats(path):
    """Read a sumstats file, inferring separator and compression."""
    compression = "infer"
    for sep in ("\t", " ", ","):
        try:
            probe = pd.read_csv(path, sep=sep, compression=compression, dtype=str, nrows=2)
            if probe.shape[1] > 1:
                return pd.read_csv(path, sep=sep, compression=compression, dtype=str, low_memory=False)
        except Exception:
            continue
    return pd.read_csv(path, sep=None, engine="python", compression=compression, dtype=str, low_memory=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",     required=True, help="Input sumstats file (TSV/CSV, optionally .gz)")
    parser.add_argument("--output",    required=True, help="Output TSV file")
    parser.add_argument("--snp",       default="NA",  help="Column for SNP/variant ID")
    parser.add_argument("--chr",       default="NA",  help="Column for chromosome")
    parser.add_argument("--pos",       default="NA",  help="Column for base-pair position")
    parser.add_argument("--a1",        default="NA",  help="Column for effect allele (A1)")
    parser.add_argument("--a2",        default="NA",  help="Column for other allele (A2)")
    parser.add_argument("--z",         default="NA",  help="Column for Z-score; computed from beta/SE if NA")
    parser.add_argument("--beta",      default="NA",  help="Column for beta; used to compute Z if --z is NA")
    parser.add_argument("--or",        default="NA",  help="Column for odds ratio; log(OR) is used as beta if --beta is NA")
    parser.add_argument("--se",        default="NA",  help="Column for SE; used to compute Z if --z is NA")
    parser.add_argument("--p",         default="NA",  help="Column for p-value")
    parser.add_argument("--n-col",     default="NA",  help="Column for sample size N")
    parser.add_argument("--n",         default="NA",  help="Fixed sample size N")
    parser.add_argument("--n-cas",     default="NA",  help="Fixed N cases")
    parser.add_argument("--n-con",     default="NA",  help="Fixed N controls")
    args = parser.parse_args()

    log_path = args.output.rsplit(".", 1)[0] + ".log"
    lf = open(log_path, "w")

    def log(msg):
        print(msg, file=lf, flush=True)
        print(msg, flush=True)

    log(f"Reading {args.input}")
    df = read_sumstats(args.input)
    log(f"Loaded {len(df):,} rows x {df.shape[1]} columns")
    log(f"Input columns: {', '.join(df.columns)}")

    out = pd.DataFrame(index=df.index)

    # ------------------------------------------------------------------
    # Simple column renames
    # ------------------------------------------------------------------
    for std, src in [("SNP",  args.snp),
                     ("CHR",  args.chr),
                     ("POS",  args.pos),
                     ("A2",   args.a2),
                     ("A1",   args.a1),
                     ("SE",   args.se),
                     ("P",    args.p)]:
        if src == "NA":
            continue
        if src not in df.columns:
            log(f"WARNING: column '{src}' not found — {std} will be absent from output")
            continue
        out[std] = df[src]

    # ------------------------------------------------------------------
    # BETA — prefer explicit beta column, fall back to log(OR)
    # ------------------------------------------------------------------
    or_col = getattr(args, "or")  # 'or' is a reserved word, access via getattr
    if args.beta != "NA":
        if args.beta not in df.columns:
            log(f"WARNING: beta column '{args.beta}' not found")
        else:
            out["BETA"] = pd.to_numeric(df[args.beta], errors="coerce")
            log(f"Using beta column '{args.beta}'")

    if "BETA" not in out.columns and or_col != "NA":
        if or_col not in df.columns:
            log(f"WARNING: OR column '{or_col}' not found")
        else:
            or_vals = pd.to_numeric(df[or_col], errors="coerce")
            out["BETA"] = np.log(or_vals)
            log(f"Computed BETA = log({or_col})")

    # ------------------------------------------------------------------
    # Z-score
    # ------------------------------------------------------------------
    if args.z != "NA":
        if args.z not in df.columns:
            log(f"WARNING: Z column '{args.z}' not found — will attempt to compute from beta/SE")
        else:
            out["Z"] = pd.to_numeric(df[args.z], errors="coerce")
            log(f"Using Z-score column '{args.z}'")

    if "Z" not in out.columns:
        if "BETA" in out.columns and "SE" in out.columns:
            beta = pd.to_numeric(out["BETA"], errors="coerce")
            se   = pd.to_numeric(out["SE"],   errors="coerce")
            out["Z"] = beta / se
            log(f"Computed Z = BETA / SE")
        else:
            log("WARNING: no Z column and beta/SE unavailable — Z will be absent from output")

    # ------------------------------------------------------------------
    # Sample sizes
    # ------------------------------------------------------------------
    if args.n_col != "NA":
        if args.n_col not in df.columns:
            log(f"WARNING: N column '{args.n_col}' not found")
        else:
            out["N"] = pd.to_numeric(df[args.n_col], errors="coerce")
            log(f"Using N column '{args.n_col}'")
    elif args.n != "NA":
        out["N"] = args.n
        log(f"Using fixed N = {args.n}")

    if args.n_cas != "NA":
        out["N_CAS"] = args.n_cas
        log(f"Using fixed N_CAS = {args.n_cas}")

    if args.n_con != "NA":
        out["N_CON"] = args.n_con
        log(f"Using fixed N_CON = {args.n_con}")

    # ------------------------------------------------------------------
    # Write in standard column order, dropping absent columns
    # ------------------------------------------------------------------
    keep = [c for c in OUTPUT_COLS if c in out.columns]
    out[keep].to_csv(args.output, sep="\t", index=False)
    log(f"Written {len(out):,} rows, columns: {', '.join(keep)} → {args.output}")
    lf.close()


if __name__ == "__main__":
    main()
