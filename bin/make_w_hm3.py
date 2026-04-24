import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rsids", required=True, help="File with HapMap3 RSIDs, one per line (or TSV with a 'SNP' column)")
parser.add_argument("--bim", required=True, nargs="+", help="PLINK .bim file(s), e.g. chr{1..22}.bim")
parser.add_argument("--out", default="w_hm3.snplist")
args = parser.parse_args()

# Load RSIDs
rsid_df = pd.read_csv(args.rsids, sep="\t")
if "SNP" in rsid_df.columns:
    hapmap_rsids = set(rsid_df["SNP"])
else:
    hapmap_rsids = set(rsid_df.iloc[:, 0])

# Load bim file(s)
bim_cols = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
bim = pd.concat(
    [pd.read_csv(f, sep="\t", header=None, names=bim_cols) for f in args.bim],
    ignore_index=True
)

# Filter to HapMap3 RSIDs and drop duplicates
out = (
    bim[bim["SNP"].isin(hapmap_rsids)][["SNP", "A1", "A2"]]
    .drop_duplicates(subset="SNP")
)

out.to_csv(args.out, sep="\t", index=False)
print(f"Written {len(out):,} SNPs to {args.out}")
