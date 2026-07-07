#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
import yaml
from sklearn.decomposition import non_negative_factorization
from sklearn.metrics import r2_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def load_hvgs(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def load_spectra(path):
    return pd.read_csv(path, sep="\t", compression="infer", index_col=0)


def reconcile_genes(adata, hvgs, spectra, output_prefix):
    available = set(adata.var_names) & set(spectra.columns)
    missing = [g for g in hvgs if g not in available]

    if missing:
        print(f"WARNING: {len(missing)}/{len(hvgs)} HVGs missing from the h5ad/spectra overlap")
        missing_file = f"{output_prefix}.missing_genes.tsv"
        with open(missing_file, "w") as fh:
            fh.write("gene\n")
            fh.write("\n".join(missing) + "\n")
        print(f"Missing genes written to {missing_file}")

    return [g for g in hvgs if g in available]


def compute_gene_std(adata, genes):
    sub = adata[:, genes]
    if sp.issparse(sub.X):
        std = np.asarray(sub.X.todense()).std(axis=0, ddof=1)
    else:
        std = np.asarray(sub.X).std(axis=0, ddof=1)
    return pd.Series(std, index=genes)


def compute_baseline_usage(adata, hvgs, spectra_tpm, gene_std, nmf_params_path):
    norm_tpm = adata[:, hvgs].copy()
    if sp.issparse(norm_tpm.X):
        sc.pp.scale(norm_tpm, zero_center=False)
        norm_tpm_hvg = np.asarray(norm_tpm.X.todense())
    else:
        norm_tpm.X = norm_tpm.X / norm_tpm.X.std(axis=0, ddof=1)
        norm_tpm_hvg = np.asarray(norm_tpm.X)

    spectra_tpm_rf = spectra_tpm.loc[:, hvgs].div(gene_std.loc[hvgs], axis=1)

    refit_nmf_kwargs = yaml.load(open(nmf_params_path), Loader=yaml.FullLoader)
    refit_nmf_kwargs.update(dict(
        n_components=spectra_tpm_rf.shape[0],
        H=spectra_tpm_rf.values.astype(norm_tpm_hvg.dtype),
        update_H=False,
    ))
    usages, _, _ = non_negative_factorization(norm_tpm_hvg, **refit_nmf_kwargs)
    baseline_usages = pd.DataFrame(usages, index=adata.obs_names, columns=spectra_tpm_rf.index)

    return baseline_usages, norm_tpm_hvg, spectra_tpm_rf


def compute_coexpression(adata, genes):
    sub = adata[:, genes].copy()
    sc.pp.log1p(sub)
    X = np.asarray(sub.X.todense()) if sp.issparse(sub.X) else np.asarray(sub.X)
    corr = np.corrcoef(X, rowvar=False)
    return pd.DataFrame(corr, index=genes, columns=genes)


def cache_coexpression(coexp_df, output_prefix):
    cache_path = f"{output_prefix}.coexpression.h5ad"
    ad.AnnData(
        X=coexp_df.values,
        obs=pd.DataFrame(index=coexp_df.index),
        var=pd.DataFrame(index=coexp_df.columns),
    ).write_h5ad(cache_path)
    print(f"Coexpression matrix cached to {cache_path}")


def load_or_compute_coexpression(coexp_path, adata, genes, output_prefix):
    if coexp_path is not None:
        print(f"Loading precomputed coexpression matrix from {coexp_path}")
        coexp = sc.read_h5ad(coexp_path)
        return pd.DataFrame(
            np.asarray(coexp.X.todense()) if sp.issparse(coexp.X) else np.asarray(coexp.X),
            index=coexp.obs_names,
            columns=coexp.var_names,
        )

    print("Computing gene-gene coexpression matrix")
    coexp_df = compute_coexpression(adata, genes)
    cache_coexpression(coexp_df, output_prefix)
    return coexp_df


def rank_genes(score_row, ntop):
    ranked = score_row.sort_values(ascending=False)
    top = ranked.iloc[:ntop]
    return [(gene, rank + 1, score) for rank, (gene, score) in enumerate(top.items())]


def ablate_gene(spectra_row, coexp_row, hvgs, gene_std, norm_tpm_hvg, absolute):
    
    # Determine the ablation weights using the co-expression
    values = coexp_row.copy()
    if absolute:
        values = values.abs()

    sum_abs_coexp = coexp_row.abs().sum()
    n_high_coexp = int((values > 0.9).sum())

    clipped = values.clip(lower=0)
    weights = 1 - clipped

    stats = dict(
        sum_abs_coexp=sum_abs_coexp,
        sum_clipped_weights=weights.sum(),
        n_high_coexp=n_high_coexp,
    )

    # Ablate the spectra
    ablated_spectra = spectra_row * weights
    ablated_spectra_hvg = ablated_spectra.loc[hvgs]
    ablated_spectra_rf = ablated_spectra_hvg / gene_std.loc[hvgs]

    ablated_usage = norm_tpm_hvg.dot(ablated_spectra_rf.values)
    return ablated_usage, stats


def ols(x, y):
    """
    Fit y = intercept + slope * x by ordinary least squares and return R^2.
    Inputs:
      x, y : array-like (same length)
    Returns:
      r2 : float
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape")
    n = x.size
    if n < 2:
        raise ValueError("need at least two points")
    mx = x.mean()
    my = y.mean()
    SS_xy = np.sum((x - mx) * (y - my))
    SS_xx = np.sum((x - mx) ** 2)
    if SS_xx == 0:
        raise ValueError("x has zero variance; slope undefined")
    slope = SS_xy / SS_xx
    intercept = my - slope * mx
    y_pred = intercept + slope * x
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - my) ** 2)
    # If y has zero variance, define R^2 as 1.0 when residuals are zero, else 0.0
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0
    r2 = 1 - ss_res / ss_tot
    return slope, intercept, r2


def run_ablations(spectra_tpm, score_df, coexp_df, hvgs, gene_std, norm_tpm_hvg,
                       spectra_tpm_rf, baseline_usages, ntop, absolute):
    results = []
    for program in spectra_tpm.index:
        print(f"Ablating genes for program {program}")

        spectra_row = spectra_tpm.loc[program]
        score_row = score_df.loc[program]
        baseline_usage_vec = baseline_usages[program].values
        
        rows = []
        for gene, rank, score in rank_genes(score_row, ntop):
            coexp_row = coexp_df.loc[gene]
            ablated_usage, stats = ablate_gene( spectra_row, coexp_row, hvgs, gene_std, norm_tpm_hvg, absolute )
            beta, intercept, r2 = ols(baseline_usage_vec, ablated_usage)
        
            rows.append(dict(
                gene=gene,
                gep=program,
                rank=rank,
                score=score,
                spectra_value=spectra_row.loc[gene],
                beta=beta,
                intercept=intercept,
                r2=r2,
                r2_inv=1-r2,
                sum_abs_coexp=stats["sum_abs_coexp"],
                sum_clipped_weights=stats["sum_clipped_weights"],
                n_high_coexp=stats["n_high_coexp"],
            ))

        results.append(pd.DataFrame(rows))
        
    return pd.concat(results, ignore_index=True)


def plot_r2_by_rank(results_df, output_pdf_path):
    with PdfPages(output_pdf_path) as pdf:
        for program, df_program in results_df.groupby("gep"):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(df_program["rank"], df_program["r2_inv"], marker="o", linewidth=1, markersize=3)
            ax.set_xlabel("Gene rank in GEP")
            ax.set_ylabel("1 - R2 (ablated vs. true usage)")
            ax.set_title(program)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"R2-by-rank plot written to {output_pdf_path}")


def run_pipeline(args):
    print(f"Loading h5ad from {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    hvgs = load_hvgs(args.hvg)
    spectra_tpm = load_spectra(args.spectra)
    score_df = load_spectra(args.spectra_scores) if args.spectra_scores else spectra_tpm

    if args.geps:
        missing_geps = [g for g in args.geps if g not in spectra_tpm.index]
        if missing_geps:
            raise ValueError(f"Requested GEP(s) not found in spectra: {missing_geps}")
        spectra_tpm = spectra_tpm.loc[args.geps]
        score_df = score_df.loc[args.geps]

    hvgs = reconcile_genes(adata, hvgs, spectra_tpm, args.output)

    print("Computing baseline program usage")
    gene_std = compute_gene_std(adata, hvgs)
    baseline_usages, norm_tpm_hvg, spectra_tpm_rf = compute_baseline_usage(
        adata, hvgs, spectra_tpm, gene_std, args.nmf_params
    )

    coexp_genes = sorted(set(adata.var_names) & set(spectra_tpm.columns))
    coexp_df = load_or_compute_coexpression(args.coexp, adata, coexp_genes, args.output)

    print(f"Running per-gene ablation for {spectra_tpm.shape[0]} programs (ntop={args.ntop})")
    results_df = run_ablations(
        spectra_tpm, score_df, coexp_df, hvgs, gene_std, norm_tpm_hvg,
        spectra_tpm_rf, baseline_usages, args.ntop, args.coexp_absolute,
    )

    tsv_path = f"{args.output}.ablation_results.tsv"
    results_df.to_csv(tsv_path, sep="\t", index=False)
    print(f"Ablation results written to {tsv_path}")

    plot_r2_by_rank(results_df, f"{args.output}.r2_by_rank.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Co-expression-aware per-gene ablation of NMF programs."
    )
    parser.add_argument("--h5ad", required=True, help="h5ad file with TPM normalized counts")
    parser.add_argument("--hvg", required=True, help="File with HVG genes, one per line")
    parser.add_argument("--spectra", required=True, help="Consensus NMF spectra TPM (H matrix)")
    parser.add_argument("--ntop", type=int, default=250, help="Max genes to compute ablation for (default 250)")
    parser.add_argument("--nmf-params", required=True, help="cNMF run parameter yaml file")
    parser.add_argument("--coexp", default=None, help="Optional precomputed gene coexpression h5ad")
    parser.add_argument("--spectra-scores", default=None,
                         help="Optional spectra scores (SD units) used for ranking genes; defaults to --spectra")
    parser.add_argument("--geps", nargs="+", default=None,
                         help="Optional subset of GEP/program names (spectra row labels) to run ablation for; defaults to all")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--coexp-absolute", action="store_true",
                         help="Use absolute co-expression values instead of ignoring values < 0")
    args = parser.parse_args()

    run_pipeline(args)
