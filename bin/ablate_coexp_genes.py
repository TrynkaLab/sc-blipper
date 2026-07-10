#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
import yaml
from tqdm import tqdm
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

#def compute_baseline_usage_dot(adata, norm_tpm_hvg, spectra_tpm_rf):
def compute_baseline_usage_dot(adata, hvgs, spectra_tpm, gene_std, nmf_params_path):

    """
    Baseline usage as the raw dot product of each cell's normalized profile
    with each program's (unablated) reweighted spectra, i.e. the same
    operation used by ablate_gene, with weights = 1. Unlike compute_baseline_usage,
    this does not jointly refit against all programs and has no non-negativity
    constraint, so it is the fair comparison point for ablate_gene's dot-product usage.
    """
    
    norm_tpm = adata[:, hvgs].copy()
    if sp.issparse(norm_tpm.X):
        sc.pp.scale(norm_tpm, zero_center=False)
        norm_tpm_hvg = np.asarray(norm_tpm.X.todense())
    else:
        norm_tpm.X = norm_tpm.X / norm_tpm.X.std(axis=0, ddof=1)
        norm_tpm_hvg = np.asarray(norm_tpm.X)
    
    spectra_tpm_rf = spectra_tpm.loc[:, hvgs].div(gene_std.loc[hvgs], axis=1)

    dot_usages = norm_tpm_hvg.dot(spectra_tpm_rf.values.T)
    baseline_usages = pd.DataFrame(dot_usages, index=adata.obs_names, columns=spectra_tpm_rf.index)

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
    n_gt_90 = int((values > 0.9).sum())
    n_gt_75 = int((values > 0.75).sum())
    n_gt_50 = int((values > 0.5).sum())

    clipped = values.clip(lower=0)
    weights = 1 - clipped

    stats = dict(
        sum_abs_coexp=sum_abs_coexp,
        sum_clipped_weights=weights.sum(),
        n_gt_90=n_gt_90,
        n_gt_75=n_gt_75,
        n_gt_50=n_gt_50,
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
                       spectra_tpm_rf, baseline_usages, ntop, absolute, geps):
    results = []
    for program in spectra_tpm.index:
        
        if geps:
            if program not in geps:
                continue
        
        spectra_row = spectra_tpm.loc[program]
        score_row = score_df.loc[program]
        baseline_usage_vec = baseline_usages[program].values

        rows = []
        ranked_genes = rank_genes(score_row, ntop)
        for gene, rank, score in tqdm(ranked_genes, desc=f"Ablating {program}"):
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
                n_gt_90=stats["n_gt_90"],
                n_gt_75=stats["n_gt_75"],
                n_gt_50=stats["n_gt_50"],
            ))

        results.append(pd.DataFrame(rows))
        
    return pd.concat(results, ignore_index=True)


def get_expression_df(adata, genes):
    sub = adata[:, genes].copy()
    sc.pp.log1p(sub)
    X = np.asarray(sub.X.todense()) if sp.issparse(sub.X) else np.asarray(sub.X)
    return pd.DataFrame(X, index=adata.obs_names, columns=genes)


def plot_true_vs_ablated_scatter(gene, ablated_usage_df, baseline_usages, output_pdf_path):
    with PdfPages(output_pdf_path) as pdf:
        for gep in ablated_usage_df.columns:
            true_usage = baseline_usages[gep].values
            ablated_usage = ablated_usage_df[gep].values
            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(true_usage, ablated_usage, s=5, alpha=0.4, edgecolor="none", rasterized=True)
            ax.set_xlabel("True GEP usage")
            ax.set_ylabel(f"Ablated GEP usage ({gene} removed)")
            ax.set_title(gep)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"True vs ablated usage scatter plots written to {output_pdf_path}")
    
            #lo = min(true_usage.min(), ablated_usage.min())
            #hi = max(true_usage.max(), ablated_usage.max())
            #ax.plot([lo, hi], [lo, hi], "k--", linewidth=1)


def plot_ablated_usage_matrix(gene, ablated_usage_df, output_pdf_path):
    order = ablated_usage_df.mean(axis=1).sort_values(ascending=False).index
    matrix = ablated_usage_df.loc[order]
    fig, ax = plt.subplots(figsize=(max(4, 0.6 * matrix.shape[1]), 8))
    im = ax.imshow(matrix.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=90)
    ax.set_yticks([])
    ax.set_xlabel("GEP")
    ax.set_ylabel("Cells")
    ax.set_title(f"Ablated GEP usage ({gene} removed)")
    fig.colorbar(im, ax=ax, label="Usage")
    fig.tight_layout()
    fig.savefig(output_pdf_path, dpi=300)
    plt.close(fig)
    print(f"Ablated usage matrix plot written to {output_pdf_path}")


def plot_coexpression_distribution(gene, coexp_row, output_pdf_path):
    values = coexp_row.drop(labels=[gene], errors="ignore")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=50, color="steelblue", edgecolor="black", linewidth=0.3)
    for threshold in (0.5, 0.75, 0.9):
        ax.axvline(threshold, color="red", linestyle="--", linewidth=1)
    ax.set_xlabel(f"Co-expression with {gene}")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"Co-expression distribution for {gene}")
    fig.tight_layout()
    fig.savefig(output_pdf_path)
    plt.close(fig)
    print(f"Co-expression distribution plot written to {output_pdf_path}")


def plot_top_partner_scatter(gene, coexp_row, adata, ntop, output_pdf_path):
    ranked = coexp_row.drop(labels=[gene], errors="ignore").abs().sort_values(ascending=False)
    partners = ranked.head(ntop).index.tolist()
    expr = get_expression_df(adata, [gene] + partners)

    fig, axes = plt.subplots(1, len(partners), figsize=(4 * len(partners), 4), squeeze=False)
    for ax, partner in zip(axes[0], partners):
        ax.scatter(expr[gene], expr[partner], s=5, alpha=0.4, edgecolor="none", rasterized=True)
        ax.set_xlabel(f"{gene} (log1p)")
        ax.set_ylabel(f"{partner} (log1p)")
        ax.set_title(f"r = {coexp_row[partner]:.2f}")
    fig.tight_layout()
    fig.savefig(output_pdf_path, dpi=300)
    plt.close(fig)
    print(f"Top co-expressed partner scatter plots written to {output_pdf_path}")


def run_gene_diagnostics(gene, geps, spectra_tpm, coexp_df, hvgs, gene_std, norm_tpm_hvg,
                          baseline_usages, adata, absolute, output_prefix):
    if gene not in coexp_df.index:
        raise ValueError(f"Requested gene '{gene}' not found in coexpression matrix")

    coexp_row = coexp_df.loc[gene]

    ablated_usages = {}
    for gep in tqdm(geps, desc=f"Ablating {gene}"):
        spectra_row = spectra_tpm.loc[gep]
        ablated_usage, _ = ablate_gene(spectra_row, coexp_row, hvgs, gene_std, norm_tpm_hvg, absolute)
        ablated_usages[gep] = ablated_usage
    ablated_usage_df = pd.DataFrame(ablated_usages, index=adata.obs_names)

    matrix_tsv_path = f"{output_prefix}.{gene}.ablated_usage_matrix.tsv"
    ablated_usage_df.to_csv(matrix_tsv_path, sep="\t")
    print(f"Ablated usage matrix written to {matrix_tsv_path}")

    plot_true_vs_ablated_scatter(
        gene, ablated_usage_df, baseline_usages, f"{output_prefix}.{gene}.true_vs_ablated_scatter.pdf"
    )
    plot_coexpression_distribution(
        gene, coexp_row, f"{output_prefix}.{gene}.coexpression_distribution.pdf"
    )
    plot_top_partner_scatter(
        gene, coexp_row, adata, 5, f"{output_prefix}.{gene}.top5_coexpressed_partners.pdf"
    )


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
   
    # Load spectra
    spectra_tpm = load_spectra(args.spectra)
    spectra_tpm.index = spectra_tpm.index.astype(str)
    
    score_df = load_spectra(args.spectra_scores) if args.spectra_scores else spectra_tpm
    score_df.index = score_df.index.astype(str)
    
    print(f"Loaded spectra TPM with geps {spectra_tpm.index} from {args.spectra}")
    
    print(f"Loading h5ad from {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    hvgs = load_hvgs(args.hvg)
    
    hvgs = reconcile_genes(adata, hvgs, spectra_tpm, args.output)

    print("Computing baseline program usage")
    gene_std = compute_gene_std(adata, hvgs)
    baseline_usages, norm_tpm_hvg, spectra_tpm_rf = compute_baseline_usage_dot(
        adata, hvgs, spectra_tpm, gene_std, args.nmf_params
    )

    coexp_genes = sorted(set(adata.var_names) & set(spectra_tpm.columns))
    coexp_df = load_or_compute_coexpression(args.coexp, adata, coexp_genes, args.output)

    if args.geps:
        missing_geps = [g for g in args.geps if g not in spectra_tpm.index]
        if missing_geps:
            raise ValueError(f"Requested GEP(s) not found in spectra: {missing_geps}")
        ngeps = len(args.geps)
    else:
        ngeps = spectra_tpm.shape[0]


    print(f"Running per-gene ablation for {ngeps} programs (ntop={args.ntop})")
    if not args.gene:
        results_df = run_ablations(
            spectra_tpm, score_df, coexp_df, hvgs, gene_std, norm_tpm_hvg,
            spectra_tpm_rf, baseline_usages, args.ntop, args.coexp_absolute, args.geps
        )

        tsv_path = f"{args.output}.ablation_results.tsv"
        results_df.to_csv(tsv_path, sep="\t", index=False)
        print(f"Ablation results written to {tsv_path}")

        plot_r2_by_rank(results_df, f"{args.output}.r2_by_rank.pdf")

    if args.gene:
        gene_geps = args.geps if args.geps else spectra_tpm.index.tolist()
        print(f"Running detailed ablation diagnostics for gene '{args.gene}' across {len(gene_geps)} GEP(s)")
        run_gene_diagnostics(
            args.gene, gene_geps, spectra_tpm, coexp_df, hvgs, gene_std, norm_tpm_hvg,
            baseline_usages, adata, args.coexp_absolute, args.output
        )


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
    parser.add_argument("--gene", default=None,
                         help="Optional single gene to produce detailed ablation diagnostics for "
                              "(true vs. ablated usage scatter, ablated usage matrix, co-expression "
                              "distribution, top-5 co-expressed partner scatter). Interacts with --geps.")
    parser.add_argument("--output", required=True, help="Output prefix")
    parser.add_argument("--coexp-absolute", action="store_true",
                         help="Use absolute co-expression values instead of ignoring values < 0")
    args = parser.parse_args()

    run_pipeline(args)
