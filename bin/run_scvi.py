#!/usr/bin/env python3
import logging
import torch
import jax
import jax.numpy as jnp
import scanpy as sc
import scvi
import anndata as ad
from scipy.sparse import hstack, issparse
import numpy as np
import argparse
import matplotlib.pyplot as plt
import random
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_pytorch_cuda():
    has_cuda=False
    cuda_available = torch.cuda.is_available()
    logging.info(f"PyTorch CUDA available: {cuda_available}")
    if cuda_available:
        has_cuda=True
        device_count = torch.cuda.device_count()
        logging.info(f"Number of CUDA devices: {device_count}")
        for i in range(device_count):
            logging.info(f"Device {i}: {torch.cuda.get_device_name(i)}")
    else:
        logging.warning("No CUDA devices available in PyTorch.")
        
    return has_cuda

def check_jax_devices():
    devices = jax.devices()
    has_jax=False
    logging.info(f"JAX devices detected: {devices}")
    # Optionally, show device types separately
    cpu_devices = jax.devices("cpu")
    logging.info(f"JAX CPU devices: {cpu_devices}")
    try:
        cuda_devices = jax.devices("cuda")
        logging.info(f"JAX CUDA devices: {cuda_devices}")
        has_jax=True
    except RuntimeError:
        logging.warning("No CUDA devices found by JAX.")
        
    return has_jax


def stdscale_quantile_celing(_adata, max_value=None, quantile_thresh=None):
    sc.pp.scale(_adata, zero_center=False, max_value=max_value)    
    if quantile_thresh is not None:
        if issparse(_adata.X):
            threshval = np.quantile(np.array(_adata.X.todense()).reshape(-1), quantile_thresh)
        else:
            threshval = np.quantile(_adata.X.reshape(-1), quantile_thresh)
            
        _adata.X[_adata.X > threshval] = threshval     
        
        
def main(args):
    
    #-------------------------------------------------------------------------
    # Partially reverse engineered from https://github.com/dylkot/cNMF/blob/main/src/cnmf/preprocess.py
    # The normal workflow for cNMF preprocssing (no Harmony) is:
    # 1. cNMF preprocessing: raw counts  > TP10K (all genes), HVG list selection
    # 2. cNMF preprocessing: raw counts[,HVG] > stdscale_quantile_celing * > harmony (or not)
    # 3. cNMF: get_norm_counts > re-scales to unit variance (if running cNMF Prepare, this has already beed one in step 2)
    # *(scales raw counts to unit variance, but doesn't center)
    #
    # Here
    # 1. cNMF preprocessing: TP10K + HVG selection (identical to cNMF without Harmony)
    # 2. Prepare the raw counts adata for scVI training (not unit normalized)
    # 3. Run scVI
    # 4. Save outputs applying stdscale_quantile_celing (in cNMF prepare this is )
    
    
    if args.seed is not None:
        seed = int(args.seed)
        scvi.settings.seed = seed
        os.environ["PYTHONHASHSEED"] = str(seed)
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

        # Optional: make PyTorch more deterministic
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    
    # Step 1. cNMF preprocessing
    # Read the anndata
    adata = ad.read_h5ad(args.input)
    
    #-------------------------------------------------------------------------
    missing_cols = []
    
    # Check batch_key
    if args.batch_key is not None and args.batch_key not in adata.obs.columns:
        raise ValueError(f"batch_key not found on .obs: '{args.batch_key}'")
    
    # Check categorical covariates
    if args.cat_covariates is not None:
        for col in args.cat_covariates:
            if col not in adata.obs.columns:
                missing_cols.append(f"categorical covariate: '{col}'")
    
    # Check continuous covariates
    if args.cont_covariates is not None:
        for col in args.cont_covariates:
            if col not in adata.obs.columns:
                missing_cols.append(f"continuous covariate: '{col}'")
    
    if missing_cols:
        error_msg = "The following columns are missing from adata.obs:\n"
        error_msg += "\n".join([f"  - {col}" for col in missing_cols])
        error_msg += f"\nAvailable columns: {list(adata.obs.columns)}"
        raise ValueError(error_msg)
    
    logging.info("All required columns validated successfully in adata.obs")
    
    #-------------------------------------------------------------------------
    # tp10k matrix
    tp10k = adata.copy()
    sc.pp.normalize_total(tp10k, target_sum=1e4, copy=False)
    
    # HVGs
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=args.variable_genes)

    # Subset counts to HVGs
    adata = adata[:, adata.var['highly_variable']]
    hvgs = list(adata.var.index)
    #-------------------------------------------------------------------------
    # Step 2. Prepare adata for scVI training
    adata = adata.copy()
    
    if (args.model is None):
        scvi.model.SCVI.setup_anndata(
            adata,
            batch_key=args.batch_key,
            categorical_covariate_keys=args.cat_covariates,
            continuous_covariate_keys=args.cont_covariates,
        )

        model = scvi.model.SCVI(adata, n_latent=args.n_latent)
        
        # Step 3. Run scVI
        model.train(
            max_epochs=args.epochs,
            validation_size          = 0.1,
            early_stopping           = args.skip_early_stopping,
            early_stopping_monitor   = "elbo_validation",
            early_stopping_warmup_epochs=25)
        
        model.save(f"{args.output_prefix}", overwrite=True)
    else:
        model = scvi.model.SCVI.load(args.model, adata=adata)
    
    # Step 4. Get the output, library_size set to latent (returns on the scale of each cell, as opposed to a fixed tpm number like 1e4)
    # The library size is arbitrary here
    # TODO: in future, remove the copy step
    denoised = adata.copy()
    denoised.X = model.get_normalized_expression(denoised, library_size="latent")
     
    # Add the latent embedding    
    denoised.obsm['scVI_latentspace'] = model.get_latent_representation(denoised)
    
    # Save the denoised data before scaling
    if not args.skip_scvi_h5ad:
        sc.write(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '.h5ad', denoised)
    

    if not args.skip_cnmf_h5ad:
        # Apply the same scaling as in cNMF preprocessing
        stdscale_quantile_celing(denoised, max_value=args.max_scaled_thresh, quantile_thresh=args.quantile_thresh)
        
        # Save the ouputs
        sc.write(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '.Corrected.HVG.Varnorm.h5ad', denoised)
        sc.write(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '.TP10K.h5ad', tp10k)
        with open(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '.Corrected.HVGs.txt', 'w') as F:
            F.write('\n'.join(hvgs))

    #-----------------------------------------------------------------------------------
    if not args.skip_plot:
        
        logging.info("Generating and saving ELBO training plot")
        elbo_train = model.history["elbo_train"]
        plt.plot(elbo_train, label="train")
        plt.legend()
        plt.tight_layout()
        plt.savefig(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '_elbo_train.png',
                    bbox_inches='tight', dpi=300)
        plt.close()
            
        logging.info("Generating and saving UMAP plots for TP10K and scVI denoised data")
        
        fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))

        # -------------------- add TP10K HVG UMAP  --------------------
        # Subset tp10k to the previously computed HVGs, log1p, run PCA + UMAP and save plot
        try:
            if len(hvgs) > 0:
                tp10k_hvg = tp10k[:, hvgs].copy()
                # normalization and dimensionality reduction
                sc.pp.log1p(tp10k_hvg)
                sc.pp.pca(tp10k_hvg, n_comps=50)
                sc.pp.neighbors(tp10k_hvg, use_rep="X_pca")
                sc.tl.umap(tp10k_hvg, min_dist=0.3)

                sc.pl.umap(
                    tp10k_hvg,
                    color=args.batch_key,
                    frameon=False,
                    show=False,
                    ax=axes[0],
                    title="TP10K HVG UMAP")
                
            else:
                logging.warning("HVG list is empty; skipping TP10K HVG UMAP.")
        except Exception as e:
            logging.warning(f"Failed to generate TP10K HVG UMAP: {e}")
        
        logging.info("TP10K HVG UMAP generated successfully.")
        # -------------------- add scVI HVG UMAP  --------------------
        sc.pp.neighbors(denoised, use_rep="scVI_latentspace")
        sc.tl.umap(denoised, min_dist=0.3)
        sc.pl.umap(
            denoised,
            color=args.batch_key,
            frameon=False,
            show=False,
            ax=axes[1],
            title="scVI HVG UMAP"
        )
      
        plt.tight_layout()
        plt.savefig(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '_prepost_umap.png',
                    bbox_inches='tight', dpi=300)
        plt.close()
        
        # Save SCVI and TP10K UMAP coordinates and cell info
        umap_df = denoised.obs.copy()
        umap_df[['scVI_UMAP_1', 'scVI_UMAP_2']] = denoised.obsm['X_umap']
        umap_df[['tp10k_UMAP_1', 'tp10k_UMAP_2']] = tp10k_hvg.obsm['X_umap']
        umap_df.to_csv(args.output_prefix + "__scvi__nl" + str(args.n_latent) + '_meta_with_umap_coords.tsv', sep='\t')
        
        
        logging.info("UMAP plots saved successfully. Done.")
    else:
        logging.info("Plot generation skipped as per user request. Done.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scVI training on single-cell RNA-seq data with cNMF preprocessing")
    parser.add_argument("--input", type=str, required=True, help="Path to input h5ad file with raw counts on .X")
    parser.add_argument("--output_prefix", type=str, required=True, help="Output prefix for saved files")
    parser.add_argument("--batch_key", type=str, required=True, help="Column name in adata.obs for batch information")
    parser.add_argument("--cat_covariates", type=str, nargs="+", default=None, help="Categorical covariate column names in adata.obs")
    parser.add_argument("--cont_covariates", type=str, nargs="+", default=None, help="Continuous covariate column names in adata.obs")
    parser.add_argument("--variable_genes", type=int, default=2000, help="Number of highly variable genes to select [default: 2000]")
    parser.add_argument("--n_latent", type=int, default=10, help="Dimensionality of latent space [default: 10]")
    parser.add_argument("--max_scaled_thresh", type=float, default=None, help="Maximum value for scaling [default: None]")
    parser.add_argument("--quantile_thresh", type=float, default=0.9999, help="Quantile threshold for ceiling [default: 0.9999]")
    parser.add_argument("--model", type=str, default=None, help="Path to pre-trained scVI model. If provided, skips training [default: None]")
    parser.add_argument("--epochs", type=int, default=800, help="Number of training epochs for scVI [default: 800]")
    parser.add_argument("--skip_early_stopping", default=False, action='store_true', help="Whether to skip early stopping during training [default: False]")
    parser.add_argument("--skip_plot", default=False, action='store_true', help="Whether to generate UMAPs and save plots [default: False]")
    parser.add_argument("--skip_scvi_h5ad", default=False, action='store_true', help="Whether to skip saving the unscaled scVI h5ad file [default: False]")
    parser.add_argument("--skip_cnmf_h5ad", default=False, action='store_true', help="Whether to skip saving the unscaled cNMF ready output files [default: False]")
    parser.add_argument("--seed", type=int, default=None, help="Seed for reproducible results [default: None (random)]")

    args = parser.parse_args()
    
    if not check_jax_devices() or not check_pytorch_cuda():
        logging.error("CUDA devices are not properly configured for PyTorch or JAX. Check your installation or GPU drivers.")
        exit(1)
    
    #if args.cat_covariates is not None:
    #    args.cat_covariates = list(args.cat_covariates.split(','))
    
    #if args.cont_covariates is not None:
    #    args.cont_covariates = list(args.cont_covariates.split(','))
    
    main(args)
    
    #args = argparse.Namespace()
    #args.input = "/lustre/scratch124/humgen/projects_v2/healthy_imm_expr/sc_blipper_testing/data/h5ad/merged/testing_cnmf/testing_cnmf_merged.h5ad"
    #args.batch_key = "orig_h5ad"
    #args.cat_covariates = ["donor_id_bioivt"]
    #args.cont_covariates = None
    #args.variable_genes = 2000
    #args.output_prefix = "/lustre/scratch124/humgen/projects_v2/healthy_imm_expr/sc_blipper_testing/testing_scvi/scvi_testing"
    #args.n_latent = 10
    #args.max_scaled_thresh = None
    #args.quantile_thresh = 0.9999
    #args.model=None
    #args.skip_plot=False
    #args.skip_scvi_h5ad=False
