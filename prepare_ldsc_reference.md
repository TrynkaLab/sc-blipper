

# Download the reference files from Zenodo

## For build 37:
```
mkdir GRCh37

cd GRCh37

# LD scores
wget https://zenodo.org/records/10515792/files/1000G_Phase3_baselineLD_v2.2_ldscores.zip
unzip 1000G_Phase3_baselineLD_v2.2_ldscores.zip

# Weights
wget https://zenodo.org/records/10515792/files/1000G_Phase3_weights_hm3_no_MHC.tgz
tar -zxvf 1000G_Phase3_weights_hm3_no_MHC.tgz
mv 1000G_Phase3_weights_hm3_no_MHC weights

### Allele frequencies
wget https://zenodo.org/records/10515792/files/1000G_Phase3_frq.tgz
tar -zxvf 1000G_Phase3_frq.tgz
mv 1000G_Phase3_frq frq

# Plink reference
wget https://zenodo.org/records/10515792/files/1000G_Phase3_plinkfiles.tgz
tar -zxvf 1000G_Phase3_plinkfiles.tgz

# Perpare hapmap 3
mkdir hapmap3
wget https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt
mv hm3_no_MHC.list.txt hm3_no_MHC.list

# Create the hm3 list for munge sumstats
python make_w_hm3.py --rsids hm3_no_MHC.list --bim plink/1000G.EUR.hg38.*.bim
mv w_hm3.snplist hapmap3/

# Sumstats (optional)
#wget https://zenodo.org/records/10515792/files/sumstats_indep107.tgz
#tar -zxvf sumstats_indep107.tgz
```

## For build 38
```
# Reference bundle
wget https://zenodo.org/records/10515792/files/GRCh38.tgz
tar -zxvf GRCh38.tgz
cd GRCh38

# Hapmap snp
wget https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt
mkdir hapmap3
mv hm3_no_MHC.list.txt hm3_no_MHC.list

# Create the hm3 list for munge sumstats (script in bin folder)
python make_w_hm3.py --rsids hm3_no_MHC.list --bim plink/1000G.EUR.hg38.*.bim
mv w_hm3.snplist hapmap3/

# Unzip the rest
tar -xzvf weights.tgz 
tar -xzvf plink_files.tgz 
tar -xzvf baselineLD_v2.2.tgz

cd ..

# Allele freqeuncies, could be computed from PLINK
wget https://zenodo.org/records/10515792/files/1000G_Phase3_frq.tgz
tar -zxvf 1000G_Phase3_frq.tgz
mv 1000G_Phase3_frq frq
```


# Install instuctions
In future, will see if this can be integrated this into the main enviroment.yml
```
conda create -n sc-blipper-ldsc python=3.9 bedtools
conda activate sc-blipper-ldsc

pip install git+https://github.com/TrynkaLab/ldsc

# This is my fork of the CBIIT version with the follwing fixed applied
# Fix 1: There was an issue with reading annots
# Fix 2: make_annot.py in some cases returns more then the number of SNPs in the bim, possibly related to multiallelic sites. Added a deduplication so the annot always has unique SNP positions that match the bim.

# This is the original, it is has a couple of small bugs that prevent it from running
#pip install git+https://github.com/CBIIT/ldsc.git

```