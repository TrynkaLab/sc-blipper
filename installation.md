# Create conda env

Create and activate conda environment
```bash
conda env create -f environment.yml
conda activate sc-blipper
```

> Note: If you want to give the env a different name: `conda env create -f environment.yml -n my_env`

Note down the install path for the conda env, we will need it later for configuring the pipeline
```bash
echo $CONDA_PREFIX  # On Unix/Linux/macOS
```

Exit conda env
```bash
conda deactivate
```

# I cannot use conda, now what?
You may be able to create a singularity/apptainer container based on the conda enviroment, As of writing I have not tested
this. Some resources below:
- https://arcdocs.leeds.ac.uk/arc3-arc4/usage/conda-containers.html
- https://stackoverflow.com/questions/76146763/create-apptainer-container-with-environment-yml-without-creating-a-new-conda-env
- https://csc-training.github.io/csc-env-eff/hands-on/singularity/singularity_extra_replicating-conda.html

Once you have your containers, you can set the parameters 

```
params.rn_container="/path/to/container"
params.scvi.container="/path/to/container-scvi"
```

Also make sure to tell nextflow to use your container engine of choice, example for singularity:
```
conda.enabled = false
singularity.enabled = true
```

More details: https://www.nextflow.io/docs/latest/container.html#singularity


# Install scVI (optional)
We will make a new enviroment for scVI to avoid conflicts and keep scVI optional so the pipeline is lighter
You will ideally need a GPU for scVI to run
```bash
conda env create -f environment_scvi.yml
conda activate sc-blipper-scvi 
```

Note down the install path for the conda env, we will need it later for configuring the pipeline
```bash
echo $CONDA_PREFIX  # On Unix/Linux/macOS
```

At this point I strongly recommend testing GPU is working. If on HPC, make sure to request a GPU node
You can test by running python and typing:
```python
import torch
print(torch.cuda.is_available()) # Should return True
import jax
print(jax.devices()) # Should show GPU devices
```

Exit conda env
```bash
conda deactivate
```

## If this doesn't work

Remove the conda enviroment you just created
```bash
conda remove -n sc-blipper-scvi --all
```

Manually create a new one
```bash
conda create -n sc-blipper-scvi python==3.10
conda activate sc-blipper-scvi
```

Follow instructions here: https://docs.scvi-tools.org/en/1.0.0/installation.html
Install Pytorch, JAX, they might need to be tweaked depending on your system and GPU
Pytorch: https://pytorch.org/get-started/locally/
JAX: https://docs.jax.dev/en/latest/installation.html#installation

Cuda 12 on linux, make sure the CUDA version matches your GPU drivers
```bash
pip3 install torch torchvision
pip3 install -U "jax[cuda12]"
pip install scvi-tools scikit-misc ipython
```

