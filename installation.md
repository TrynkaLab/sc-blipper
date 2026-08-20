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


# Verify GPU dependencies

The shared `environment.yml` includes the CPU, CUDA cNMF, and scVI dependencies.
If you plan to run GPU cNMF or scVI, test the environment on a GPU node:

```python
import jax
import scvi
import torch

print(torch.cuda.is_available())  # Should return True
print(jax.devices())              # Should include a GPU device
print(scvi.__version__)
```
