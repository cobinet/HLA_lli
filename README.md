Each folder has the pipeline for a different HLA typing software, where each folder is named after the typer used by that pipeline.

# Requirements

You need a functioning nextflow installation and nextflow-compatible container software (i.e. docker, singularity) to run this project.

Every pipeline is expecting a folder `1-Input` with the paired-end reads in the working path.

# Launch in HPC
As the pipeline is both memory and cpu intensive it's recommended to launch the nextflow process wrapped onto a sbatch, so that it does not run on the login node. for further information consult this [guide](https://lescailab.unipv.it/guides/eos_guide/use_nextflow.html). Please make sure to estimate how long is the process going to take and modify accordingly thi `launch_nf.job` file.
```bash
sbatch launch_nf.job main.nf nextflow.conf
```
Supported time formats are  "min", "min:sec", "hours:min:sec", "days-hours", "days-hours:min" and "days-hours:min:sec". Also for testing you could use:
```bash
srun --mem=2g --pty /bin/bash
```