Each folder has the pipeline for a different HLA typing software, where each folder is named after the typer used by that pipeline.

# Requirements

You need a functioning nextflow installation and nextflow-compatible container software (i.e. docker, singularity) to run this project.

Every pipeline is expecting a folder `1-Input` with the paired-end reads in the working path.
