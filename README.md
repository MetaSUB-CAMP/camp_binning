# CAMP Binning


[![Documentation Status](https://img.shields.io/readthedocs/camp_binning)](https://camp-documentation.readthedocs.io/en/latest/binning.html) ![Version](https://img.shields.io/badge/version-0.9.0-brightgreen)

## Overview

This module is designed to function as both a standalone MAG binning pipeline as well as a component of the larger CAMP metagenome analysis pipeline. As such, it is both self-contained (ex. instructions included for the setup of a versioned environment, etc.), and seamlessly compatible with other CAMP modules (ex. ingests and spawns standardized input/output config files, etc.). 

As far as the binning procedure goes, the design philosophy is just to replicate the functionality of [MetaWRAP](https://github.com/bxlab/metaWRAP)(one of the original ensemble methods) with i) better dependency conflict management and ii) improved integration with new binning algorithms. 

Currently, the binning algorithms MetaBAT2, CONCOCT, SemiBin, and MaxBin2 are wrapped along with the bin refinement tool DAS Tool.

## Installation

1. Clone repo from [Github](<https://github.com/MetaSUB-CAMP/camp_binning>).
```Bash
git clone https://github.com/MetaSUB-CAMP/camp_binning
```

2. Set up the conda environment using `configs/conda/binning.yaml`. 
```Bash
# Create and activate conda environment 
cd camp_binning
conda env create -f configs/conda/binning.yaml
conda activate binning
```

3. The conda version of MaxBin2 doesn't seem to work, so the best way to add it to the module is to install it separately. 
```Bash    
cd bin/
wget https://sourceforge.net/projects/maxbin2/files/latest/download
tar -xf download
spack load gcc@6.3.0 # This is only necessary for HPCs with extremely old gcc's 
cd MaxBin-2.2.7/src
make
./autobuild_auxiliary
wget https://github.com/loneknightpy/idba/releases/download/1.1.3/idba-1.1.3.tar.gz
tar -xf idba-1.1.3.tar.gz
cd idba-1.1.3/
./configure --prefix=/path/to/bin/MaxBin-2.2.7/auxiliary/idba-1.1.3 # IDBA-UD was not included in the auxiliary build
make
# Optional: Export or add the following to ~/.bashrc
export PATH=$PATH:/path/to/bin/MaxBin-2.2.7:/path/to/bin/MaxBin-2.2.7/auxiliary/FragGeneScan_1.30:/path/to/bin/MaxBin-2.2.7/auxiliary/hmmer-3.1b1/src:/path/to/bin/MaxBin-2.2.7/auxiliary/bowtie2-2.2.3:/path/to/bin/MaxBin-2.2.7/auxiliary/idba-1.1.3/bin
```

4. Update the relevant parameters (ex. `ext/` (the location of external tools and scripts) and `maxbin2_script` (the location of the MaxBin2 script)) in `test_data/parameters.yaml`.

5. Make sure the installed pipeline works correctly. With 40 threads and a maximum of 100 GB allocated, the test dataset should finish in approximately 35 minutes.
```Bash
# Run tests on the included sample dataset
python /path/to/camp_binning/workflow/binning.py test
```

## Using the Module

**Input**: `/path/to/samples.csv` provided by the user.

**Output**: 1) An output config file summarizing the locations of 2) the MAGs generated by MetaBAT2, CONCOCT, and VAMB. See `test_data/test_out.tar.gz` for a sample output work directory.

- `/path/to/work/dir/binning/final_reports/samples.csv` for ingestion by the next module (ex. quality-checking)
- `/path/to/work/dir/binning/*/sample_name/bins/`, where `*` refers to a binning or refinement tool, which contains FastAs of MAGs inferred by that tool

### Module Structure
```
└── workflow
    ├── Snakefile
    ├── binning.py
    ├── utils.py
    ├── __init__.py
    └── ext/
        └── scripts/
```
- `workflow/binning.py`: Click-based CLI that wraps the `snakemake` and other commands for clean management of parameters, resources, and environment variables.
- `workflow/Snakefile`: The `snakemake` pipeline. 
- `workflow/utils.py`: Sample ingestion and work directory setup functions, and other utility functions used in the pipeline and the CLI.
- `ext/`: External programs, scripts, and small auxiliary files that are not conda-compatible but used in the workflow.

### Running the Workflow

1. Make your own `samples.csv` based on the template in `configs/samples.csv`.
    - `ingest_samples` in `workflow/utils.py` expects Illumina reads in FastQ (may be gzipped) form 
    - `samples.csv` requires either absolute paths or paths relative to the directory that the module is being run in

2. Update the relevant parameters in `configs/parameters.yaml`.

3. Update the computational resources available to the pipeline in `configs/resources.yaml`. 

#### Command Line Deployment

To run CAMP on the command line, use the following, where `/path/to/work/dir` is replaced with the absolute path of your chosen working directory, and `/path/to/samples.csv` is replaced with your copy of `samples.csv`. 
    - The default number of cores available to Snakemake is 1 which is enough for test data, but should probably be adjusted to 10+ for a real dataset.
    - Relative or absolute paths to the Snakefile and/or the working directory (if you're running elsewhere) are accepted!
    - The parameters and resource config YAMLs can also be customized.
```Bash
python /path/to/camp_binning/workflow/binning.py \
    (-c number_of_cores_allocated) \
    (-p /path/to/parameters.yaml) \
    (-r /path/to/resources.yaml) \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
```

#### Slurm Cluster Deployment

To run CAMP on a job submission cluster (for now, only Slurm is supported), use the following.
    - `--slurm` is an optional flag that submits all rules in the Snakemake pipeline as `sbatch` jobs. 
    - In Slurm mode, the `-c` flag refers to the maximum number of `sbatch` jobs submitted in parallel, **not** the pool of cores available to run the jobs. Each job will request the number of cores specified by threads in `configs/resources/slurm.yaml`.
```Bash
sbatch -J jobname -o jobname.log << "EOF"
#!/bin/bash
python /path/to/camp_binning/workflow/binning.py --slurm \
    (-c max_number_of_parallel_jobs_submitted) \
    (-p /path/to/parameters.yaml) \
    (-r /path/to/resources.yaml) \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
EOF
```

#### Finishing Up

1. After checking over `/path/to/work/dir/binning/*/sample_name/bins/` and making sure you have bins that seem reasonable given the input dataset, you can delete all intermediate files to save space. 
```Bash
python /path/to/camp_binning/workflow/binning.py cleanup \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
```

2. If for some reason the module keeps failing, CAMP can print a script containing all of the remaining commands that can be run manually. 
```Bash
python /path/to/camp_binning/workflow/binning.py --dry_run \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
```

## Credits

- This package was created with [Cookiecutter](https://github.com/cookiecutter/cookiecutter>) as a simplified version of the [project template](https://github.com/audreyr/cookiecutter-pypackage>).
- Free software: MIT
- Documentation: https://camp-documentation.readthedocs.io/en/latest/binning.html


