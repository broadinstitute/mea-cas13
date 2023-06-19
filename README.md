# MEAs for Cas13 diagnostic guide design &nbsp;&middot;&nbsp; [![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE) 

The model-guided exploration algorithms (MEAs) implemented in this repository design diagnostic probes that optimize well-defined functions over viral sequence diversity. 

To learn more, please see our [preprint](blank link) that describes the algorithms implemented here and computat

## Table of Contents
* [Setting up the MEAs](#setting-up-the-meas)
* [Using the model-guided exploration algorithms](#using-the-model-guided-exploration-algorithms)
* [Output format](#output-format)
* [Example usage](#example-usage)
    * [Multi-Target Detection](#multi-target-detection)
    * [Variant Identification](#variant-identification)
* [Summary of contents](#summary-of-contents)
* [License](#license)

## Setting up the MEAs

The algorithms developed in this project depend on many Python packages. 

We've included a script, [`setup.sh`](./setup.sh) that will download the required Python dependencies, and install the mea Python package.

To setup the MEA package, we suggest cloning the repository and running the setup script as follows:
```bash
# Clone the repository
git clone https://github.com/broadinstitute/mea-cas13.git
cd mea-cas13

# Create and activate a conda environment
conda create --name mea-package-env python=3.7.10
conda activate mea-package-env

# Run the setup script
./setup.sh
```

Alternatively, if you'd like to do this manually, you can create a virtual environment and use pip to install all of the packages listed in [`requirements.txt`](./requirements.txt).

After completing this setup, you can design guides using the the `design_guides.py` script as described in the below sections.

## Using the model-guided exploration algorithms

Once the required package dependencies are installed, you can use the algorithms to design diagnostic probes.

The main program for designing diagnostic assays is [`design_guides.py`](./bin/design_guides.py). 

[`design_guides.py`](./design_guides.py) requires two four positional arguments:
```bash
python design_guides.py [objective] [exploration_algorithm] [fasta_path] [results_path] ...
```

* `objective` indicates which objective guides will be designed for. `mult` designs guides to achieve optimal activity over sequence diversity (multi-target detection objective). `diff` designs guides to optimally differentiate between a provided on-target and off-target set (variant identification objective).
* `exploration_algorithm` indicates which exploration algorithm will be used to design the diagnostic guide sequences.
`wgan-am` runs the WGAN-AM algorithm, `evolutionary` runs the evolutionary algorithm, and `both` runs both algorithms
* `fasta_path` is a path to the input .FASTA file(s) that the guides will be designed for
    * For the `mult` objective, this should be a path to a single .FASTA file. 
    * For the `diff` objective, this should be a path to a directory containing two or more .FASTA files that are all aligned to one another (e.g. position 1 in A.fasta corresponds to position 1 in B.fasta and C.fasta). The algorithms will design a guide for all the possible on-target and off-target pairs. (e.g. guides will be designed with A.fasta as the on-target and B.fasta and C.fasta as the off-targets, another set of guides will be designed with B.fasta as the on-target and A.fasta and C.fasta as the off-targets, and a set of guides  be designed with C.fasta as the on-target and A.fasta and B.fasta as the off-targets)
* `results` is a path to a folder where the results are stored. The format of the output files is described [below](#output-format).

`design_guides.py` also has several optional arguments, as detailed below:
* `--use_range` is a path to a .tsv file that specifies  the ranges of genomic sites that should be considered. The .tsv file should have to columns of integers: the first column is the start of the included range, and the second column is the end of the included range.
* `--num_cpu` tells the script how many CPUs it can parallelize jobs over. `--num_cpu` is automatically set to 1, as most machines only have 1 CPU, but if this program is being run on a machine with more CPUs, then increasing `--num_cpu` will enable multiple genomic sites to be designed for at once, enabling faster performance.

## Output format

The algorithm creates one file per explorer run, stored at `results_path/results_compiled_wgan_am.tsv/` or `results_path/results_compiled_evolutionary.tsv/` that summarizes the guide designs with the highest fitness across all the genomic sites considered

All compiled results files have the following columns:
* `algo` is the algorithm used to design the guides.
* `fitness` is the fitness of the guide computed for the design objective specified (either `mult` or `diff`).
* `guide_sequence` is the sequence of the guide that was designed by the algorithm. *NOTE: This is in the frame of the protospacer, so it must be reverse-complemented, converted from a DNA to RNA sequence, and have the LwCas13a direct repeat added on in order to serve as a functional guide RNA sequence.*
* `start_pos` indicates the positions in the target where the guide binds. More specifically, the guide's spacer binds the positions `start_pos` to `start_pos + 28` in the target.
* `shanon_entropy` is the entropy at that genomic site.
* `G_PFS` indicates whether or not the genomic site has a 'G' at the protospacer-flanking site.

The output files of guides designed for the `mult` objective have the following additional columns:
* `perc_highly_active` indicates the percentage of target in the target set that the guide is determined to be active on. In order to be considered highly active, a guide must be classified as active by the classifier model and have a predicted activity within the top quartile by the regression model.
* `shannon_entropy` is the average shannon entropy across the 28 nucleotide guide-binding region of the target set. A higher shannon entropy indicates that the genomic site has a higher degree of genomic diversity.

The output files of guides designed for the `diff` objective have the following additional columns:
* `on_target_name` is the name of the on-target set. This is taken from the name of the .fasta file that served as the input.
* `off_target_name` is the name of the off-target set. This is also taken from the name of the .fasta file that served as the input.
* `mean_on_target_act` is the mean activity of the generated guide against the on-target set. The activity is computed as described in the methods section of our manuscript.
* `mean_off_target_act` is the mean activity of the generated guide against the off-target set.

The algorithms also create additional files in the `results_path`:
* `results_path/processed_sites.pkl` is a pickled file containing information for all the genomic sites that the algorithm will design guides across.
* `results_path/results_by_site/` contains subdirectores `results_path/results_by_site/evolutionary/`, `results_path/results_by_site/wgan_am/`, or both. Each .tsv file in these subdirectories holds the results for a genomic site that the algorithms were run against. The top designs across all the genomic sites that the algorithms considered are compiled by `design_guides.py` and are available in `results_path/results_compiled_wgan_am.tsv/` and `results_path/results_compiled_evolutionary.tsv/`


## Example usage

Below is an example of how the [`design_guides.py`](./design_guides.py) script can be run on genomic data to automatically design optimal diagnostic guides.

### Multi-Target Detection

This repository includes an example alignment of parainfluenza virus 4 genomes in the file [`PIV4.fasta`](./examples/input/PIV4.fasta). To design guides based on these genomes using the evolutionary exploration algorithm, run the below script:

```bash
python design_guides.py multi evolutionary ./examples/input/PIV4.fasta ./PIV4_example_output/ --num_cpu 1 --use_range ./examples/input/PIV4_example_range.tsv
```
This will use the evolutionary algorithm to design diagnostic guides for the genomic sequences in [`PIV4.fasta`](./examples/input/PIV4.fasta) using 1 CPU. Only the genomic sites in [`ranges.tsv`](./examples/input/PIV4_example_range.tsv) or the positions between 2100 to 2200 in the allignment (corresponding to the beginning of the V gene) will considered. The results will be saved to `./PIV4_example_output/`. The runtime of the algorithms scales linearly with the number of genomic sites considered, if you are planning to run the algorithms on tens of thousands of positions, we recommend using a machine with many CPUs for the jobs to be parallelized over.

### Variant Identification

This repository includes example .fasta files for the region surrounding the E484K SNP in SARS-CoV-2 `./examples/input/E484K_variant_identification`. Specifically, [`E484.fasta`](./examples/input/E484K_variant_identification/E484.fasta) contains the WT sequence and [`E484K.fasta`](./examples/input/E484K_variant_identification/E484K.fasta) contains the sequence with the mutation.

To design guides for this variant identification task using the WGAN-AM exploration algorithm, run the below script:
 
```bash
python design_guides.py diff wgan_am ./examples/input/E484K_variant_identification/ ./E484K_example_output/ --num_cpu 1
```
This will run the WGAN-AM algorithm to design diagnostic guides that can identify the targeted SNP, and will save the per-site and compiled results to `./E484K_example_output/`.

## Summary of contents
Below is a summary of this repository's contents:
* `design_guides.py`: Main python script that should be run by users. This script, which is further described below, designs diagnostic probes using the model-based exploration algorithms.
* `requirements.txt`: List of the pip package dependencies that are required to be installed for the algorithms to run.
* `examples/input/`: Contains sample aligned FASTA files of viral genomes and .tsvs for whitelisted genomic regions that can be used to test the algorithms.
* `mea/utils/`: Scripts that enable the import of genomic sequences from FASTA files, the manipulation of DNA/RNA sequences, the training and import of the WGAN, and the usage of the predictive models.
* `mea/explorers/`: Scripts that implement both the WGAN-AM algorithm and the evolutionary algorithm explorers.
* `mea/models/`: Scripts that evaluate the fitness of guides designed for two objective functions: multi-target detection (cas13_mult) and variant identification (cas13_diff).

## License
The MEA design software is licensed under the terms of the [MIT license](./LICENSE).
