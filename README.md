# SUVA

## Introduce

SUVA (Splice sites Usage Variation Analysis) is a new method  that defined AS by splice site usage without prior annotations. It is powerful to analyse complex AS events.

Website： https://github.com/ablifedev/SUVA/

## Download

Download current release: [SUVA_v2.0](https://github.com/ablifedev/SUVA/archive/SUVA_v2.0.tar.gz)


## Installation

The program only works on linux platform. The detailed requirements are as followed:  

Required software:

* Perl (>=5.10, https://www.perl.org/get.html)

* Python3 (>=3, https://www.python.org/downloads/)

* R (>=3.5.0, https://cloud.r-project.org/)

Required python packages:

* HTSeq

After all of the requirments being installed, you can install the program as followed:   

1. Unzip the file `SUVA_v2.0.tar.gz`:
    `tar -zxf SUVA_v2.0.tar.gz`


## Usage

You can follow the example instruction in example data folder. Simply, after you make a config file, you can run the pipeline with the command:

    `perl suva.pl -c suva_samplepair.config -g gff_file -a geneanno_file -p fdr_threshold_for_filter -t cpu_number`

For comparision of groups with replicates, you can use the template config file: suva_rep.config

