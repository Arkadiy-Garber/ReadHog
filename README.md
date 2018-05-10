# ReadHog
Recruitment fidelity of metagenome reads against a set of assembled genomes

## Dependencies

-bowtie2

-SAMtools

This program accepts a folder of raw metagenome sequence reads and a folder of assembled genomes.
The script will map reads to all the assembled genomes at the same time, and to spearately to each genome.
The script will also take as input a tab-delimited file in which every row contains a set of genomes to which reads will be mapped to at the same time. This is optional.

The output will consist of a CSV file in which the number of reads recruited to each genome under each condition is reported as separate columns.
