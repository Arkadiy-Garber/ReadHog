# ReadHog
Recruitment fidelity of metagenome reads against a set of assembled genomes

Thanks to Sean McAllister for giving me the idea to write this, and providing feedback.

## Dependencies

-bowtie2

-SAMtools

You have assembled you contigs, and binned them into metagenome-assembled genomes (MAGs). Now, you would like to know the relative abundance of each MAG, or measure the expression of certain genes, and so you map the metagenome or metatranscriptome reads to your MAGs. Standard procedure emplyoed by most researcher is to map reads to all MAGs simultaneously, using a mapping software like Bowtie2 or BWA. These mapping programs employ scorring strategies that figure out which part of the assembled DNA sequences are the best match to a raw read. In some cases, for example, when a microbial community includes organisms that are highly-related to each other (for example, different strain of the same species), there may be more than one bin that recruits a read with the same score. This may lead to innaccuracies with regard to haw many reads are recruited by particular MAG or gene.

To address these potential issues ReadHog will: 1) map reads to all bins simultaneously, 2) map reads to all bins individually, and 3) map reads to various combinations of bins (provided by the user in a TSV file). The third condition is optional.

This program accepts a folder of raw metagenome sequence reads (in fastq format) and a folder of assembled genomes. Right now, the script only accepts unpaired reads, which all have to be in one folder; these can be single-end sequencing results, or paired-end reads with overlap that have  been combined using a program like FLASH.

The script will map reads to all the assembled genomes simulataneously, and separately to each genome.
The script will also take as input a tab-delimited file (see samples folder for the required format of this file) in which every row contains a set of genomes to which reads will be mapped to at the same time..

The output will consist of a CSV file in which the number of reads recruited to each genome under each mapping condition is reported as separate columns.

### Example:

    python3 ReadHog.py -reads_directory reads/ -reads_ext fastq -bin_directory bins/ -bin_ext fasta -output_csv readMappingInfo.csv -bin_comparisons binComparisons.txt -num_threads 8
