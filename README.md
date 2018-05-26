# ReadHog
Recruitment fidelity of metagenome reads against a set of assembled genomes

Thanks to Sean McAllister for giving me the idea to write this, and providing feedback.

## Dependencies

-bowtie2

-SAMtools

Let's say you are analyzing a metagenome; you have already assembled you contigs, and binned them into discreet metagenome-assembled genomes (MAGs). Not you would like to know the relative abundance of each MAG, and so you map the metagenome reads to your MAGs. Sounds simple enough, right? Most will probably map reads to all the MAGs simultaneously, using a mapping software like Bowtie2 or BWA. These mapping programs employ scorring strategies that figure out which part of the assembled DNA sequences are the best match to a raw read. In some cases, for example, when a microbial community includes a set of organisms that are highly-related to each other, there may be more than one bin that recruits a read with the same score. I dont know how these mapping programs decide which bin will get that read (e.g. is it random? is it based on some arbitrary order of bins?). So Sean and I decided to create this script which is essentially a bowtie2 wrapper. This program will map the raw metagenome reads that you provide to your bins under several conditions: 1) it will map reads to all bins simultaneously, 2) it will map reads to all bins individually, and 3) to various combinations of bins (provided by the user in a TSV file). The third condition is optional.

This program accepts a folder of raw metagenome sequence reads and a folder of assembled genomes. Right now, the script only accepts unpaired reads, which all have to be in one folder; these can be single-end sequencing results, or paired-end reads that have already been combined using a program like FLASH.

The script will map reads to all the assembled genomes at the same time, and separately to each genome.
The script will also take as input a tab-delimited file in which every row contains a set of genomes to which reads will be mapped to at the same time; this is optional.

The output will consist of a CSV file in which the number of reads recruited to each genome under each mapping condition is reported as separate columns.

### Example:

    
