#!/usr/bin/env python3
# !/bin/sh
# Author: Arkadiy Garber
from collections import defaultdict
import os
import re
import textwrap
import argparse


parser = argparse.ArgumentParser(
    prog="readhog.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    Developed by Arkadiy Garber^1 and Sean M. McAllister^1;
    ^1University of Delaware, Geological Sciences
    Please send comments and inquiries to arkg@udel.edu
    '''))

parser.add_argument('-reads_directory', help='directory with read files in fastq format (required)')
parser.add_argument('-reads_ext', help="the extention used for each fastq read file (do not include the period) "
                                     "(default = fastq)", default="fastq")
parser.add_argument('-bin_directory', help="directory with bins (required)")
parser.add_argument('-bin_ext', help="the extention used for each bin file (do not include the period) "
                                    "(default = fasta)", default="fasta")
parser.add_argument('-output_csv', help="basename for output file with read mapping information"
                                 " (default = read_mapping_info)", default="read_mapping_info.csv")
parser.add_argument('-bin_comparisons', help="tab-delimited file for which combination of bins to map reads to "
                                             "(optional)", default="NA")
parser.add_argument('-num_threads', help="the number of threads to use (default = 1)", default=1)
parser.add_argument('-mode', help="paired (if you are providing a directory with paired-end reads) or single (if your folder"
                                  " has only single-end reads (default = single). If you are providing a directory with "
                                  "paired-end reads, all forward-end read file names must end with either _1.fq or _1.fastq; "
                                  " all reverse-end file names must end with _2.fq of _2.fastq. All other files in the reads directory that "
                                  "do not end with _1.fq, _1.fastq, _2.fq, or _2.fastq will be assumed to be unpaired "
                                  "reads", default="single")

args = parser.parse_args()


def sum(list):
    count = 0
    for i in list:
        count += int(i)
    return count


def lastItem(iterable):
    if len(iterable) > 1:
        x = ''
        for i in iterable:
            x = i
        return x
    else:
        return "unpaired"


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                seq = ''
            else:
                header = i[1:]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    return Dict


if args.mode == "single":
    print("Concatenating reads: single mode")
    outfile = open(args.reads_directory + "/combined_reads.fastq", "w")
    reads = os.listdir(args.reads_directory)
    for i in reads:
        ext = lastItem(i.split("."))
        if ext == args.reads_ext and i != "combined_reads.fastq":
            readFile = open(args.reads_directory + "/" + i, "r")
            for line in readFile:
                outfile.write(line)
    outfile.close()


if args.mode == "paired":
    print("Concatenating reads: paired mode...")
    print("processing forward reads...")
    outfile1 = open(args.reads_directory + "/combined_reads_1.fastq", "w")
    reads = sorted(os.listdir(args.reads_directory))
    for i in reads:
        ext = lastItem(i.split("_"))
        extension = lastItem(i.split("."))
        if re.findall(r'1', ext) and extension == args.reads_ext and i != "combined_reads_1.fastq":
            print("---" + i)
            readFile = open(args.reads_directory + "/" + i, "r")
            for line in readFile:
                outfile1.write(line)
    outfile1.close()

    outfile2 = open(args.reads_directory + "/combined_reads_2.fastq", "w")
    reads = sorted(os.listdir(args.reads_directory))
    print("processing reverse reads...")
    for i in reads:
        ext = lastItem(i.split("_"))
        extension = lastItem(i.split("."))
        if re.findall(r'2', ext) and extension == args.reads_ext and i != "combined_reads_1.fastq" and i != "combined_reads_2.fastq":
            print("---" + i)
            readFile = open(args.reads_directory + "/" + i, "r")
            for line in readFile:
                outfile2.write(line)
    outfile2.close()

    outfile = open(args.reads_directory + "/combined_reads_unpaired.fastq", "w")
    reads = os.listdir(args.reads_directory)
    print("processing unpaired reads...")
    for i in reads:
        ext = lastItem(i.split("_"))
        extension = lastItem(i.split("."))
        if not re.findall(r'1', ext) and extension == args.reads_ext and not re.findall(r'2', ext) and i != "combined_reads_1.fastq" and \
                        i != "combined_reads_2.fastq" and i != "combined_reads_unpaired.fastq":
            print("---" + i)
            readFile = open(args.reads_directory + "/" + i, "r")
            for line in readFile:
                outfile.write(line)
    outfile.close()


# numReads = 0
# print("Done. Counting reads")
# reads = open(args.reads_directory + "/combined_reads.fastq", "r")
# for i in reads:
#     if re.match(r'^@', i):
#         numReads += 1


print("Concatenating bins...")
bins = os.listdir(args.bin_directory)
outfile = open(args.bin_directory + "/combined_contigs.fasta", "w")
CombinedBinsDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in bins:
    ext = lastItem(i.split("."))
    if ext == args.bin_ext and i != "combined_contigs.fasta":
        bin = open(args.bin_directory + "/" + i, "r")
        for line in bin:
            outfile.write(line)
outfile.close()


print("Building bowtie index from the combined bins' contigs")
os.system("bowtie2-build --threads " + str(args.num_threads) + " -f " + args.bin_directory + "/combined_contigs.fasta "
          + args.bin_directory + "/combined_contigs.fasta --quiet")

print("Mapping reads to combined contigs")
if args.mode == "single":
    os.system("bowtie2 -t --threads " + str(args.num_threads) + " -x " + args.bin_directory + "/combined_contigs.fasta"
              + "-U " + args.reads_directory + "/combined_reads.fastq" + " "
              "-S combined_contigs.sam --no-unal --quiet --reorder")

if args.mode == "paired":
    os.system("bowtie2 -t --threads " + str(args.num_threads) + " -x " + args.bin_directory +
              "/combined_contigs.fasta -1 " + args.reads_directory + "/combined_reads_1.fastq -2 " +
              args.reads_directory + "/combined_reads_2.fastq -U " + args.reads_directory +
              "/combined_reads_unpaired.fastq -S combined_contigs.sam --no-unal --quiet --reorder")


print("Running SAMtools...")
os.system("samtools view -bS -o combined_contigs.bam combined_contigs.sam")
print("sorting...")
os.system("samtools sort -@ " + args.num_threads + " -o combined_contigs.bam.sorted combined_contigs.bam")
print("indexing...")
os.system("samtools index -b combined_contigs.bam.sorted combined_contigs.bam.sorted.bai")
print("summarizing results...")
os.system("samtools idxstats combined_contigs.bam.sorted >> samtools.idxstats.out")

print("Sorting output...")
bins = os.listdir(args.bin_directory)
binDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in bins:
    ext = lastItem(i.split("."))
    if ext == args.bin_ext and i != "combined_contigs.fasta":
        bin = open(args.bin_directory + "/" + i, "r")
        for line in bin:
            if re.match(r'^>', line):
                binDict[line.rstrip()[1:]] = i

print("Calculating reads/bin from combined contigs map\n...\n")
CovDictCombined = defaultdict(list)
idxTotal = open("samtools.idxstats.out")
for i in idxTotal:
    ls = i.rstrip().split("\t")
    contig = ls[0]
    if contig != "*":
        bin = binDict[contig]
        CovDictCombined[bin].append(ls[2])

print("Mapping reads to each individial bin\n...\n")
CovDictInd = defaultdict(list)
bins = os.listdir(args.bin_directory)
for i in bins:
    ext = lastItem(i.split("."))
    if ext == args.bin_ext:
        print("Building bowtie index for " + i)
        os.system("bowtie2-build --threads " + str(args.num_threads) + " -f " + args.bin_directory + "/" + i + " "
                  + args.bin_directory + "/" + i + " --quiet")
        print("Mapping reads to " + i)

        if args.mode == "single":
            os.system("bowtie2 -t --threads " + str(args.num_threads) + " -x " + args.bin_directory + "/" + i + " "
              "-U " + args.reads_directory + "/combined_reads.fastq" + " "
              "-S " + i + ".sam --no-unal --quiet --reorder")

        if args.mode == "paired":
            os.system("bowtie2 -t --threads " + str(args.num_threads) + " -x " + args.bin_directory + "/" + i +
                      " -1 " + args.reads_directory + "/combined_reads_1.fastq -2 " + args.reads_directory +
                      "/combined_reads_2.fastq -U " + args.reads_directory + "/combined_reads_unpaired.fastq -S " + i +
                      ".sam --no-unal --quiet --reorder")

        os.system("samtools view -bS -o " + i + ".bam " + i + ".sam")
        os.system("samtools sort -@ " + args.num_threads + " -o " + i + ".bam.sorted " + i + ".bam")
        os.system("samtools index -b " + i + ".bam.sorted " + i + ".bam.sorted.bai")
        os.system("samtools idxstats " + i + ".bam.sorted >> " + i + ".samtools.idxstats.out")
        idx = open(i + ".samtools.idxstats.out")
        print("Calculating reads for " + i)
        for line in idx:
            ls = line.rstrip().split("\t")
            CovDictInd[i].append(ls[2])
        print("")

print("Finished read recruitment to individual bins")

if args.bin_comparisons == "NA":
    pass
else:
    CovDictComp = defaultdict(lambda: defaultdict(list))
    count = 0
    print("Beginning individual bin comparisons from provided file: " + args.bin_comparisons)
    comparisons = open(args.bin_comparisons, "r")
    for line in comparisons:
        count += 1
        string = "bin_comparison_" + str(count) + ".fasta"
        ls = line.rstrip().split("\t")
        BINDICT = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
        for bin in ls:
            os.system("less " + args.bin_directory + "/" + bin + " >> " + args.bin_directory + "/" + string)
            binFile = open(args.bin_directory + "/" + bin, "r")
            for i in binFile:
                if re.match(r'^>', i):
                    header = i.rstrip()[1:]
                    BINDICT[header] = bin
        os.system("bowtie2-build --threads " + str(args.num_threads) + " -f " + args.bin_directory + "/" + string + " "
                  + args.bin_directory + "/" + string + " --quiet")
        if args.mode == "single":
            os.system("bowtie2 -t --threads " + str(args.num_threads) + " -x " + args.bin_directory + "/" + string +
                      " -U " + args.reads_directory + "/combined_reads.fastq -S " + string +
                      ".sam --no-unal --quiet --reorder")

        if args.mode == "paired":
            os.system("bowtie2 -t --threads " + str(args.num_threads) + " -x " + args.bin_directory + "/" + string +
                      " -1 " + args.reads_directory + "/combined_reads_1.fastq -2 " + args.reads_directory +
                      "/combined_reads_2.fastq -U " + args.reads_directory + "/combined_reads_unpaired.fastq -S " +
                      string + ".sam --no-unal --quiet --reorder")

        os.system("samtools view -bS -o " + string + ".bam " + string + ".sam")
        os.system("samtools sort -@ " + args.num_threads + " -o " + string + ".bam.sorted " + string + ".bam")
        os.system("samtools index -b " + string + ".bam.sorted " + string + ".bam.sorted.bai")
        os.system("samtools idxstats " + string + ".bam.sorted >> " + string + ".samtools.idxstats.out")
        idx = open(string + ".samtools.idxstats.out")
        print("Calculating reads for " + string)
        for LINE in idx:
            LS = LINE.rstrip().split("\t")
            if LS[0] != "*":
                bin = BINDICT[LS[0]]
                CovDictComp[string][bin].append(LS[2])
        print("")


print("writing final file...")
out = open(args.output_csv, "w")
if args.bin_comparisons == "NA":
    out.write("bin" + "," + "read_proportion_when_mapped_together" + "," + "read_proportion_when_mapped_separately" + "\n")
else:
    out.write("bin" + "," + "read_proportion_when_mapped_together" + "," + "read_proportion_when_mapped_separately" + ",")
    for i in CovDictComp.keys():
        out.write(i.split(".")[0] + ",")
    out.write("\n")

# print("")
# print("total number of reads in the provided fastq file: " + str(numReads))
# print("")

for i in CovDictCombined.keys():
    bin = i

    CovListCombined = CovDictCombined[i]
    totalReadsCombined = sum(CovListCombined)
    CovListInd = CovDictInd[i]
    totalReadsInd = sum(CovListInd)
    print(bin + "recruited " + str(totalReadsCombined) + " reads during mapping to all concatenated bins")
    print(bin + "recruited " + str(totalReadsInd) + " reads during mapping only to this bin")
    if args.bin_comparisons == "NA":
        out.write(bin + "," + str(totalReadsCombined) + "," + str(totalReadsInd) + "\n")

    else:
        out.write(bin + "," + str(totalReadsCombined) + "," + str(totalReadsInd) + ",")
        for j in CovDictComp.keys():
            if len(CovDictComp[j][bin]) == 0:
                out.write("NA" + ",")
            else:
                CovListComp = CovDictComp[j][bin]
                totalReadsComp = sum(CovListComp)
                out.write(str(totalReadsComp) + ",")
                print(bin + "recruited " + str(totalReadsComp) + " reads during bin_comparison_" + str(j.split(".")[0]) + "\n")
    out.write("\n")
    print("")


os.system("mkdir ReadHog-auxilary_files")
os.system("mv " + args.reads_directory + "/combined_reads.fastq ReadHog-auxilary_files")
os.system("mv " + args.bin_directory + "/combined_contigs.fasta* ReadHog-auxilary_files")
if args.bin_comparisons != "NA":
    os.system("mv " + args.bin_directory + "/bin_comparison* ReadHog-auxilary_files")
for i in bins:
    ext = lastItem(i.split("."))
    if ext == args.bin_ext:
        os.system("mv " + args.bin_directory + "/" + i + ".* ReadHog-auxilary_files")
        os.system("mv " + i + ".bam ReadHog-auxilary_files")
        os.system("mv " + i + ".sam ReadHog-auxilary_files")
        os.system("mv " + i + ".samtools.idxstats.out ReadHog-auxilary_files")
        os.system("mv " + i + ".bam.sorted* ReadHog-auxilary_files")


print("Thank you for using ReadHog!")
