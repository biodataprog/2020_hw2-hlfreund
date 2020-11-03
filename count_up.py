#!/usr/bin/env python3

# Part 2 of HW 2

# Download this file
# Count up and print out the number genes (gene feature)
# Compute the total length of the genes (length is the END - START)
# Use the FASTA file to compute the total length of genome (by adding up the length of each sequence in the file). Recall I lectured on a basic code to read in a FASTA file - you can also see that code template here
# Print out the percentage of the genome which is coding

### Create objects with your filenames
gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"

fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa"
fasta_gz="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

# ---------------------------------------------------------------------------

import os,gzip,itertools,csv,re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>' # define that the line to be returned is the line with the ">"

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence


# utilizing UNIX OS to use curl and getfiles here
if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa")

# ---------------------------------------------------------------------------

### First the GFF File ... ###

### Below we count the number of genes (gene feature) - aka each row
### Compute total lenght of genes (lengh = End - start)
### Also going to count CDS (coding) regions and their length for later
gen_count = 0
tot_gen_len = 0
cds_count = 0
tot_cds_len = 0

with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue # skip rows that start with #
        if row[2] == ("gene"): # specifically counting genes
            gen_count += 1
            tot_gen_len += (int(row[4])-int(row[3]))
        if row[2] == ("CDS"): # specifically counting coding regions
            cds_count += 1
            tot_cds_len += (int(row[4])-int(row[3]))
    
    print("The total number of genes is {}.".format(gen_count))
    print("The total length of genes is {} nucleotides.".format(tot_gen_len))
    print("The total number of CDS regions is {}.".format(cds_count))
    print("The total gene length of coding regions is {} nucleotides.".format(tot_cds_len))
    

# ---------------------------------------------------------------------------

### Then the FASTA file ... ###

### Use FASTA file to compute total length of genome
## add up length of each sequence in file

import itertools
import gzip
import sys
import re
      
# based on post here
# https://drj11.wordpress.com/2010/02/22/python-getting-fasta-with-itertools-groupby/
            
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'
    
# this function reads in fasta file and returns pairs of data
# where the first item is the ID and the second is the sequence
# it isn't that efficient as it reads it all into memory
# but this is good enough for our project
def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

# ---------------

### If getting FASTA file from command line, do the next 7 lines
#if len(sys.argv) <= 1:
#    print("expecting on cmdline argument: fasta_parser.py FILE.fasta")
#    exit()

# ---------------

# here is my program
# get the filename from the cmdline
#filename = sys.argv[1]

# ---------------

### Unzip fasta file
#with open(fasta_gz, 'rb') as fd:   This is to unzip a gz file
#gzip_fd = gzip.GzipFile(fileobj=fd)
#destinations = pd.read(gzip_fd)

# ---------------

### Already unzipped fasta file, so continuing with fa file
filename=os.path.basename(fasta)

with open(filename,"r") as f:
    pairs = aspairs(f)
    seqs  = dict(pairs)
#    seqs = dict(aspairs(f))
            
### Iterate through the sequence(s) in FASTA file
#n=0
#for seqid in seqs:
#    print("id is ",seqid, "seq is ",seqs[seqid])
    
genome_length=0

for k,v in seqs.items():
    #print( "id is ",k,"seq is",v)
    #n += 1
    genome_length+=len(v)
    #print(seq_length)
print("The total length of the genome is {} nucleotides.".format(genome_length)) # in this case, there is only one genome in this file - so measuring length of whole genome

# ---------------------------------------------------------------------------

### Calculate percentage of genome that is a coding region
gen_cds_perc = (float(tot_cds_len)/float(genome_length)) * 100

print("The percentage of the genome that is coding is {:.4}%".format(gen_cds_perc))
