#!/usr/bin/env python3

import os, gzip, itertools, re
import collections
import pandas

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

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

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

######### Using Bash's curl to get the fasta files if they are not in directory

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

######### Salmonella first ! ###########

# Remember file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"

S_num_genes_1=0

with gzip.open(file1,"rt") as fh: # opening fasta file, reading line by line, and counting # of genes based on ">"
    for line in fh:
        if ">" in line:
            S_num_genes_1 += 1

S_totalgene_length=0
S_GC_counts=0
S_codon ={}
S_condon_num = 0
bases = ['A','T','G','C']
S_codons = collections.defaultdict(int)


with gzip.open(file1,"rt") as fh:
    paired_seqs = aspairs(fh)
    seqs = dict(paired_seqs)
    for id,seq in seqs.items():
        S_totalgene_length += len(seqs[id])
        for base in seq:
            if base == "G" or base == "C":
                S_GC_counts += 1
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                cod = (base1 + base2 + base3)
                S_codon[cod] = 0
                set(S_codon)
    for id in seqs:
        for i in range(0, len(seqs[id]), 3): ## remember seqs[id] = seq (the value for the key id)
            S_codon[seqs[id][i:i+3]] += 1
            #print(codon)
    for c in S_codon:
        S_codons[c] += S_codon[c]


S_G_C_freq = (float(S_GC_counts/S_totalgene_length)) * 100

S_codon_freqs={}
S_totalcodons = 0

for c in S_codon:
    S_totalcodons += S_codon[c]

for c in S_codons:
    if c in S_codon_freqs:
        S_codon_freqs[c] += float(S_codons[c]/S_totalcodons) * 100
    else:
        S_codon_freqs[c] = float(S_codons[c]/S_totalcodons) * 100

#print(S_codon_freqs)

######### Mycobacteria next ! ###########
# Remember file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

M_num_genes_1=0

with gzip.open(file2,"rt") as fh:  # opening fasta file, reading line by line, and counting # of genes based on ">"
    for line in fh:
        if ">" in line:
            M_num_genes_1 += 1

M_totalgene_length=0

M_GC_counts=0

M_codon ={}
M_condon_num = 0
bases = ['A','T','G','C']
M_codons = collections.defaultdict(int)


with gzip.open(file2,"rt") as fh:
    paired_seqs = aspairs(fh)
    seqs = dict(paired_seqs)
    for id,seq in seqs.items():
        M_totalgene_length += len(seqs[id])
        for base in seq:
            if base == "G" or base == "C":
                M_GC_counts += 1
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                cod = (base1 + base2 + base3)
                M_codon[cod] = 0
                set(M_codon)
    for id in seqs:
        for i in range(0, len(seqs[id]), 3): ## remember seqs[id] = seq (the value for the key id)
            M_codon[seqs[id][i:i+3]] += 1
            #print(codon)
    for c in M_codon:
        M_codons[c] += M_codon[c]

M_G_C_freq = (float(M_GC_counts/M_totalgene_length)) * 100

M_codon_freqs={}
M_totalcodons = 0

for c in M_codon:
    M_totalcodons += M_codon[c]

for c in M_codons:
    if c in M_codon_freqs:
        M_codon_freqs[c] += float(M_codons[c]/M_totalcodons) * 100
    else:
        M_codon_freqs[c] = float(M_codons[c]/M_totalcodons) * 100

#print(M_codon_freqs)

########## Combining the Salmonella and Mycobacteria codon frequencies ##############
## combining each dictionary filled with codon frequencies
result = {}
for key in (S_codon_freqs.keys() | M_codon_freqs.keys()):
    if key in S_codon_freqs:
        result.setdefault(key, []).append(S_codon_freqs[key])
    if key in M_codon_freqs:
        result.setdefault(key, []).append(M_codon_freqs[key])

######## Printing out all of the results ##########

print("The number of genes in Salmonella enterica subsp is {}".format(S_num_genes_1))
print("The total length of all the genes in Salmonella enterica is {}".format(S_totalgene_length))
print("The G+C% in Salmonella enterica CDS sequences is {:.4}%".format(S_G_C_freq))
print("The total number of codons in Salmonella enterica CDS sequences is {}".format(S_totalcodons))

print()

print("The number of genes in Mycobacterium tuberculosis is {}".format(M_num_genes_1))
print("The total length of all the genes in Mycobacterium tuberculosis  is {}".format(M_totalgene_length))
print("The G+C% in Mycobacterium tuberculosis  CDS sequences is {:.4}%".format(M_G_C_freq))
print("The total number of codons in Mycobacterium tuberculosis  CDS sequences is {}".format(M_totalcodons))

print()

all_GC = ((M_GC_counts + S_GC_counts)/(S_totalgene_length + M_totalgene_length)) * 100
print("The total G+C% in both species' CDS sequences is {:.4}%".format(all_GC))

print()

print("Codon Frequency in Sp1 Frequency in Sp2")
for c, v in result.items():
    print("{}   {:.3}%  {:.3}%".format(c, v[0], v[1]))
