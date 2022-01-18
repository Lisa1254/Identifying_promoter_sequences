#Assignment 3: BINF6410, Lisa Hoeg

import sys
import os
import re
from random import randint

print("Welcome to promoter search program.\nRequired arguments are directory, promoter file, gene file.\nOutputs will include files of gene counts for each promoter, and notable promoters for each gene set.\n")

in_dir = sys.argv[1]
in_prom = sys.argv[2]
in_genes = sys.argv[3]

#Define functions: 

#Function for finding which genes on which chromosome
def chr_genes(filename):
    genes_list = [re.search("chromosome.([0-9]+|[MP]t)", filename).group()]
    file_open = open(filename)
    file_lines = file_open.readlines()
    file_open.close()
    for line in file_lines:
        if (line[0] != "#"):
            data = line.split("\t")
            if (data[2] == "gene"):
                if (re.search("Zm00001d[0-9]+", data[8]).group()):
                    gene_id = re.search("Zm00001d[0-9]+", data[8]).group()
                    gene_data = gene_id + ":" + data[3] + ":" + data[4] + ":" + data[6]
                    genes_list.append(gene_data)
    return(genes_list)

#Function for finding location of gene for genes of interest
def get_gene_loc(genes_search, chr_dict):
    gene_list = []
    find_data = {}
    for chr in chr_dict:
        for each in chr_dict[chr]:
                data = each.split(":")
                for gene in genes_search:
                    if gene == data[0]:
                        if data[3] == "+":
                            temp = [gene, data[1], data[3]]
                        else:
                            temp = [gene, data[2], data[3]]
                        gene_list.append(temp)
        find_data[chr] = gene_list
    return(find_data)

#Function to import sequences
def chr2seq(dir_path):
    seq_list = []
    chr_list = []
    for name in os.listdir(dir_path):
        if (re.search(".fa(sta|)", name)):
            temp_ch_name = re.search("chromosome.([0-9]+|[MP]t)", name).group()
            chr_list.append(temp_ch_name)
            seq_file = open(in_dir+name)
            seq_lines = seq_file.readlines()
            seq_file.close()
            temp_seq = ""
            for q in range(1,len(seq_lines)):
                temp_data = seq_lines[q].rstrip("\n")
                temp_seq += temp_data
            seq_list.append(temp_seq)
    seq_dict = {}
    for w in range(0,len(chr_list)):
        seq_dict[chr_list[w]] = seq_list[w]
    return(seq_dict)

#Nested functions to move search window for promoter:
def move_cs(start, end):
    temp_start+=50
    temp_end+=50
    temp_seq = dict_chr_seq[chr][temp_start:temp_end]
    temp_seq = temp_seq.lower()
    temp_ns = temp_seq.count("N") + temp_seq.count("n")

def move_ncs(start, end):
    temp_start+=-50
    temp_end+=-50
    temp_seq = dict_chr_seq[chr][temp_start:temp_end]
    temp_seq = temp_seq[::-1]
    temp_seq = temp_seq.replace("A", "t")
    temp_seq = temp_seq.replace("C", "g")
    temp_seq = temp_seq.replace("G", "c")
    temp_seq = temp_seq.replace("T", "a")
    temp_ns = temp_seq.count("N") + temp_seq.count("n")

#Function to identify promoter sequences for genes of interest
def get_prom_seq(chr_data_loc, dict_chr_seq):
    gene_seqs_dict = {}
    out_seq_n = open(out_n_file,"w")
    for chr in chr_data_loc:
        for each in chr_data_loc[chr]:
            if each[2] == "+":
                temp_start = int(each[1])
                temp_end = int(each[1])+499
                temp_seq = dict_chr_seq[chr][temp_start:temp_end]
                temp_ns = temp_seq.count("N") + temp_seq.count("n")
                while (temp_ns > 100):
                    move_cs(temp_start, temp_end)
                temp_match = 0
                while (temp_match == 0):
                    for prom in prom_list:
                        if (temp_seq.find(prom)>-1):
                            temp_match+=1
                    move_cs(temp_start, temp_end)          
            else:
                temp_start = int(each[1])-499
                temp_end = int(each[1])
                temp_seq = dict_chr_seq[chr][temp_start:temp_end]
                temp_ns = temp_seq.count("N") + temp_seq.count("n")
                while (temp_ns > 100):
                    move_ncs(temp_start, temp_end)
                temp_match = 0
                while (temp_match == 0):
                    for prom in prom_list:
                        if (temp_seq.find(prom)>-1):
                            temp_match+=1
                    move_ncs(temp_start, temp_end)          
            gene_seqs_dict[each[0]] = temp_seq
    out_seq_n.close()
    return(gene_seqs_dict)

#Function to search for each promoter within sequence, add to count; 
# searches for longest promoters first before smaller ones, 
# as small ones are often contained within the longer sequences
def count_promoters(count_genes, counts_dict):
    prom_order = sorted(list(counts_dict.keys()), key=len)[::-1]
    for gene in count_genes:
        for prom in counts_dict:
            if (re.search(prom, count_genes[gene])):
                counts_dict[prom] += 1
                break

#Function to write counts to file
def write_order_promoters(pcounts, outfile):
    out_open_all = open(outfile+"_all.txt", "w")
    out_open_all.write("Counts by promoter for list of genes\n\n")
    out_open_note = open(outfile+"_note.txt", "w")
    out_open_note.write("Promoters with notably high or low counts\n\n")
    note_zero = []
    note_low = []
    note_high = []
    for promoter in pcounts:
        out_open_all.write(promoter + "\t" + str(pcounts[promoter]) + "\n")
        if (pcounts[promoter] == 0):
            note_zero.append(promoter)
        elif (pcounts[promoter] < 5):
            note_low.append(promoter)
        elif (pcounts[promoter] > 50):
            note_high.append(promoter)
    out_open_note.write("Promoters with high counts:\n")
    for promoter in note_high:
        out_open_note.write(promoter + "\t" + str(pcounts[promoter]) + "\n")
    out_open_note.write("Promoters with low counts:\n")
    for promoter in note_low:
        out_open_note.write(promoter + "\t" + str(pcounts[promoter]) + "\n")
    out_open_note.write("Promoters with no counts:\n")
    for promoter in note_low:
        out_open_note.write(promoter + "\n")
    out_open_all.close()
    out_open_note.close()

##Run program

#Get chromosome list
gff_list = []
for name in os.listdir(in_dir):
    if (name.find("gff")>-1):
        gff_list.append(name)

#Get genes on a chromosomes (uses function defined above)
dict_chr_genes = {}
for gff in gff_list:
    temp_list = chr_genes(in_dir+gff)
    dict_chr_genes[temp_list[0]] = temp_list[1:]

#Get supplied genes
genes_file = open(in_dir+in_genes)
genes_lines = genes_file.readlines()
genes_file.close()
genes_list_given = []
for each in genes_lines:
    temp_gene = each.rstrip('\n')
    genes_list_given.append(temp_gene)

#Make list of random genes
genes_list_random = []
while (len(genes_list_random) < len(genes_list_given)):
    rand_n1 = randint(0,len(dict_chr_genes)-1)
    rand_ch = list(dict_chr_genes.keys())[rand_n1]
    rand_n2 = randint(0,len(dict_chr_genes[rand_ch])-1)
    rand_g = dict_chr_genes[rand_ch][rand_n2]
    rand_g = rand_g.split(":")
    rand_g = rand_g[0]
    if ((rand_g not in genes_list_random) & (rand_g not in genes_list_given)):
        genes_list_random.append(rand_g)

#Find loc for genes of interest (uses function defined above)
chr_data = get_gene_loc(genes_list_given, dict_chr_genes)
chr_data_r = get_gene_loc(genes_list_random, dict_chr_genes)

#Import sequences (uses function defined above)
chr_seq_dict = chr2seq(in_dir)

#Set up dictionary to add locations and counts for each promoter
prom_file = open(in_dir+in_prom)
prom_lines = prom_file.readlines()
prom_file.close()
prom_list = []
prom_c_dict = {}
prom_c_dict_r = {}
for each in prom_lines:
    temp_prom = each.rstrip('\n')
    temp_prom = temp_prom.lower()
    prom_list.append(temp_prom)
    prom_c_dict[temp_prom] = 0
    prom_c_dict_r[temp_prom] = 0

#Identify sequence for promoter (uses function defined above)
gene_seqs = get_prom_seq(chr_data, chr_seq_dict, in_dir+"found_ns_given.txt")
gene_seqs_r = get_prom_seq(chr_data_r, chr_seq_dict, in_dir+"found_ns_random.txt")

#Count promoters (uses function defined above)
count_promoters(gene_seqs, prom_c_dict)
count_promoters(gene_seqs_r, prom_c_dict_r)

# Write files of count numbers (uses function defined above)
write_order_promoters(prom_c_dict, in_dir+"given_counts")
write_order_promoters(prom_c_dict_r, in_dir+"random_counts")

# The results seem to show that the genes list given and the random list of genes 
# have the same two highly expressed promoters (aaacca and aaag)
# It is possible that those are the most common promoters for housekeeping genes.
# Both gene sets had a lot of promoters that had no expression within the selected genes.
