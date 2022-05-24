#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Mon Jan 17 2022

1.read_fasta and write_fasta are two functions used to interleaved fasta files
2.removing mitogenome reference sequence from file to keep only transcriptome contigs
3.genes from mitogenome reference are added at the contig file
4.mafft alignment of contig + genes file against mitogenome reference sequence using addfragments option

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os
import re
import sys
import fileinput

def read_fasta(fasta):
    with open(fasta, 'r') as fast:
        headers, sequences = [], []
        for line in fast:
            if line.startswith('>'):
                head = line.replace('>','').strip()
                headers.append(head)
                sequences.append('')
            else :
                seq = line.strip()
                if len(seq) > 0:
                    sequences[-1] += seq
    return (headers, sequences)

def write_fasta(headers, sequences, fasta):
    with open(fasta, 'w') as fast:
        for i in range(len(headers)):
            fast.write('>' + headers[i] + '\n' + sequences[i] + '\n')

def removing_ref(contig_interleaved):

	input_file=open(contig_interleaved, 'r')
	contigs_only=open(contig_interleaved.replace('_interleaved.fasta', '_contigs_only.fasta'), 'w')

	count_line=0

	for line in input_file:
	    if count_line ==0:
	        name_ref = line.strip()
	        count_line += 1 
	    elif count_line ==1:
	        seq_ref = line.strip()
	        count_line+=1
	    elif count_line%2==0:
	        name_contig=line.strip()
	        count_line+=1
	    elif count_line%2!=0:
	        seq_contig=line.strip()     
	        count_line+=1
	        contigs_only.write(name_contig + '\n' + seq_contig + '\n')

	input_file.close()
	contigs_only.close()

def merging_files(contigs_only, genes, merged_file):

    input_genes=open(genes, 'r')
    input_contig=open(contigs_only, 'r')
    output_file=open(merged_file, 'a')

    for contig in input_contig:
    	if contig.strip().startswith('>'):
    		seq_contig=next(input_contig)
    		name_contig=contig
    		output_file.write(name_contig+seq_contig)
    	else:
    		pass
    for genes in input_genes:
    	if genes.strip().startswith('>'):
    		seq_genes=next(input_genes)
    		name_genes=genes
    		output_file.write(name_genes+seq_genes)
    	else:
    		pass

    input_genes.close()
    input_contig.close()
    output_file.close()

def mafft(split, merged_file):

	variables=dict(
	contigs=merged_file,
	ref='./mitogenomes_interleaved/' + split[0],
	final_alignment= merged_file.replace('.fasta', '_aligned.fasta')
	)

	commands="""
	mafft --thread 16 --localpair --addfragments {contigs} {ref} > {final_alignment}
	""".format(**variables)

	command_list= commands.split('\n')
	for line in command_list:
		os.system(line)

def main():

	ref_file=open('ref_ALL_contigs_mafft_alignments_mito_genes.txt', 'r')
	merged_dir='./mito_genes_contigs_merged'
	mkdir(merged_dir)
	mafft_dir='./mafft_alignments_ref_mito_genes_contigs'
	mkdir(mafft_dir)

	for line in ref_file:
		split=line.strip().split('\t')
		ref=split[0]
		contig='./mafft_alignments_sorted_cleaned_realigned/kept_samples/'+split[1]
		genes='./mitogenomes_genes/'+split[2]
		headers, sequences=read_fasta(contig)
		contig_interleaved=contig.replace('_first_alignment_sorted_contigs_final_alignment_cleaned_contigs_aligned.fasta', '_interleaved.fasta')
		write_fasta(headers, sequences, contig_interleaved)
		contigs_only=contig_interleaved.replace('_interleaved.fasta', '_contigs_only.fasta')
		removing_ref(contig_interleaved)
		merged_file=contigs_only.replace('.fasta', '_mito_genes_merged.fasta')
		merging_files(contigs_only, genes, merged_file)
		alignment_file=merged_file.replace('.fasta', '_aligned.fasta')
		mafft(split, merged_file)
		os.system('rm '+contig_interleaved)
		os.system('mv '+alignment_file+' '+mafft_dir)
		os.system('mv '+merged_file+' '+merged_dir)

if __name__ == '__main__':
    main()
