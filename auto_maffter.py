#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Mon Jan 10 2022

1.read_fasta and write_fasta are two functions used to interleaved fasta files
2.first mafft alignment of hit contigs sequences from blastn and blastp against mitogenome reference
3.calculating number of nucleotides differences between reference sequence and each contigs and keeping only contig sequences that have less than 50% of differences
4.mitogenome reference sequence is removed to keep only contigs
5.final mafft alignement with kept contigs against mitogenome reference

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

def mafft_1(split):

	variables=dict(
	contigs='./ALL_hit_contigs_merged/'+split[1],
	ref='./mitogenomes_interleaved/' + split[0],
	first_alignment= split[1].replace('.fasta', '_first_alignment.fasta')
	)

	commands="""
	mafft --thread 16 --localpair --addfragments {contigs} {ref} > {first_alignment}
	""".format(**variables)

	command_list= commands.split('\n')
	for line in command_list:
		os.system(line)

def PID_calculation(first_alignment_interleaved):

	inputFH=open(first_alignment_interleaved, 'r')
	outputfile_PID=open(first_alignment_interleaved.replace('.fasta', '_sorted.fasta'), 'a')
	summary_file=open(first_alignment_interleaved.replace('.fasta', '_summary.txt'), 'w')

	summary_file.write('ref'+ '\t'+'contig'+ '\t'+ 'genomic_distance'+'\n')

	count_line=0
	PID=50

	for line in inputFH:
	    if count_line ==0:
	        name_ref = line.strip()
	        count_line += 1 
	    elif count_line ==1:
	        seq_ref = line.strip()
	        outputfile_PID.write(name_ref + '\n'+ seq_ref+'\n')
	        count_line+=1
	    elif count_line%2==0:
	        name_contig=line.strip()
	        count_line+=1
	    elif count_line%2!=0:
	        seq_contig=line.strip()     
	        count_line+=1

	        counter_dist = 0
	        counter_distGAP=0
	        comp_dist=0
	        comp_distGAP=0

	        for i in range(len(seq_contig)):
	            if line[i] == '-' and seq_ref[i] != '-':
	                counter_distGAP +=1
	                comp_distGAP +=1
	            elif line[i] != '-' and seq_ref[i] =='-':
	                counter_distGAP += 1
	                comp_distGAP += 1
	            elif line[i] !='-' and seq_ref[i] == line[i]:
	                comp_dist+=1
	            elif line[i] != '-' and seq_ref[i] != line[i]:
	                counter_dist +=1
	                comp_dist +=1
	               
	        dist_nuc = counter_dist*100/comp_dist
	        dist_gap = (counter_distGAP+counter_dist)*100/(comp_dist+comp_distGAP)
	        summary_file.write(str(name_ref).replace('>', '') + '\t'+ str(name_contig).replace('>', '') + '\t'  + str(dist_nuc)+'\n')
	        if dist_nuc < PID:
	            outputfile_PID.write(name_contig + '\n'+seq_contig+ '\n')

	inputFH.close()
	outputfile_PID.close()
	summary_file.close()

def removing_ref(outputfile_PID):

	inputfile_PID=open(outputfile_PID, 'r')
	outputfile_PID_contigs_only=open(outputfile_PID.replace('.fasta', '_contigs.fasta'), 'w')

	count_line=0

	for line in inputfile_PID:
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
	        outputfile_PID_contigs_only.write(name_contig + '\n' + seq_contig + '\n')

	inputfile_PID.close()
	outputfile_PID_contigs_only.close()

def mafft_2(split, outputfile_PID_contigs_only):

	variables=dict(
	contigs=outputfile_PID_contigs_only,
	ref='./mitogenomes_interleaved/' + split[0],
	final_alignment= outputfile_PID_contigs_only.replace('.fasta', '_final_alignment.fasta')
	)

	commands="""
	mafft --thread 16 --localpair --addfragments {contigs} {ref} > {final_alignment}
	""".format(**variables)

	command_list= commands.split('\n')
	for line in command_list:
		os.system(line)

def main():

	ref_file = open('ref_ALL_contigs_mafft_alignments_last_files.txt', 'r')
	mafft_folder='./mafft_alignments_sorted_last_files'
	mkdir(mafft_folder)
	summary_PID=mafft_folder+'/summary_PID'
	mkdir(summary_PID)

	for line in ref_file:
		split=line.strip().split('\t')
		ref=split[0]
		contig=split[1]
		first_alignment=split[1].replace('.fasta', '_first_alignment.fasta')
		mafft_1(split)
		headers, sequences = read_fasta(first_alignment)
		first_alignment_interleaved=first_alignment.replace('_aligned.fasta', '_interleaved.fasta')
		write_fasta(headers, sequences, first_alignment_interleaved)
		PID_calculation(first_alignment_interleaved)
		summary_file=first_alignment_interleaved.replace('.fasta', '_summary.txt')
		outputfile_PID=first_alignment_interleaved.replace('.fasta', '_sorted.fasta')
		removing_ref(outputfile_PID)
		outputfile_PID_contigs_only=outputfile_PID.replace('.fasta', '_contigs.fasta')
		final_alignment=outputfile_PID_contigs_only.replace('.fasta', '_final_alignment.fasta')
		mafft_2(split, outputfile_PID_contigs_only)
		os.system('rm '+ first_alignment+ ' ' +first_alignment_interleaved+ ' ' +outputfile_PID_contigs_only+ ' '+outputfile_PID)		
		os.system('mv '+ final_alignment + ' ' + mafft_folder)
		os.system('mv ' + summary_file + ' '+ summary_PID)

if __name__ == '__main__':
    main()