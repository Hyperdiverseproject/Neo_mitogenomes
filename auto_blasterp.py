#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Wen Sep 22 2021

1.read_fasta and write_fasta are two functions used to interleaved fasta files
2.blastp contigs used as reference against mitogenome sequence
3.parsing outputfile from blastn and extracting hit contig names
4.searching and extracting sequences of hit contigs in original tanscritome file

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os
import re

def read_fasta(fasta, fasta2):
    with open(fasta, 'r') as fast:
        headers, sequences = [], []
        for line in fast:
            if line.startswith('>'):
                head = line.replace('>','').strip()
                head2= head.replace('|', ' ').strip()
                headers.append(head2)
                sequences.append('')
            else :
                seq = line.strip()
                if len(seq) > 0:
                    sequences[-1] += seq
    dict_transcript={}
   
    for i in range(len(headers)):
    	dict_transcript[headers[i].split(' ')[1]] = sequences[i]
    		
    with open(fasta2, 'r') as fast2:
        headers2, sequences2 = [], []
        for line in fast2:
            if line.startswith('>'):
                head3 = line.replace('>','').strip()
                headers2.append(head3)
                sequences2.append('')
            else :
                seq = line.strip()
                if len(seq) > 0:
                    sequences2[-1] += seq
    dict_transcript2={}
   
    for i in range(len(headers2)):
    	dict_transcript2[headers2[i].split(' ')[0]] = sequences2[i]

    return (headers, sequences, dict_transcript, headers2, sequences2, dict_transcript2)

def write_fasta(headers, sequences, fasta, headers2, sequences2, fasta2):
    with open(fasta, 'w') as fast:
        for i in range(len(headers)):
            fast.write('>' + headers[i] + '\n' + sequences[i] + '\n')
    with open(fasta2, 'w') as fast2:
        for i in range(len(headers2)):
            fast2.write('>' + headers2[i] + '\n' + sequences2[i] + '\n')

def blastp(split):

		variables=dict(
		transcri_ORFfinder=split[0],
		mito=split[1],
		accession=split[2],
		output=split[0][:-6]+ "_" +split[1][:-6] + "_blastp.fasta",
		path_db=split[0]
		)

		commands="""
		makeblastdb -in {transcri_ORFfinder} -parse_seqids -dbtype prot
		blastp -query {mito}  -outfmt 7 -out {output} -db {path_db} -num_threads 8 
		""".format(**variables)

		command_list= commands.split('\n')

		for line in command_list:
			os.system(line)

def parse_blastp(output, headers, sequences):
	blastp_results=open(output, 'r')
	hit_transcriptome=[]
	for line in blastp_results:
		line =line.strip()
		if line[0] != '#':
			split2 =line.split('\t')
			split3=re.split('ORF\d*_', split2[1])
			split4=split3[1].split(':')
			if float(split2[2]) > 30:
				if split4[0] not in hit_transcriptome: #remove duplicate
					hit_transcriptome.append(split4[0])

	return(hit_transcriptome)

def search_transcript(hit_transcriptome, dict_transcript2, hit_contig_name):

	sorted_hit_contig=open(hit_contig_name,'w')

	for hit in hit_transcriptome:
		if hit in dict_transcript2.keys():
			print(hit)
			sorted_hit_contig.write('>'+ hit+ '\n' + dict_transcript2[hit] + '\n')
		else:
			print('losers')

def main():

	ref_file = open('ref_all_transcriptomes_mitogenomes_blastp.txt', 'r')

	for line in ref_file:
		split=line.strip().split('\t')
		transcri_ORFfinder=split[0]
		transcri=split[3]
		mito=split[1]
		accession=split[2]
		output=split[0][:-6]+ "_" +split[1][:-6] + "_blastp.fasta"
		output_folder=split[0][:-6]+'_blastp_output_'+accession
		headers, sequences, dict_transcript, headers2, sequences2, dict_transcript2 = read_fasta(transcri_ORFfinder, transcri)
		transcri_out=transcri_ORFfinder.replace('.fasta', '_interleaved.fasta')
		transcri_out2=transcri.replace('.fasta', '_interleaved.fasta')
		write_fasta(headers, sequences, transcri_out, headers2, sequences2, transcri_out2)
		blastp(split)
		hit_transcriptome = parse_blastp(output, headers, sequences)
		print('blaspt_parsed')
		hit_contig_name= split[1][:-6]+'_'+transcri_ORFfinder.replace('.fasta', '_sorted_blastp_hit_contig.fasta')
		search_transcript(hit_transcriptome,dict_transcript2, hit_contig_name)
		print('hit_contig_found')
		os.system('rm *.phr *.pin *.pog *.psd *.psi *.psq *_interleaved.fasta')
		mkdir(output_folder)
		os.system('cp '+ output + ' ' + hit_contig_name + ' ' + output_folder)

if __name__ == '__main__':
    main()