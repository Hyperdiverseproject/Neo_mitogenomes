#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Mon Oct 18 2021

1.read_fasta and write_fasta are two functions used to interleaved fasta files
2.blastn contigs used as reference against mitogenome sequence
3.parsing outputfile from blastn and extracting hit contig names
4.searching and extracting sequences of hit contigs in original tanscritome file 

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os

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
    dict_transcript={}
   
    for i in range(len(headers)):
    	dict_transcript[headers[i].split(' ')[0]] = sequences[i]
    		
    return (headers, sequences, dict_transcript)

def write_fasta(headers, sequences, fasta):
    with open(fasta, 'w') as fast:
        for i in range(len(headers)):
            fast.write('>' + headers[i] + '\n' + sequences[i] + '\n')

def blastn(split):

		variables=dict(
		transcri=split[0],
		mito=split[1],
		accession=split[2],
		output=split[1][:-6] + "_blastn.fasta",
		path_db=split[0]
		)

		commands="""
		makeblastdb -in {transcri} -parse_seqids -dbtype nucl
		blastn -query {mito}  -outfmt 7 -out {output} -db {path_db} -num_threads 8 -evalue 1e-10
		""".format(**variables)

		command_list= commands.split('\n')

		for line in command_list:
			os.system(line)

def parse_blastn(output, headers, sequences):
	blastn_results=open(output, 'r')
	hit_transcriptome=[]
	for line in blastn_results:
		line =line.strip()
		if line[0] != '#':
			split2 =line.split('\t')
			if float(split2[2]) > 50:
				if split2[1] not in hit_transcriptome: #remove duplicate
					hit_transcriptome.append(split2[1])
	return(hit_transcriptome)

def search_transcript(hit_transcriptome, dict_transcript, hit_contig_name):

	sorted_hit_contig=open(hit_contig_name,'w')

	for hit in hit_transcriptome:
		if hit in dict_transcript.keys():
			print(hit)
			sorted_hit_contig.write('>'+ hit+ '\n' + dict_transcript[hit] + '\n')
		else:
			print('losers')

def main():

	ref_file = open('ref_all_transcriptomes_mitogenomes_blastn.txt', 'r') #we used an external text file with all file names

	for line in ref_file:
		split=line.strip().split('\t')
		transcri=split[0]
		mito=split[1]
		accession=split[2]
		output=split[1][:-6]+"_blastn.fasta"
		output_folder=split[0][:-6]+'_blastn_output_'+accession
		headers, sequences, dict_transcript = read_fasta(transcri)
		transcri_out=transcri.replace('.fasta', '_interleaved.fasta')
		write_fasta(headers, sequences, transcri_out)
		blastn(split)
		hit_transcriptome = parse_blastn(output, headers, sequences)
		hit_contig_name= split[2]+'_'+transcri.replace('.fasta', '_sorted_blastn_hit_contig.fasta')
		search_transcript(hit_transcriptome,dict_transcript, hit_contig_name)
		os.system('rm *.nhr *.nin *.nog *.nsd *.nsi *.nsq *_interleaved.fasta')
		mkdir(output_folder)
		os.system('cp '+ output + ' ' +  hit_contig_name + ' ' + output_folder)

if __name__ == '__main__':
    main()

