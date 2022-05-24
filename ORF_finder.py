#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Wen Sep 22 2021

performing ORFfinder translation on transcriptome files

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os

def get_args():
    parser = argparse.ArgumentParser(description='launch ORFfinder on transcriptomes')
    parser.add_argument('data_dir', help = 'path to data directory')
    
    return parser.parse_args()

def ORF_finder(transcri, transcri2, output_dir):

		variables=dict(
		transcri=transcri,
		output=transcri2.split('.fasta')[0]+'_ORF_finder.fasta',
		output_dir=output_dir
		)

		commands="""
		transcriptomes/ORFfinder -in {transcri} -g 1 -s 0 -ml 90 -out {output_dir}/{output} -outfmt 0 2>ORF.err 
		""".format(**variables)

		command_list= commands.split('\n')
		for line in command_list:
			os.system(line)

def main():
	args = get_args()
	input_dir = args.data_dir
	output_dir = input_dir+'_ORF_finder'
	mkdir(output_dir)
	list_file= listdir(input_dir)
	
	for file in list_file:
		if file != "ORFfinder":
			print("processing "+ file)
			transcri = input_dir+ '/' + file
			transcri2 = file     
			ORF_finder(transcri, transcri2, output_dir)
			print("Done")

if __name__ == '__main__':
    main()