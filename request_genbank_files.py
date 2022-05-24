#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Wen Sep 22 2021

extracting genbank files of mitogenome reference using the accession numbers

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os

def esearch_efetch(accession_numbers, output_dir):

		variables=dict(
		accession=accession_numbers,
		output= accession_numbers+ '.gb',
		output_dir=output_dir
		)

		commands="""
		esearch -db nuccore -query {accession} | edirect -fetch -format gp > {output_dir}/{output}
		""".format(**variables)

		command_list= commands.split('\n')
		print(command_list)

		for line in command_list:
			os.system(line)

def main():

	output_dir = 'mitogenomes_genbank_files'
	mkdir(output_dir)
	ref_file = open('accession_numbers.txt', 'r')

	for line in ref_file:
		accession_numbers=line.replace('\n', '')
		esearch_efetch(accession_numbers, output_dir)

if __name__ == '__main__':
    main()