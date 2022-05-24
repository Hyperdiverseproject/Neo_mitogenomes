#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Thu Dec 16 2021

1.mapping reads against hit contigs file with bowtie2
2.sort sam file and transform it into bam file wirth n option
3.fixmate bam file
4.sort bam file by coordinates
5.index bam file
6.calculate depth of coverage

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os
import sys
import multiprocessing

def mapping(split, output_folder):

		variables=dict(
		read1='/media/thomas/DATA/transcriptomes_raw_data_mitogenomes/Transcriptomes_Raw_Data/' + split[0],
		read2='/media/thomas/DATA/transcriptomes_raw_data_mitogenomes/Transcriptomes_Raw_Data/' + split[1],
		contig_blast='/media/thomas/Data2/Boulot/These/Bioinfo/Mitogenomes/09-remapping_blast_hit_contigs/ALL_contigs/ALL_hit_contigs_merged/' + split[2],
		sample=split[3],
		build_bowtie2=split[3]+'_build',
		out_paired=split[3]+'_out_paired',
		outfile_sort=split[3]+'_sorted',
		output_directory=output_folder
		)

		commands="""
		bowtie2-build {contig_blast} {build_bowtie2}
		bowtie2 -x {build_bowtie2} -1 {read1} -2 {read2} --local --very-sensitive-local --no-discordant -p 15 -S {out_paired}.sam 2> {sample}_paired.stderr
		echo "bowtie2 done"
		samtools sort -@ 16 -n -o {outfile_sort}.bam {out_paired}.sam
		echo "file sorted with -n option"
		samtools fixmate -m {outfile_sort}.bam {outfile_sort}_fixmate.bam
		echo "fixmate done"
		samtools sort -@ 16 -o {outfile_sort}_fixmate_coordinate.bam {outfile_sort}_fixmate.bam
		samtools index {outfile_sort}_fixmate_coordinate.bam
		echo "file sorted vanilla"
		samtools coverage -o {sample}_blastn_mapping_coverage_table.txt {outfile_sort}_fixmate_coordinate.bam
		echo "coverage done"
		rm *.bt2 *.bai {out_paired}.sam {outfile_sort}.bam {outfile_sort}_fixmate_coordinate.bam {outfile_sort}_fixmate.bam
		mv {sample}* {output_directory}
		""".format(**variables)

		command_list= commands.split('\n')
		for line in command_list:
			os.system(line)

def main():

	ref_file = open('ref_reads_ALL_contigs.txt', 'r')

	for line in ref_file:
		split=line.strip().split('\t')
		output_folder=split[3]+'_ALL_contigs_mapping_samtools_coverage_no_markdup'
		mkdir(output_folder)
		mapping(split, output_folder)

if __name__ == '__main__':
    main()