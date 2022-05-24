#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Tue Nov 23 2021

this script was write to merged hit contigs from blastn and blastp output files

author: Thomas Lemarcis

'''

import argparse
from os import listdir
from os import mkdir
import os
import re
import sys
import fileinput

def get_args():

    parser = argparse.ArgumentParser(description='merging hit blastn and blastp contig files')
    parser.add_argument('blastn_dir', help='path to contigs directory')
    parser.add_argument('blastp_dir', help='path to mitogenomes directory')
    
    return parser.parse_args()

def merging_files(list_merged_dir, merged_dir, list_blastp, blastp_dir):

    for blastn in list_merged_dir:
        blastn_file=blastn
        blastn_ID=blastn.split('_')[0:3]
   
        for blastp in list_blastp:
            blastp_file=blastp
            blastp_split=blastp.split('_')
            blastp_ID=[blastp_split[i] for i in (0,2,3)]
            
            input_FH=open(blastp_dir+'/'+blastp_file, 'r')
            output_FH=open(merged_dir+'/'+blastn_file, 'a')
            output_FH_read=open(merged_dir+'/'+blastn_file, 'r')

            list_name=[]
            for line in output_FH_read:
                if line.startswith('>'):
                    list_name.append(line)

            for line in input_FH:
                if line.strip().startswith('>'):
                    seq=next(input_FH)
                    name_seq=line

                    if blastn_ID == blastp_ID:
                        if name_seq not in list_name:
                            output_FH.write(name_seq+seq)
                        else:
                            pass

def main():

    args=get_args()
    blastn_dir=args.blastn_dir
    blastp_dir=args.blastp_dir
    list_blastn=listdir(blastn_dir)
    list_blastp=listdir(blastp_dir)
    merged_dir=blastn_dir+'_'+blastp_dir+'_merged'
    mkdir(merged_dir)

    for file in list_blastn:
        os.system('cp '+blastn_dir+'/'+file+' '+merged_dir+'/'+file.replace('_clustered_sorted_blastn_hit_contig.fasta', '_blastn_blastp_hit_contigs_merged.fasta'))

    list_merged_dir=listdir(merged_dir)
    merging_files(list_merged_dir, merged_dir, list_blastp, blastp_dir)

if __name__ == '__main__':
    main()