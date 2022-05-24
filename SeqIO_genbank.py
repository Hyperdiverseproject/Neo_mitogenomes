#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''

created on Wen Sep 22 2021

extracting protein sequences from genbank files

author: Thomas Lemarcis

'''

from Bio import SeqIO
import argparse
from os import listdir
from os import mkdir
import os

def get_args():
    parser = argparse.ArgumentParser(description='parse_genbank_files')
    parser.add_argument('data_dir', help = 'path to data directory')
    
    return parser.parse_args()

def get_cds_feature_with_qualifier_value(gb_file, name, value):
    """Function to look for CDS feature by annotation value in sequence record.
    
    e.g. You can use this for finding features by locus tag, gene ID, or protein ID.
    """
    genome_record = SeqIO.read(gb_file, "genbank")

    for feature in genome_record.features:
        if feature.type == "CDS" and value in feature.qualifiers.get(name, []):

            return feature

    return None

def parse_genbank_files(gb_file, proteins, protein_tags):

    genome_record = SeqIO.read(gb_file, "genbank")

    with open(proteins, "w") as aa_output:
        for tag in protein_tags:
            cds_feature = get_cds_feature_with_qualifier_value(gb_file, "protein_id", tag)
            gene_sequence = cds_feature.extract(genome_record.seq)
            name_gene=cds_feature.qualifiers.get("gene")[0]
            name=genome_record.name
            protein_sequence = cds_feature.qualifiers.get("translation")[0]
            assert protein_sequence == cds_feature.qualifiers["translation"][0]
            aa_output.write(">%s_%s_%s\n%s\n" % (name,name_gene, tag, protein_sequence))

    print("Done")
    aa_output.close()

def main():

    args = get_args()
    input_dir = args.data_dir
    output_dir = input_dir+'_sorted'
    mkdir(output_dir)
    list_file= listdir(input_dir)
    
    for file in list_file:
        protein_tags=[]
        gb_file = input_dir +'/'+ file
        gb_file_read=open(gb_file, 'r')
        proteins= output_dir+ '/' + file.split('.gb')[0]+'_proteins.fasta'
        for line in gb_file_read:
            if line.strip().startswith('/protein_id'):
                line_clean=line.replace('\n','')
                line_clean2= line_clean.replace('"','')
                split=line_clean2.split('=')
                protein_tags.append(split[1])
        parse_genbank_files(gb_file, proteins, protein_tags)
        gb_file_read.close()

if __name__ == '__main__':
    main()