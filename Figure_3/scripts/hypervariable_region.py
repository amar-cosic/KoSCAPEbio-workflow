import os
from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import gzip
import argparse 
import pathlib
import sys#delete later
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import os
import subprocess as sb
import glob
import shutil
import re

import logging


from skbio import DNA
from skbio.alignment import local_pairwise_align_ssw

output_path="./"
records = list(SeqIO.parse(f"{output_path}/expanded_mock_comunity.fasta", "fasta"))
def find_primer_positions(primer, sequence, alignment_file):
    # Create DNA objects for your sequence and primer
    sequence_dna = DNA(sequence)
    primer_dna = DNA(primer)

    # Perform the alignment
    alignment, score, start_end_positions = local_pairwise_align_ssw(sequence_dna, primer_dna)

    # Check if any alignment is found
    if alignment is None:
        print(f"No alignments found for primer {primer} in sequence {sequence[:10]}...")  # Print first 10 characters of the sequence for identification
        return None, None

    # Write the alignment to the file
    with open(alignment_file, "a") as file:
        file.write(str(alignment))

    # Get the start and end positions
    start = start_end_positions[0][0]
    end = start_end_positions[0][1]

    print(f"Alignment found for primer {primer} in sequence {sequence[:10]}... Start: {start}, End: {end}")

    return start, end

primer_341F = "CCTACGGGNGGCWGCAG"  
primer_805R = "ATTAGAWACCCBNGTAGTCC"  

primer_515F = "GTGCCAGCMGCCGCGGTAA"  
primer_926R = "AAACTYAAAKGAATTGACGG" 

primer_806R = "ATTAGAWACCCBDGTAGTCC" 
                

alignment_file =f"{output_path}/alignments.txt" 
v3_v4_records = []
v4_v5_records = []    
v4_records = []
v3_v5_records=[]
for record in records:
    seq_str = str(record.seq)
    print(f"Working on sequence: {seq_str[:10]}...") 
    
    # find positions of primers
    print("Finding primer 341F...")
    start_341F, end_341F = find_primer_positions(primer_341F, seq_str, alignment_file)
    
    print("Finding primer 805R...")
    start_805R, end_805R = find_primer_positions(primer_805R, seq_str, alignment_file)
    
    print("Finding primer 515F...")
    start_515F, end_515F = find_primer_positions(primer_515F, seq_str, alignment_file)
    
    print("Finding primer 926R...")
    start_926R, end_926R = find_primer_positions(primer_926R, seq_str, alignment_file)

    print("Finding primer 806R...")
    start_806R, end_806R = find_primer_positions(primer_806R, seq_str, alignment_file)


    # extract the regions between primers
    print("Extracting regions between primers...")
    v3_v4_seq = record[end_341F:end_805R]
    v4_v5_seq = record[end_515F:end_926R]
    v4_seq = record[end_515F:start_806R]
    v3_v5_seq = record[end_341F:end_926R]

    # append to records lists if length is less than or equal to 600 bp
    if len(v3_v4_seq) <= 500:
        v3_v4_records.append(v3_v4_seq)
    else:
        print(f"V3-V4 sequence for record {record.id} was longer than 600bp and has been skipped")

    if len(v4_v5_seq) <= 500:
        v4_v5_records.append(v4_v5_seq)
    else:
        print(f"V4-V5 sequence for record {record.id} was longer than 600bp and has been skipped")
    
    if len(v4_seq) <= 350:
        v4_records.append(v4_seq)
    else:
        print(f"V4 sequence for record {record.id} was longer than 600bp and has been skipped")

    if len(v3_v5_seq) <= 700:
        v3_v5_records.append(v3_v5_seq)
    else:
        print(f"V3-V5 sequence for record {record.id} was longer than 600bp and has been skipped")


# Write the extracted regions to separate files
SeqIO.write(v3_v4_records, f"{output_path}/V3_V4_duplicated.fasta", "fasta")
SeqIO.write(v4_v5_records, f"{output_path}/V4_V5_duplicated.fasta", "fasta")
SeqIO.write(v4_records, f"{output_path}/V4_duplicated.fasta", "fasta")
SeqIO.write(v3_v5_records, f"{output_path}/V3_V5_duplicated.fasta", "fasta")


from Bio import SeqIO

def remove_duplicates(input_file, output_file):
    seen_ids = set()
    unique_records = []

    for record in SeqIO.parse(input_file, "fasta"):
        if record.id not in seen_ids:
            seen_ids.add(record.id)
            unique_records.append(record)

    with open(output_file, "w") as output_handle:
        SeqIO.write(unique_records, output_handle, "fasta")

# List of file names
file_names = ["V3_V4_duplicated.fasta", "V4_V5_duplicated.fasta", "V4_duplicated.fasta", "V3_V5_duplicated.fasta"]

for file_name in file_names:
    input_file = os.path.join(output_path, file_name)
    output_file = os.path.join(output_path, file_name.replace("_duplicated.fasta", ".fasta"))
    
    remove_duplicates(input_file, output_file)
