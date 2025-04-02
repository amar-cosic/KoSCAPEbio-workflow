def split_sequences_by_nucleotide(input_file, output_file_without_n, output_file_with_n):
    with open(input_file, 'r') as infile, \
         open(output_file_without_n, 'w') as outfile_without_n, \
         open(output_file_with_n, 'w') as outfile_with_n:
        
        write_seq = None
        current_header = ''
        current_sequence = []
        
        for line in infile:
            if line.startswith('>'):
                if current_header:
                    sequence = ''.join(current_sequence)
                    if 'N' in sequence:
                        outfile_with_n.write(current_header)
                        outfile_with_n.write(sequence + '\n')
                    else:
                        outfile_without_n.write(current_header)
                        outfile_without_n.write(sequence + '\n')
                
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line.strip())
        if current_header:
            sequence = ''.join(current_sequence)
            if 'N' in sequence:
                outfile_with_n.write(current_header)
                outfile_with_n.write(sequence + '\n')
            else:
                outfile_without_n.write(current_header)
                outfile_without_n.write(sequence + '\n')

input_file = 'V3_V4_duplicated.fasta'
output_file_without_n = 'V3_V4_duplicated_without_N.fasta'
output_file_with_n = 'V3_V4_duplicated_with_N.fasta'

split_sequences_by_nucleotide(input_file, output_file_without_n, output_file_with_n)

print(f"Sequences have been split into {output_file_without_n} and {output_file_with_n} based on the presence of 'N' in the sequences.")
