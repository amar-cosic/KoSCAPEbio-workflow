from Bio import SeqIO

def extract_positions(fasta_file, positions):
    sequences = SeqIO.parse(fasta_file, "fasta")
    for record in sequences:
        extracted_sequence = ""
        previous_pos = None

        for pos in positions:
            if previous_pos is not None and pos - previous_pos > 1:
                extracted_sequence += "xx"
            extracted_sequence += record.seq[pos-1]
            previous_pos = pos
        print(f">{record.id}")
        print(extracted_sequence)
positions = [89, 99, 100, 101, 102, 103, 104, 114, 116, 117, 118, 119, 120, 121, 133, 269, 297, 386]
fasta_file = "aligned_V3-V4_final.fasta"
extract_positions(fasta_file, positions)
