from Bio import SeqIO
import math
import numpy as np
def calculate_shannon_entropy(sequences):
    num_sequences = len(sequences)
    length_of_sequences = len(sequences[0])
    entropy_list = []
    position_info_list = []
    for i in range(length_of_sequences):
        column = [sequences[j][i] for j in range(num_sequences)]
        counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'Others': 0, 'Gaps': 0}
        for char in column:
            if char == '-':
                counts['Gaps'] += 1
                continue
            if char in counts:
                counts[char] += 1
            else:
                counts['Others'] += 1
        column_no_gaps = [char for char in column if char != '-']
        if len(column_no_gaps) == 0:
            entropy_list.append(0)
            position_info_list.append([counts['A'], counts['G'], counts['C'], counts['T'], counts['Others'], counts['Gaps']])
            continue
        frequency = {}
        for char in column_no_gaps:
            if char in frequency:
                frequency[char] += 1
            else:
                frequency[char] = 1
        entropy = 0
        for freq in frequency.values():
            p = freq / len(column_no_gaps)
            entropy += -p * math.log(p, 2)
        entropy_list.append(entropy)
        position_info_list.append([counts['A'], counts['G'], counts['C'], counts['T'], counts['Others'], counts['Gaps']])

    return entropy_list, position_info_list

def main(fasta_file):
    sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    entropy_values, position_info = calculate_shannon_entropy(sequences)
    print("Position\tEntropy\tA\tG\tC\tT\tOthers\tGaps")
    for index, (entropy, counts) in enumerate(zip(entropy_values, position_info)):
        print(f"Position {index + 1}\t{entropy:.6f}\t{counts[0]}\t{counts[1]}\t{counts[2]}\t{counts[3]}\t{counts[4]}\t{counts[5]}")
main("aligned_V3-V4.fasta")
