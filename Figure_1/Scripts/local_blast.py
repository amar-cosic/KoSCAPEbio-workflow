import os
import subprocess
import pandas as pd
def prepare_db(db_path):
    cmd = f"makeblastdb -in {db_path} -dbtype nucl -out blast_db/custom_db"
    subprocess.run(cmd, shell=True, check=True)
def run_blast(query, db, output):
    cmd = f"blastn -query {query} -db {db} -out {output} -outfmt 6"
    subprocess.run(cmd, shell=True, check=True)
def parse_blast_file(file_path):
    blast_data = {}
    with open(file_path, 'r') as blast_file:
        for line in blast_file:
            qseqid, sseqid, pident, length = line.strip().split('\t')[:4]
            length = int(length)
            pident = float(pident)
            blast_data[sseqid] = {'genome': qseqid, 'pident': pident, 'length': length}
    return blast_data
def parse_fasta_headers(fasta_path):
    headers = []
    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                header = line.strip().split(None, 1)[0][1:]
                headers.append(header)
    return headers
def process_blast_directory(blast_dir):
    all_data = {}
    for filename in os.listdir(blast_dir):
        if filename.endswith('.tsv'):
            file_path = os.path.join(blast_dir, filename)
            file_data = parse_blast_file(file_path)
            genome_key = filename.replace('_vs_custom_db.tsv', '')
            all_data[genome_key] = file_data
    return all_data
def create_summary_table(genomes_dir, blast_data, headers, output_file):
    genome_names = {f.replace('.fasta', ''): f for f in os.listdir(genomes_dir) if f.endswith('.fasta')}
    column_names = ['sample_name'] + headers
    rows = []
    for genome_key, genome_file in genome_names.items():
        row = [genome_file] + [None] * len(headers)
        genome_results = blast_data.get(genome_key, {})
        
        for header in headers:
            if header in genome_results:
                row[headers.index(header) + 1] = genome_results[header]['pident']
        rows.append(row)
    df = pd.DataFrame(rows, columns=column_names)
    df.to_csv(output_file, sep='\t', index=False)
def main():
    db_path = 'blast_db/custom_db.fasta'
    genomes_dir = 'oxy_test/Genomes/fasta'
    blast_dir = 'blast_results'
    fasta_path = 'blast_db/custom_db.fasta'
    output_file = 'combined_blast_results.tsv'
    prepare_db(db_path)
    os.makedirs(blast_dir, exist_ok=True)
    for filename in os.listdir(genomes_dir):
        if filename.endswith('.fasta'):
            genome_path = os.path.join(genomes_dir, filename)
            output_path = os.path.join(blast_dir, f"{os.path.splitext(filename)[0]}_vs_custom_db.tsv")
            run_blast(genome_path, 'blast_db/custom_db', output_path)
    blast_data = process_blast_directory(blast_dir)
    headers = parse_fasta_headers(fasta_path)
    create_summary_table(genomes_dir, blast_data, headers, output_file)
if __name__ == "__main__":
    main()
