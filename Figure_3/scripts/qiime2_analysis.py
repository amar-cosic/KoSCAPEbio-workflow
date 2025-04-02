import argparse
import subprocess
import os
import sys
import qiime2
import shutil

def import_data(manifest_file, output_dir, seq_type):
    # QIIME 2 import command
    if seq_type == "paired-end":
        input_type = "SampleData[PairedEndSequencesWithQuality]"
        input_format = "PairedEndFastqManifestPhred33V2"
    elif seq_type == "single-end":
        input_type = "SampleData[SequencesWithQuality]"
        input_format = "SingleEndFastqManifestPhred33V2"
    else:
        sys.exit("Invalid sequence type. Use 'paired-end' or 'single-end'.")

    import_command = f"qiime tools import \
                    --type '{input_type}' \
                    --input-path {manifest_file} \
                    --output-path {output_dir}/demux.qza \
                    --input-format {input_format}"

    subprocess.run(import_command, shell=True, check=True)

def summarize_data(input_demux, output_dir):
    summary_command = f"qiime demux summarize \
                    --i-data {input_demux} \
                    --o-visualization {output_dir}/demux-summary.qzv"

    subprocess.run(summary_command, shell=True, check=True)


def denoise_data(input_demux, output_dir, seq_type, trim_left, len_left, trim_right=None, len_right=None):
    if seq_type == "single-end":
        denoise_command = f"""
            qiime dada2 denoise-single \
                --i-demultiplexed-seqs {input_demux} \
                --p-trim-left {trim_left} \
                --p-trunc-len {len_left} \
                --o-representative-sequences {output_dir}/rep-seqs.qza \
                --o-table {output_dir}/table.qza \
                --o-denoising-stats {output_dir}/denoising-stats.qza
        """
    elif seq_type == "paired-end":  # paired-end
        denoise_command = f"""
            qiime dada2 denoise-paired \
                --i-demultiplexed-seqs {input_demux} \
                --p-trim-left-f {trim_left} \
                --p-trim-left-r {trim_right} \
                --p-trunc-len-f {len_left} \
                --p-trunc-len-r {len_right} \
                --o-representative-sequences {output_dir}/rep-seqs.qza \
                --o-table {output_dir}/table.qza \
                --o-denoising-stats {output_dir}/denoising-stats.qza
        """
    else:
        sys.exit("Invalid sequence type. Use 'paired-end' or 'single-end'.")

    subprocess.run(denoise_command, shell=True, check=True)


def visualize_denoising_stats(input_stats, output_dir):
    stats_command = f"qiime metadata tabulate \
                     --m-input-file {input_stats} \
                     --o-visualization {output_dir}/denoising-stats.qzv"

    subprocess.run(stats_command, shell=True, check=True)

def export_table_and_rep_seqs(input_table, input_rep_seqs, output_dir):
    table_export_dir = os.path.join(output_dir, 'exported_table')
    table = qiime2.Artifact.load(input_table)
    table.export_data(table_export_dir)

    biom_file = os.path.join(table_export_dir, 'feature-table.biom')
    tsv_file = os.path.join(output_dir, 'table.tsv')
    convert_command = f"biom convert -i {biom_file} -o {tsv_file} --to-tsv"
    
    convert_result = subprocess.run(convert_command, shell=True)
    
    if convert_result.returncode == 0:
        print("Successfully exported and converted feature table to TSV.")
    else:
        print("Error: Failed to export and convert feature table to TSV.")

    # Clean up the export folder
    shutil.rmtree(table_export_dir)

    # Export representative sequences to FASTA format
    rep_seqs = qiime2.Artifact.load(input_rep_seqs)
    rep_seq_dir=os.path.join(output_dir, 'rep-seqs')
    rep_seqs.export_data(os.path.join(output_dir, 'rep-seqs'))

    rename_command=f"mv {rep_seq_dir}/dna-sequences.fasta {rep_seq_dir}/rep_seq.fasta"
    rename_result=subprocess.run(rename_command,shell=True)
    move_command=f"mv {rep_seq_dir}/rep_seq.fasta {output_dir}"
    
    move_result=subprocess.run(move_command,shell=True)
    if move_result.returncode ==0 and rename_result.returncode==0:
        print("Successfully exported representative sequences to FASTA.")
    else:
        print("Error: Failed to export representative sequences to FASTA.")

    shutil.rmtree(rep_seq_dir)

    if convert_result.returncode == 0:
        print("QIIME 2 data export and conversion completed successfully.")
    else:
        print("QIIME 2 data export and conversion encountered errors.")

def classify_sequences(classifier, input_reads, output_classification):
    classify_command = f"""
    qiime feature-classifier classify-sklearn \
    --i-classifier {classifier} \
    --i-reads {input_reads} \
    --o-classification {output_classification}
    """
    subprocess.run(classify_command, shell=True, check=True)
    print("Sequence classification completed.")

def tabulate_taxonomy(input_taxonomy, output_visualization):
    tabulate_command = f"""
    qiime metadata tabulate \
    --m-input-file {input_taxonomy} \
    --o-visualization {output_visualization}
    """
    subprocess.run(tabulate_command, shell=True, check=True)
    print("Taxonomy tabulation completed.")

def create_taxa_bar_plots(table, taxonomy, metadata, output_visualization):
    bar_plot_command = f"""
    qiime taxa barplot \
    --i-table {table} \
    --i-taxonomy {taxonomy} \
    --m-metadata-file {metadata} \
    --o-visualization {output_visualization}
    """
    subprocess.run(bar_plot_command, shell=True, check=True)
    print("Taxa bar plot visualization completed.")



def main():
    parser = argparse.ArgumentParser(description="Run QIIME2 Analysis")
    parser.add_argument("-manifest", required=True, help="Path to the manifest file")
    parser.add_argument("-seq_type", choices=["paired-end", "single-end"], required=True, help="Sequence type")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-skip_import", action='store_false', default=True, help="Skip the import step")
    parser.add_argument("-trim_left", "--trim_left", type=int, default=40, help="Position at which forward read sequences should be trimmed (default: 40)")
    parser.add_argument("-len_left", "--len_left", type=int, default=244, help="Position at which forward read sequences should be truncated (default: 244)")
    parser.add_argument("-trim_right", "--trim_right", type=int, default=40, help="Position at which reverse read sequences should be trimmed (default: 40, paired-end only)")
    parser.add_argument("-len_right", "--len_right", type=int, default=244, help="Position at which reverse read sequences should be truncated (default: 244, paired-end only)")
    parser.add_argument("-classifier", required=True, help="Path to the Qiime classifier")
    parser.add_argument("-metadata", required=False, help="Path to the metadata file for taxa barplot")


    args = parser.parse_args()


    # Create the output directory if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    input_demux = os.path.join(args.output, "demux.qza")

    if args.skip_import:
        import_data(args.manifest, args.output, args.seq_type)

        summarize_data(input_demux, args.output)

    denoise_data(input_demux, args.output, args.seq_type, args.trim_left, args.len_left, args.trim_right, args.len_right)

    input_stats = os.path.join(args.output, "denoising-stats.qza")
    visualize_denoising_stats(input_stats, args.output)

    export_table_and_rep_seqs(os.path.join(args.output, "table.qza"),
                              os.path.join(args.output, "rep-seqs.qza"),
                              args.output)

    classify_sequences(args.classifier, os.path.join(args.output, "rep-seqs.qza"), os.path.join(args.output, "taxonomy.qza"))

    tabulate_taxonomy(os.path.join(args.output, "taxonomy.qza"), os.path.join(args.output, "taxonomy.qzv"))

    create_taxa_bar_plots(os.path.join(args.output, "table.qza"), os.path.join(args.output, "taxonomy.qza"), args.metadata, os.path.join(args.output, "taxa-bar-plot.qzv"))

    print("QIIME 2 data analysis pipeline completed.")

if __name__ == "__main__":
    main()
