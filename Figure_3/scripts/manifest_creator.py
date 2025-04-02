import os
import sys

def create_manifest(input_folder, paired_end=False, type_manifest_creation=False):
    manifest_content = "sample-id\t"
    if paired_end:
        manifest_content += "forward-absolute-filepath\treverse-absolute-filepath\n"
    else:
        manifest_content += "absolute-filepath\n"
    file_paths = {}

    for filename in os.listdir(input_folder):
        full_path = os.path.abspath(os.path.join(input_folder, filename))

        if filename.endswith('.fastq'):
            if paired_end:
                parts = filename.split('_')
                base_name = '_'.join(parts[:-1])
                read_direction = parts[-1][0]

                if base_name not in file_paths:
                    file_paths[base_name] = {'1': '', '2': ''}
                file_paths[base_name][read_direction] = full_path
            else:
                base_name = filename.split('.fastq')[0]
                file_paths[base_name] = full_path
        elif type_manifest_creation and filename.endswith('.txt'):
            base_name = filename.split('__')[0]
            file_paths[base_name] = full_path

    for sample_id, paths in file_paths.items():
        if paired_end:
            manifest_content += f"{sample_id}\t{paths['1']}\t{paths['2']}\n"
        else:
            manifest_content += f"{sample_id}\t{paths}\n"

    return manifest_content

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("Usage: python script.py /absolute/path/to/folder [-s or -p]. -s for single_end / -p for paired_end. If you want manifest for srst2 type -sr as the third argument. Example script.py /absolute/path/to/folder -s -sr")
        sys.exit(1)

    input_folder = sys.argv[1]
    paired_end = "-p" in sys.argv
    type_manifest_creation = "-sr" in sys.argv

    manifest_content = create_manifest(input_folder, paired_end, type_manifest_creation)

    output_file = os.path.join(input_folder, 'manifest.tsv')
    with open(output_file, 'w') as file:
        file.write(manifest_content)

    print(f"Manifest file created at {output_file}")
