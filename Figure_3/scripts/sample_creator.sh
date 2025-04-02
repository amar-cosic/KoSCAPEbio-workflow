#!/bin/bash

if [ "$1" == "-s" ]; then
    MODE="single_end"
    READ_MODE="" 
    echo "Running in Single-End mode."
elif [ "$1" == "-p" ]; then
    MODE="paired_end"
    READ_MODE="-p"
    echo "Running in Paired-End mode."
else
    echo "Usage: $0 [-s | -p]"
    echo "-s: Run in Single-End mode"
    echo "-p: Run in Paired-End mode"
    exit 1
fi

BASE_DIR="mock_fastq"

MODE_DIR="$BASE_DIR/$MODE"

declare -a SUB_DIRS=("V3_V4" "V4" "V4_V5" "V3_V5")

if [ ! -d "$BASE_DIR" ]; then
    echo "Creating base directory: $BASE_DIR"
    mkdir "$BASE_DIR"
fi

if [ ! -d "$MODE_DIR" ];then
    echo "Creating mode-specific directory: $MODE_DIR"
    mkdir "$MODE_DIR"
fi

for subdir in "${SUB_DIRS[@]}"; do
    FULL_PATH="$MODE_DIR/$subdir"
    if [ ! -d "$FULL_PATH" ]; then
        echo "Creating subdirectory: $FULL_PATH"
        mkdir "$FULL_PATH"
    else
        echo "Subdirectory already exists: $FULL_PATH"
    fi
done

echo "Directory structure setup completed."


NUM_SAMPLES=5

for i in $(seq 1 $NUM_SAMPLES); do
    for subdir in "${SUB_DIRS[@]}"; do
        FULL_PATH="$MODE_DIR/$subdir"
        INPUT_FILE="hypervariable_regions/${subdir}.fasta"
        
        if [ ! -f "$INPUT_FILE" ]; then
            echo "Error: Input file $INPUT_FILE does not exist."
            exit 1
        fi

        echo "Generating reads for sample $i in $subdir"
        art_illumina -ss MSv3 -i "$INPUT_FILE" $READ_MODE -l 250 -f 200 -m 600 -s 10 -o "${FULL_PATH}/sample-${i}_"
    done
done


for mode_subdir in "$BASE_DIR"/*; do
    if [ -d "$mode_subdir" ]; then
        for subdir in "$mode_subdir"/*; do
            if [ -d "$subdir" ]; then
                echo "Renaming in $subdir"
                for file in "$subdir"/*.fq; do
                    if [ -f "$file" ]; then
                        mv "$file" "${file%.fq}.fastq"
                    fi
                done
            fi
        done
    fi
done

echo "Renaming completed."
