#!/bin/bash

# =============================================================================
# Download African E. coli genomes using NCBI datasets CLI
# - Skips already downloaded genomes
# - Retries up to 3 times on failure
# - Extracts FASTA (.fna) to a single folder
# =============================================================================

# Input TSV file (first column = Assembly accession)
INPUT_TSV="ecoli_clean.tsv"

# Output directory for FASTA files
GENOME_DIR="genomes"

# Temporary download directory
TMP_DIR="tmp_downloads"

# Log file
LOG_FILE="download.log"

# Number of retries per genome
MAX_RETRIES=3

# Delay between retries (seconds)
RETRY_DELAY=10

# Delay between different genomes (be kind to NCBI)
GENOME_DELAY=2

# Create output directories
mkdir -p "$GENOME_DIR" "$TMP_DIR"

# Clear log
> "$LOG_FILE"

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Count total genomes
total=$(tail -n +2 "$INPUT_TSV" | wc -l)
current=0

log "Starting download of $total genomes..."

# Loop through all lines except header
tail -n +2 "$INPUT_TSV" | while IFS=$'\t' read -r acc biosample taxonomy bioproject country; do
    # Trim whitespace
    acc=$(echo "$acc" | xargs)
    if [ -z "$acc" ]; then
        continue
    fi

    current=$((current + 1))
    target_fasta="${GENOME_DIR}/${acc}.fna"

    # Skip if already downloaded
    if [ -f "$target_fasta" ]; then
        log "[$current/$total] SKIP: $acc (already exists)"
        continue
    fi

    log "[$current/$total] DOWNLOAD: $acc (country: $country)"

    attempt=1
    success=0
    while [ $attempt -le $MAX_RETRIES ] && [ $success -eq 0 ]; do
        # Temporary zip file
        zip_file="${TMP_DIR}/${acc}.zip"

        # Download genome (only genome FASTA)
        datasets download genome accession "$acc" --include genome --filename "$zip_file" >> "$LOG_FILE" 2>&1

        # Check if download succeeded and zip is valid
        if [ -f "$zip_file" ] && unzip -t "$zip_file" > /dev/null 2>&1; then
            # Extract the .fna file(s)
            # The zip contains a folder structure like ncbi_dataset/data/GCA_xxx/GCA_xxx.fna
            # We'll find the .fna file and move it to the target location
            fna_file=$(unzip -l "$zip_file" | grep -E '\.fna$' | awk '{print $NF}' | head -1)
            if [ -n "$fna_file" ]; then
                unzip -q "$zip_file" "$fna_file" -d "$TMP_DIR"
                # Move it to the genome directory
                mv "${TMP_DIR}/${fna_file}" "$target_fasta"
                log "  SUCCESS: $acc -> $target_fasta"
                success=1
            else
                log "  ERROR: No .fna file found in zip for $acc"
            fi
            # Clean up zip
            rm -f "$zip_file"
        else
            log "  ATTEMPT $attempt FAILED for $acc"
            rm -f "$zip_file"
            sleep "$RETRY_DELAY"
        fi

        attempt=$((attempt + 1))
    done

    if [ $success -eq 0 ]; then
        log "  FATAL: All $MAX_RETRIES attempts failed for $acc"
        # Optionally, record failed accession for later retry
        echo "$acc" >> failed_accessions.txt
    fi

    # Be polite to NCBI between genomes
    sleep "$GENOME_DELAY"
done

log "Download process completed. Check $LOG_FILE for details."
log "Failed accessions (if any) are in failed_accessions.txt"
