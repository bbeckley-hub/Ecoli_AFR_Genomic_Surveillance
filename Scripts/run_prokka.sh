#!/bin/bash

# Prokka Annotation Script 
# Uses 14 CPUs efficiently 
# Now processes .fna files and skips already completed annotations

# Configuration
CPUS=14
OUTPUT_DIR="prokka_annotations"
LOG_FILE="prokka_annotation.log"
SUMMARY_FILE="prokka_summary.tsv"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Start log
echo "PROKKA ANNOTATION STARTED: $(date)" | tee -a "$LOG_FILE"
echo "Available CPUs: $CPUS" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"

# Initialize summary file
echo -e "Sample\tStatus\tContigs\tCDS\trRNA\ttRNA\ttmRNA\tDate" > "$SUMMARY_FILE"

# Counter variables
total_files=0
processed=0
skipped=0
failed=0

# Function to check if Prokka already completed successfully
is_prokka_complete() {
    local sample_dir="$1"
    local sample_name="$2"
    [[ -f "$sample_dir/${sample_name}.gff" && -f "$sample_dir/${sample_name}.gbk" && -f "$sample_dir/${sample_name}.faa" ]]
}

# Function to get assembly stats
get_assembly_stats() {
    local fasta_file="$1"
    local contig_count=$(grep -c '^>' "$fasta_file" 2>/dev/null || echo "0")
    local total_length=$(grep -v '^>' "$fasta_file" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
    echo "$contig_count:$total_length"
}

# Process each .fna file
for fna_file in *.fna; do
    if [[ ! -f "$fna_file" ]]; then
        continue
    fi
    
    ((total_files++))
    
    # Extract sample name (remove .fna extension)
    sample_name="${fna_file%.fna}"
    sample_output="$OUTPUT_DIR/$sample_name"
    
    echo "Processing: $sample_name" | tee -a "$LOG_FILE"
    
    # Check if already processed
    if is_prokka_complete "$sample_output" "$sample_name"; then
        echo "  ✓ SKIPPING - Already completed" | tee -a "$LOG_FILE"
        echo -e "$sample_name\tSKIPPED\t-\t-\t-\t-\t-\t$(date +%Y-%m-%d)" >> "$SUMMARY_FILE"
        ((skipped++))
        continue
    fi
    
    # Get assembly stats before processing
    IFS=':' read -r contig_count total_length <<< "$(get_assembly_stats "$fna_file")"
    echo "  Assembly: $contig_count contigs, ${total_length} bp" | tee -a "$LOG_FILE"
    
    # Remove incomplete previous runs
    if [[ -d "$sample_output" ]]; then
        echo "  Removing incomplete previous run..." | tee -a "$LOG_FILE"
        rm -rf "$sample_output"
    fi
    
    # Run Prokka with optimized parameters (no --memory option)
    echo "  Starting Prokka annotation..." | tee -a "$LOG_FILE"
    
    # Set memory limits for Java tools that Prokka uses
    export _JAVA_OPTIONS="-Xmx14g"
    export PROKKA_TMPDIR="/tmp"
    
    # Run Prokka with optimized settings
    prokka \
        --outdir "$sample_output" \
        --prefix "$sample_name" \
        --cpus "$CPUS" \
        --compliant \
        --centre X \
        --force \
        --locustag "$sample_name" \
        --genus "Bacteria" \
        --species "species" \
        --strain "$sample_name" \
        --addgenes \
        --increment 10 \
        --gcode 11 \
        --mincontiglen 200 \
        "$fna_file" \
        2>&1 | sed 's/^/  [prokka] /' | tee -a "$LOG_FILE"
    
    # Check if Prokka completed successfully
    if [[ $? -eq 0 ]] && is_prokka_complete "$sample_output" "$sample_name"; then
        # Get annotation statistics
        if [[ -f "$sample_output/${sample_name}.gff" ]]; then
            cds_count=$(grep -c "CDS" "$sample_output/${sample_name}.gff" 2>/dev/null || echo "0")
            rrna_count=$(grep -c "rRNA" "$sample_output/${sample_name}.gff" 2>/dev/null || echo "0")
            trna_count=$(grep -c "tRNA" "$sample_output/${sample_name}.gff" 2>/dev/null || echo "0")
            tmrna_count=$(grep -c "tmRNA" "$sample_output/${sample_name}.gff" 2>/dev/null || echo "0")
        else
            cds_count=0; rrna_count=0; trna_count=0; tmrna_count=0
        fi
        
        echo "  ✓ SUCCESS: CDS=$cds_count, rRNA=$rrna_count, tRNA=$trna_count" | tee -a "$LOG_FILE"
        echo -e "$sample_name\tSUCCESS\t$contig_count\t$cds_count\t$rrna_count\t$trna_count\t$tmrna_count\t$(date +%Y-%m-%d)" >> "$SUMMARY_FILE"
        ((processed++))
    else
        echo "  ✗ FAILED: Prokka annotation failed" | tee -a "$LOG_FILE"
        echo -e "$sample_name\tFAILED\t$contig_count\t0\t0\t0\t0\t$(date +%Y-%m-%d)" >> "$SUMMARY_FILE"
        ((failed++))
        
        # Clean up failed run
        if [[ -d "$sample_output" ]]; then
            rm -rf "$sample_output"
        fi
    fi
    
    echo "  ---------------------------------------" | tee -a "$LOG_FILE"
    
    # Progress update every 5 samples
    if (( processed % 5 == 0 )); then
        echo "PROGRESS: $processed processed, $skipped skipped, $failed failed" | tee -a "$LOG_FILE"
    fi
done

# Final summary
echo "==========================================" | tee -a "$LOG_FILE"
echo "PROKKA ANNOTATION COMPLETED: $(date)" | tee -a "$LOG_FILE"
echo "TOTAL FILES: $total_files" | tee -a "$LOG_FILE"
echo "SUCCESSFULLY PROCESSED: $processed" | tee -a "$LOG_FILE"
echo "SKIPPED (already done): $skipped" | tee -a "$LOG_FILE"
echo "FAILED: $failed" | tee -a "$LOG_FILE"
echo "==========================================" | tee -a "$LOG_FILE"

# Generate final report
echo "" | tee -a "$LOG_FILE"
echo "Annotation summary saved to: $SUMMARY_FILE" | tee -a "$LOG_FILE"
echo "Detailed log saved to: $LOG_FILE" | tee -a "$LOG_FILE"

# Show summary table
echo "" | tee -a "$LOG_FILE"
echo "SUMMARY TABLE:" | tee -a "$LOG_FILE"
echo "==============" | tee -a "$LOG_FILE"
column -t -s $'\t' "$SUMMARY_FILE" | tee -a "$LOG_FILE"
