#!/bin/bash
cd ~/Africa/genomes/Sts/busco_results/

rm -f summary.tsv
echo -e "genome\tC(%)\tS(%)\tD(%)\tF(%)\tM(%)\tscaffolds\tcontigs\ttotal_length\tscaffold_N50\tcontig_N50" > summary_busco.tsv

for dir in *_busco/; do
    genome="${dir%_busco/}"
    summary_file="${dir}short_summary.specific.enterobacterales_odb12.${genome}_busco.txt"
    [ ! -f "$summary_file" ] && continue

    # Use grep without ^ to catch leading whitespace
    line=$(grep "C:" "$summary_file" | head -1)

    C=$(echo "$line" | sed -n 's/.*C:\([0-9.]*\)%.*/\1/p')
    S=$(echo "$line" | sed -n 's/.*S:\([0-9.]*\)%.*/\1/p')
    D=$(echo "$line" | sed -n 's/.*D:\([0-9.]*\)%.*/\1/p')
    F=$(echo "$line" | sed -n 's/.*F:\([0-9.]*\)%.*/\1/p')
    M=$(echo "$line" | sed -n 's/.*M:\([0-9.]*\)%.*/\1/p')

    scaffolds=$(grep "Number of scaffolds" "$summary_file" | awk '{print $1}')
    contigs=$(grep "Number of contigs" "$summary_file" | awk '{print $1}')
    total_length=$(grep "Total length" "$summary_file" | awk '{print $1}')
    scaffold_N50=$(grep "Scaffold N50" "$summary_file" | awk '{print $1}')
    contig_N50=$(grep "Contigs N50" "$summary_file" | awk '{print $1}')

    echo -e "${genome}\t${C}\t${S}\t${D}\t${F}\t${M}\t${scaffolds}\t${contigs}\t${total_length}\t${scaffold_N50}\t${contig_N50}" >> summary_busco.tsv
done

echo "Done! Check summary_busco.tsv"
