cd ~/Downloads/PHYLO_A_B_C

for f in Phylo-*.pdf; do
  base="${f%.pdf}"

  # PNG at 600 dpi
  gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r600 \
     -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
     -sOutputFile="${base}.png" "$f"

  # TIFF at 600 dpi with LZW compression
  gs -dNOPAUSE -dBATCH -sDEVICE=tiff24nc -r600 \
     -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
     -sCompression=lzw \
     -sOutputFile="${base}.tiff" "$f"

done

ls -lh
