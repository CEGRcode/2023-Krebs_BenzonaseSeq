
GETENCFF=../bin/get_ENCFF-Peak_from_ENCODEsearch.py

# URL filters:
# Remove "&format=json&limit=all" to view search on ENCODE portal
SEARCH_URL="https://www.encodeproject.org/search/?type=File&file_format_type=narrowPeak&file_type=bed+narrowPeak&assay_title=TF+ChIP-seq&biosample_ontology.term_name=K562&biosample_ontology.term_name=A549&biosample_ontology.term_name=MCF-7&biosample_ontology.term_name=HepG2&status=released&output_type=conservative+IDR+thresholded+peaks&assembly=hg19&format=json&limit=all"
OUTPUT=hg19_peak_encff.tsv
python $GETENCFF -i $SEARCH_URL -o $OUTPUT
