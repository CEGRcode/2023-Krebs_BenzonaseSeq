PARSE=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/220802_NucCall/job/convert_scIDX_to_strandless.pl
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc

mkdir -p $OUTPUT

CPU=0

for file in $INPUT/*.tab; do
	var=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F"." '{print $1}')
        set -- $var
        echo $1

	perl $PARSE $file $OUTPUT/$1_NUC.tab &

        # Multi-thread to 4 cores
        let CPU++
        if [[ $CPU -eq 4 ]]; then
                wait
                CPU=0
        fi


done
