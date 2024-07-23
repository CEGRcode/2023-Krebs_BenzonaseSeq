# Depends on download from PEGR script alias (`symlink-bamfiles`)
#  check with Olivia to configure this for your computer,
#  otherwise call the python script yourself by editing this script directly

### CHANGE ME
MYSAMPLES=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing/mysamples.txt
BAMDIR=/path/to/2023-Krebs_BenzonaseSeq/data/BAM
###

cd $BAMDIR
symlink-bamfiles $MYSAMPLES
