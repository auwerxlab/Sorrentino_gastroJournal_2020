NJOBS=16 ## number of threads/jobs to use/run in parallel
RAWFOLDER="./Data/F18FTSEUHT0029_Giovanni_Ece_0319/HUMvbwE"																		
STARGENOMENAME="mouse_release95"

# Step 2: Quality Control
# directory to contain the output of the QC
mkdir Data/QC
# SAMPLEID contains the name of each of the samples for which to run QC
# to get these, 1- Find the names in the raw data folder, 2- pipe it into xargs to process each input
# 3- for each input use tr which deletes given -d the given string from the input 
SAMPLEID="$(find $RAWFOLDER/*fq.gz | xargs -n 1 basename | tr -d "fq.gz")"
# loop for each sample ID and perform QC
for SAMPLE in $SAMPLEID
do
	SAMPLENAME="$SAMPLE.fq.gz"
	echo "DOING QC FOR $SAMPLE AND $SAMPLENAME"
	# generate a folder for this sample
	mkdir ./Data/QC/$SAMPLE
	# fastqc input file -0 output file name
	/usr/local/HTS/FastQC/fastqc $RAWFOLDER/$SAMPLENAME -o ./Data/QC/$SAMPLE                   	
	# Check if this is the fastqc directory on my project location
done