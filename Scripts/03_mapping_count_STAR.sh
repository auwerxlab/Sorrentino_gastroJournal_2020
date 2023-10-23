NJOBS=16 ## number of threads/jobs to use/run in parallel
RAWFOLDER="./Data/F18FTSEUHT0029_Giovanni_Ece_0319/HUMvbwE"                                                                                                                                             
STARGENOMENAME="mouse_release95"

# Step 3: Mapping to reference genome
# Gene Levels + Transcript Levels + Gene Counts
mkdir ./Data/mapped
SAMPLEID="$(find $RAWFOLDER/*_1.fq.gz | xargs -n 1 basename | tr -d "fq.gz" )"
for SAMPLE in $SAMPLEID
do
		# removes the last occurence of '%_1' from SAMPLE name
        SAMPLE=${SAMPLE%_1}
        echo "Mapping sample $SAMPLE"
        # Ask him why he is doing this
        R1=(`ls -d -1 $RAWFOLDER/$SAMPLE"_1.fq.gz"`)
        R2=(`ls -d -1 $RAWFOLDER/$SAMPLE"_2.fq.gz"`)
        STAR --runThreadN $NJOBS \
        --genomeDir ./Data/stargenome/$STARGENOMENAME \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesIn $R1 $R2 \
        --outFileNamePrefix ./Data/mapped/$SAMPLE
done