NJOBS=16 ## number of threads/jobs to use/run in parallel
RAWFOLDER="./Data/F18FTSEUHT0029_Giovanni_Ece_0319/HUMvbwE"																	
STARGENOMENAME="mouse_release95"


# Step 1: user supplied the reference genome sequences (FASTA files) and annotations (GTF file), 
# from which STAR generate genome indexes that are utilized in the 2nd (map-ping) step. 
# The genome indexes are saved to disk and need only be generated once for each genome/annotation combination. 

# Make directory called stargenome 
# (location in which to place the genome generated for STAR which uses a non-compressed 
# suffix array - Reason for the memory requirement for STAR however contributes to its speed)
mkdir Data/stargenome
# go inside that directory
cd Data/stargenome
# generate a directory with the name specified at the start of the script
mkdir $STARGENOMENAME
# run star to generate the reference genome SA
STAR --runThreadN $NJOBS \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 33524399488 \
--genomeDir $STARGENOMENAME \
--genomeFastaFiles ../../Data/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa \
--sjdbGTFfile ../../Data/GRCm38/release95/Mus_musculus.GRCm38.95.gtf \
--sjdbOverhang 99						# Ideally, this length should be equal to the ReadLength-1
cd ../..