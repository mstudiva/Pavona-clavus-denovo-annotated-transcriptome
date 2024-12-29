## De novo transcriptome assembly pipeline for Pavona clavus, version December 29, 2024
# Adapted by Michael Studivan (studivanms@gmail.com) based on repos by Misha Matz (https://github.com/z0on/annotatingTranscriptomes.git), Eli Meyer (https://github.com/Eli-Meyer/sequence_utilities.git; https://github.com/Eli-Meyer/transcriptome_utilities.git), and  Brian Strehlow (https://github.com/bstrehlow/Transcriptome_assembly.git) for use on the FAU KoKo HPC


#------------------------------
## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

# log onto cluster
ssh mstudiva@koko-login.hpc.fau.edu


#------------------------------
## Installing scripts and setting up the workspace

# switch to home directory
cd

# unless you have done it in the past, make directory called bin,
# all your scripts should go in there:
mkdir bin

# switch to bin:
cd bin

# clone github repositories
git clone https://github.com/mstudiva/Pavona-clavus-denovo-annotated-transcriptome.git
git clone https://github.com/mstudiva/annotatingTranscriptomes.git

# move files from subdirectories to bin/:
mv Pavona-clavus-denovo-annotated-transcriptome/* .
mv annotatingTranscriptomes/* .

chmod +x ~/bin/bs

# Installing Trinity
cd ~/bin/
wget https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg


#------------------------------
## Transcriptome assembly using Trinity

# Combine all good, trimmed reads into a single fastq
srun cat *.trim > Pclavus_reads.fastq
# Move this to ~/annotate/

# Trinity
echo '#!/bin/bash' > trinity.sh
echo '#SBATCH --partition=mediumq7' >> trinity.sh
echo '#SBATCH -N 1' >> trinity.sh
echo '#SBATCH --exclusive' >> trinity.sh
echo '#SBATCH --mem=200GB' >> trinity.sh
echo 'module load singularity/3.4.1' >> trinity.sh
echo 'singularity exec -e ~/bin/trinityrnaseq.v2.15.2.simg Trinity --seqType fq --single `pwd`/Pclavus_reads.fastq --CPU 20 --max_memory 200G --output `pwd`/Pclavus_trinity' >> trinity.sh
sbatch -o trinity.o%j -e trinity.e%j trinity.sh

# Or upload the concatenated reads to Galaxy using their Trinity pipeline: https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Ftrinity%2Ftrinity%2F2.15.1%2Bgalaxy1&version=latest

mv Pclavus_trinity.Trinity.fasta Pclavus_Trinity.fasta
mv Galaxy7-[Trinity_on_data_6__Assembled_Transcripts].fasta Pclavus_Galaxy.fasta
mv *.fasta ..

conda activate bioperl
echo "seq_stats.pl Pclavus_Trinity.fasta > seqstats_Pclavus_Trinity.txt" > seq_stats
echo "seq_stats.pl Pclavus_Galaxy.fasta > seqstats_Pclavus_Galaxy.txt" >> seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Pclavus_Trinity.fasta
-------------------------
250657 sequences.
346 average length.
12396 maximum length.
180 minimum length.
N50 = 334
86.7 Mb altogether (86740838 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

Pclavus_Galaxy.fasta
-------------------------
250614 sequences.
346 average length.
12396 maximum length.
177 minimum length.
N50 = 334
86.7 Mb altogether (86737086 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## First cleaning step:
# Assemblies include many small contigs that are unlikely to provide significant matches, so for analyses based on sequence homology we consider only contigs ≥300 bp.

conda activate bioperl
echo "perl ~/bin/noshorts.pl Pclavus_Trinity.fasta 300" > noshorts
echo "perl ~/bin/noshorts.pl Pclavus_Galaxy.fasta 300" >> noshorts
launcher_creator.py -j noshorts -n noshorts -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch noshorts.slurm

echo "seq_stats.pl Pclavus_Trinity_noshorts.fasta > seqstats_Pclavus_Trinity_noshorts.txt" > seq_stats2
echo "seq_stats.pl Pclavus_Galaxy_noshorts.fasta > seqstats_Pclavus_Galaxy_noshorts.txt" >> seq_stats2
launcher_creator.py -j seq_stats2 -n seq_stats2 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats2.slurm

Pclavus_Trinity_noshorts.fasta
-------------------------
124595 sequences.
445 average length.
12396 maximum length.
300 minimum length.
N50 = 422
55.4 Mb altogether (55428670 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

Pclavus_Galaxy_noshorts.fasta
-------------------------
124645 sequences.
445 average length.
12396 maximum length.
300 minimum length.
N50 = 422
55.5 Mb altogether (55454713 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Second cleaning step:
Remove ribosomal/mitochondrial contaminants using SILVA and NCBI databases and available references for the target species

# microbial ribosomal RNA
# Go to https://www.arb-silva.de/no_cache/download/archive/current/Exports/ and copy the links for the LSURef_NR99 and SSURef_NR99 compressed fasta files
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
gunzip *.gz
# Concatenate both into a single fasta file for blasting
cat SILVA_138.2_LSURef_NR99_tax_silva.fasta SILVA_138.2_SSURef_NR99_tax_silva.fasta > SILVA_SSU_LSU_combined.fasta

# eukaryotic ribosomal RNA
# Go to https://www.arb-silva.de/search/ and search for your species; if it's not available, pick a related species
# Add it to your cart using the checkbox on the left, then select Download in the top right
# Choose FASTA without gaps, and tar.gz
cp ~/bin/arb-silva.de_2024-12-27_id1369219.tgz .
tar -vxf arb-silva.de_2024-12-27_id1369219.tgz

# concatenate ribosomal RNA assemblies
cat SILVA_SSU_LSU_combined.fasta arb-silva.de_2024-12-27_id1369219_tax_silva_trunc.fasta > Pclavus_rRNA.fasta

# mitochondrial RNA
# https://www.ncbi.nlm.nih.gov/nuccore/NC_008165.1
cp ~/bin/Pclavus_mitoRNA.fasta .

conda activate bioperl
# Running a perl script that blasts the transcriptome against inputted contamination sequences to generate a final, clean transcriptome
echo '#!/bin/bash' > contam.sh
echo '#SBATCH --partition=mediumq7' >> contam.sh
echo '#SBATCH -N 1' >> contam.sh
echo '#SBATCH --mem=200GB' >> contam.sh
echo 'conda activate bioperl' >> contam.sh
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> contam.sh
echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=Pclavus_Trinity_noshorts.fasta contam=rRNA,Pclavus_rRNA.fasta contam=Mt,Pclavus_mitoRNA.fasta table=Pclavus_Trinity_contamination.txt passed=Pclavus_Trinity_final.fasta" >> contam.sh
sbatch -o contam.o%j -e contam.e%j contam.sh

echo 'conda activate bioperl' > contam2.sh
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> contam2.sh
echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=Pclavus_Galaxy_noshorts.fasta contam=rRNA,Pclavus_rRNA.fasta contam=Mt,Pclavus_mitoRNA.fasta table=Pclavus_Galaxy_contamination.txt passed=Pclavus_Galaxy_final.fasta" >> contam2
launcher_creator.py -j contam2 -n contam2 -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch --mem=200GB contam2.slurm

Pclavus_Trinity_final.fasta
-------------------------
124595 sequences in input file
1256 sequences look like contaminants
        rRNA    967
        Mt	289
123339 sequences passed all tests
-------------------------

Pclavus_Galaxy_final.fasta
-------------------------
124645 sequences in input file
268 sequences look like contaminants
        Mt	268
        rRNA    0
124645 sequences passed all tests
-------------------------

echo "seq_stats.pl Pclavus_Trinity_final.fasta > seqstats_Pclavus_Trinity_final.txt" > seq_stats3
echo "seq_stats.pl Pclavus_Galaxy_final.fasta > seqstats_Pclavus_Galaxy_final.txt" >> seq_stats3
launcher_creator.py -j seq_stats3 -n seq_stats3 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats3.slurm


#------------------------------
## Final assembly!

Pclavus_Trinity_final.fasta
-------------------------
123339 sequences.
445 average length.
12396 maximum length.
300 minimum length.
N50 = 421
54.8 Mb altogether (54842903 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

Pclavus_Galaxy_final.fasta
-------------------------
124645 sequences.
445 average length.
12396 maximum length.
300 minimum length.
N50 = 422
55.5 Mb altogether (55454713 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

# Now proceed with separation of host/symbiont contigs in the script 'Pclavus_transcriptome_hostsym_separation_README'
