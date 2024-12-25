## De novo transcriptome assembly pipeline for Pavona clavus, version December 24, 2024 (Christmas Eve edition!)
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
git clone https://github.com/mstudiva/SCTLD-Panama-omics.git
git clone https://github.com/mstudiva/annotatingTranscriptomes.git
git clone https://github.com/mstudiva/tag-based_RNAseq.git
git clone https://github.com/hrivera28/Oculina_arbuscula_transcriptome.git

# move files from subdirectories to bin/:
mv SCTLD-Panama-omics/* .
mv annotatingTranscriptomes/* .
mv tag-based_RNAseq/* .
mv Oculina_arbuscula_transcriptome/Sarahs_scripts/* .

chmod +x ~/bin/bs

# Installing Trinity
cd ~/bin/
wget https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/trinityrnaseq.v2.15.2.simg


#------------------------------
## Transcriptome assembly using Trinity

# Combine all good reads into a single fastq
srun cat *.trim > Pclavus_reads.fastq
# Move this to ~/annotate/

# Trinity
echo '#!/bin/bash' > trinity.sh
echo '#SBATCH --partition=longq7' >> trinity.sh
echo '#SBATCH -N 1' >> trinity.sh
echo '#SBATCH --exclusive' >> trinity.sh
echo '#SBATCH --mem=200GB' >> trinity.sh
echo 'module load singularity/3.4.1' >> trinity.sh
echo 'singularity exec -e ~/bin/trinityrnaseq.v2.15.2.simg Trinity --seqType fq --single `pwd`/Pclavus_reads.fastq --CPU 20 --max_memory 200G --output `pwd`/Pclavus_trinity' >> trinity.sh
sbatch -o trinity.o%j -e trinity.e%j trinity.sh

# If any of the assemblies fail in the chrysalis step, find the output directory for each of the error files and delete them, or move them to your backup directory. They should look like this: "Pclavus_trinity/read_partitions/Fb_4/CBin_4670/c467359.trinity.reads.fa.out"
mv Pclavus_trinity/read_partitions/Fb_4/CBin_4670/c467359.trinity.reads.fa.out Pclavus_trinity/read_partitions/Fb_3/CBin_3088/c309109.trinity.reads.fa.out Pclavus_trinity/read_partitions/Fb_3/CBin_3690/c369317.trinity.reads.fa.out Pclavus_trinity/read_partitions/Fb_3/CBin_3701/c370414.trinity.reads.fa.out temp_backup/.

mv Pclavus_trinity.Trinity.fasta Cvarians_Trinity.fasta
echo "seq_stats.pl Cvarians_Trinity.fasta > seqstats_Cvarians_Trinity.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Cvarians_Trinity.fasta
-------------------------
812356 sequences.
801 average length.
26720 maximum length.
182 minimum length.
N50 = 1343
650.6 Mb altogether (650572316 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## First cleaning step:
# Assemblies include many small contigs that are unlikely to provide significant matches, so for analyses based on sequence homology we consider only contigs ≥500 bp.

echo "perl ~/bin/noshorts.pl Cvarians_Trinity.fasta 500" > noshorts
launcher_creator.py -j noshorts -n noshorts -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch noshorts.slurm

echo "seq_stats.pl noshorts_Cvarians_Trinity.fasta > seqstats_noshorts_Cvarians_Trinity.txt" > seq_stats2
launcher_creator.py -j seq_stats2 -n seq_stats2 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats2.slurm

noshorts_Cvarians_Trinity.fasta
-------------------------
346223 sequences.
1453 average length.
26720 maximum length.
500 minimum length.
N50 = 1873
503.1 Mb altogether (503073528 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Second cleaning step:
Remove ribosomal/mitochondrial contaminants using SILVA and NCBI databases and available references for the target species

# ribosomal RNA, including microbial sequences
# Go to https://www.arb-silva.de/no_cache/download/archive/current/Exports/ and copy the links for the LSURef_NR99 and SSURef_NR99 compressed fasta files
# wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
# wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
# gunzip *.gz
# Concatenate both into a single fasta file for blasting
# cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVA_SSU_LSU_combined.fasta

# Go to https://www.arb-silva.de/search/ and search for your species; if it's not available, pick a related species
# Add it to your cart using the checkbox on the left, then select Download in the top right
# Choose FASTA without gaps, and tar.gz
# Pclavus varians is available in this GitHub repo as 'arb-silva.de_2021-12-01_id1089989'
# Originally from Redmond et al. (2013) doi: 10.1093/icb/ict078
cp ~/bin/arb-silva.de_2021-12-01_id1089989.tgz .
tar -vxf arb-silva.de_2021-12-01_id1089989.tgz

# mitochondrial RNA
# Closest relative available (Pclavus varians) is available in this GitHub repo as 'Cvarians_mitoRNA'
# Originally from Plese et al. (2021) doi: 10.1016/j.ympev.2020.107011
cp ~/bin/Cvarians_mitoRNA.fasta .

# Running a perl script that blasts the transcriptome against inputted contamination sequences to generate a final, clean transcriptome
echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=noshorts_Cvarians_Trinity.fasta contam=rRNA,arb-silva.de_2021-12-01_id1089989_tax_silva.fasta contam=Mt,Cvarians_mitoRNA.fasta table=Cvarians_contamination.txt passed=Cvarians.fasta" > contam
launcher_creator.py -j contam -n contam -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch contam.slurm

Cvarians.fasta
-------------------------
346223 sequences in input file
300 sequences look like contaminants
	rRNA	223
	Mt	77
345923 sequences passed all tests
-------------------------

echo "seq_stats.pl Cvarians.fasta > seqstats_Cvarians.txt" > seq_stats3
launcher_creator.py -j seq_stats3 -n seq_stats3 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats3.slurm


#------------------------------
## Final assembly!

Cvarians.fasta
-------------------------
345923 sequences.
1453 average length.
26720 maximum length.
500 minimum length.
N50 = 1874
502.7 Mb altogether (502749784 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

# Now proceed with separation of host/symbiont contigs in the script 'Cvarians_transcriptome_hostsym_separation_README'
