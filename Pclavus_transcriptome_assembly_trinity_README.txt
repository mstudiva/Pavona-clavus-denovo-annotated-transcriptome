## De novo transcriptome assembly pipeline for Pavona clavus, version December 31, 2024
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
echo "perl ~/bin/noshorts_v2.pl Pclavus_Trinity.fasta 300" > noshorts
echo "perl ~/bin/noshorts_v2.pl Pclavus_Galaxy.fasta 300" >> noshorts
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
# Running a perl script that blasts the transcriptome against inputted contamination sequences to generate a clean transcriptome
echo 'conda activate bioperl' > contam
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> contam
echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=Pclavus_Trinity_noshorts.fasta contam=rRNA,Pclavus_rRNA.fasta contam=Mt,Pclavus_mitoRNA.fasta table=Pclavus_Trinity_contamination.txt passed=Pclavus_Trinity_clean.fasta" >> contam
launcher_creator.py -j contam -n contam -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch --mem=200GB contam.slurm

# ExcludeFasta_v2.pl is failing due to a Bioperl error, so run it on its own
srun perl ~/bin/ExcludeFasta_v2.pl Pclavus_Trinity_contamination.txt Pclavus_Trinity_clean.fasta > Pclavus_Trinity_clean2.fasta

echo 'conda activate bioperl' > contam2
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> contam2
echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=Pclavus_Galaxy_noshorts.fasta contam=rRNA,Pclavus_rRNA.fasta contam=Mt,Pclavus_mitoRNA.fasta table=Pclavus_Galaxy_contamination.txt passed=Pclavus_Galaxy_clean.fasta" >> contam2
launcher_creator.py -j contam2 -n contam2 -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch --mem=200GB contam2.slurm

srun perl ~/bin/ExcludeFasta_v2.pl Pclavus_Galaxy_contamination.txt Pclavus_Galaxy_clean.fasta > Pclavus_Galaxy_clean2.fasta

Pclavus_Trinity_clean.fasta
-------------------------
124595 sequences in input file
1256 sequences look like contaminants
        rRNA    967
        Mt	289
123339 sequences passed all tests
-------------------------

Pclavus_Galaxy_clean.fasta
-------------------------
124645 sequences in input file
1233 sequences look like contaminants
        rRNA    965
        Mt	268
123412 sequences passed all tests
-------------------------

echo "seq_stats.pl Pclavus_Trinity_clean.fasta > seqstats_Pclavus_Trinity_clean.txt" > seq_stats3
echo "seq_stats.pl Pclavus_Galaxy_clean.fasta > seqstats_Pclavus_Galaxy_clean.txt" >> seq_stats3
launcher_creator.py -j seq_stats3 -n seq_stats3 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats3.slurm

Pclavus_Trinity_clean.fasta
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

Pclavus_Galaxy_clean.fasta
-------------------------
123412 sequences.
445 average length.
12396 maximum length.
300 minimum length.
N50 = 422
54.9 Mb altogether (54871003 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

# Now need to separate host/symbiont contigs into separate assemblies


#------------------------------
## Identify the most likely origin of each sequence by comparison to a protein DB from a single close relative and one or more databases of likely contaminants

# Find the most closely related genome through Uniprot: https://www.uniprot.org/proteomes
# For Pavona clavus, the closest relative is Pocillopora meandrina: https://www.uniprot.org/taxonomy/46732
# Closest zoox relative is Symbiodinium microadriaticum: https://www.uniprot.org/taxonomy/2951
# Download the proteome assemblies through 'Browse all ##,### entries' then download as '(FASTA) canonical & isoform'
# Finally, scp to your working directory

gunzip *.gz
mv uniprotkb_taxonomy_id_46732_2024_12_30.fasta Pmeandrina.fasta
mv uniprotkb_taxonomy_id_2951_2024_12_30.fasta Smicroadriaticum.fasta

cat Pmeandrina.fasta| grep '>' | wc -l
# 31995 -matches number in Uniprot
cat Smicroadriaticum.fasta| grep '>' | wc -l
# 43302 -matches number in Uniprot

module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls
# Making a blast database for each reference
echo "makeblastdb -in Pmeandrina.fasta -dbtype prot" >mdb
echo "makeblastdb -in Smicroadriaticum.fasta -dbtype prot" >>mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch mdb.slurm

# Blasting against reference proteomes
conda activate bioperl
echo 'conda activate bioperl' > origin
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> origin
echo "perl ~/bin/CompareContamSeq_blast+.pl -q Pclavus_Trinity_clean.fasta -s 45 -t Pmeandrina.fasta -c Smicroadriaticum.fasta" >> origin
launcher_creator.py -j origin -n origin -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch --mem=200GB origin.slurm

echo 'conda activate bioperl' > origin2
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> origin2
echo "perl ~/bin/CompareContamSeq_blast+.pl -q Pclavus_Galaxy_clean.fasta -s 45 -t Pmeandrina.fasta -c Smicroadriaticum.fasta" >> origin2
launcher_creator.py -j origin2 -n origin2 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch --mem=200GB origin2.slurm

Pclavus_Trinity_clean.fasta
-------------------------
123339 sequences input.
33735 of these matched Pmeandrina.fasta more closely than any contaminants.
11462 matched contaminants more closely than Pmeandrina.fasta.
78142 matched none of the supplied DB (nomatch.screened.fasta).
-------------------------

Pclavus_Galaxy_clean.fasta
-------------------------

-------------------------


#------------------------------
## Attempting to identify additional host/symbiont sequences in the no match assembly based on taxonomic ID of each sequence's best match in NCBI's nucleotide (nt) database

# NOTE: nt and taxdb databases downloaded 30 December 2024
mkdir ~/annotate/ncbi/nt
cd ~/annotate/ncbi/nt

echo 'conda activate blast_env' > get_nt
echo 'update_blastdb.pl --decompress nt --passive' >> get_nt
launcher_creator.py -j get_nt -n get_nt -q mediumq7 -t 12:00:00 -e studivanms@gmail.com
sbatch get_nt.slurm

srun update_blastdb.pl taxdb
tar -xzf taxdb.tar.gz

# Now need to update your blastdb path to include the location of taxonomy database files
nano ~/.bashrc
# Add the following text:
export BLASTDB="$HOME/annotate/ncbi/nt:$BLASTDB"
# Save, then source to make the path active
source ~/.bashrc

# Split the no match assembly into 80 chunks to parallelize and decrease computing time per chunk (shooting for <1000 sequences per chunk)
splitFasta.pl nomatch.screened.fasta 80

conda create -n blast_env perl-bioperl blast openssl
conda activate blast_env
conda install bioconda::perl-bioperl
conda install conda-forge::mamba
conda install bioconda::perl-digest-md5
conda install bioconda::blast

conda activate blast_env
module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls
echo 'conda activate bioperl' > bl_nomatch
echo 'module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls' >> bl_nomatch
# Create list of commands for blasting each subset chunk
for i in subset*nomatch*.fasta; do
  printf "blastn -query \"%s\" -db ~/annotate/ncbi/nt/nt -evalue 0.0001 -num_threads 4 -max_target_seqs 5 -outfmt \"6 qseqid sseqid evalue pident stitle staxids sscinames scomnames sblastnames sskingdoms salltitles stitle\" -out \"%s.br\"\n" "$i" "$i";
done > bl_nomatch
launcher_creator.py -j bl_nomatch -n bl_nomatch -q shortq7 -t 6:00:00 -e studivanms@gmail.com -N 5
sbatch bl_nomatch.slurm

# check blast progress
cat subset*.br | wc -l
# found matches in blast database for 29914 out of 152449 sequences

# generate combined blast report
cat subset*nomatch*.br > allblast.br

# scp the allblast.br file to your computer and run the taxonomizr.R script


#------------------------------
## Once back from R script, scp output .txt files to your HPC directory

# Use grep to generate a fasta file of the host/symbiont blast matches in the nomatch assembly
grep -w -A 1 -f nomatch_symbiont.txt nomatch.screened.fasta --no-group-separator > nomatch_symbiont.fasta

# The length should match the taxonomy matches in the taxonomizr.R script
cat nomatch_symbiont.fasta| grep '>' | wc -l
# 75

grep -w -A 1 -f nomatch_host.txt nomatch.screened.fasta --no-group-separator > nomatch_host.fasta
cat nomatch_host.fasta| grep '>' | wc -l
# 507

# Combine the host/symbiont nomatch assemblies with the original target/contam assemblies
cat Smicroadriaticum_trunc.screened.fasta nomatch_symbiont.fasta > Gerakladium.fasta

cat target.screened.fasta nomatch_host.fasta > Pclavus.fasta

echo "seq_stats.pl Pclavus.fasta > seqstats_Pclavus.txt" > seq_stats
echo "seq_stats.pl Gerakladium.fasta > seqstats_Gerakladium.txt" >> seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Pclavus.fasta
-------------------------
119824 sequences.
1901 average length.
26720 maximum length.
60 minimum length.
N50 = 2655
227.8 Mb altogether (227777238 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

Gerakladium.fasta
-------------------------
28670 sequences.
1375 average length.
21103 maximum length.
400 minimum length.
N50 = 1672
39.4 Mb altogether (39418006 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## GC content with BBMap package

scp the transcriptomes to your computer and follow the BBMap installation instructions here: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/

sh bbmap/stats.sh in=Pclavus.fasta
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2706	0.2300	0.2288	0.2706	0.0000	0.0000	0.0000	0.4588	0.0569

sh bbmap/stats.sh in=Gerakladium.fasta
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2041	0.2978	0.2926	0.2056	0.0000	0.0000	0.0000	0.5904	0.0403


#------------------------------
## Transcriptome completeness with gVolante/BUSCO

# Upload the transcriptomes to gVolante here: https://gvolante.riken.jp/analysis.html
# Use a cutoff length of '1', sequence type of 'Coding/transcribed (nucleotide)', ortholog search pipeline of 'BUSCO v5', and ortholog set of 'Metazoa' for host and 'Eukaryota' for symbiont

Pclavus.fasta
-------------------------
Completeness Assessment Results:
	Total # of core genes queried:    954
	# of core genes detected
		Complete:    887 (92.98%)
		Complete + Partial:    904 (94.76%)
	# of missing core genes:    50 (5.24%)
	Average # of orthologs per core genes:    4.16
	% of detected core genes that have more than 1 ortholog:    83.54
	Scores in BUSCO format:    C:93.0%[S:15.3%,D:77.7%],F:1.8%,M:5.2%,n:954
-------------------------

Gerakladium.fasta
-------------------------
Completeness Assessment Results:
	Total # of core genes queried:    255
	# of core genes detected
		Complete:    47 (18.43%)
		Complete + Partial:    63 (24.71%)
	# of missing core genes:    192 (75.29%)
	Average # of orthologs per core genes:    1.17
	% of detected core genes that have more than 1 ortholog:    8.51
	Scores in BUSCO format:    C:18.5%[S:16.9%,D:1.6%],F:6.3%,M:75.2%,n:255
-------------------------


# Now follow the script 'Pclavus_transcriptome_annotation_README' for generating annotation files for differential expression analysis
