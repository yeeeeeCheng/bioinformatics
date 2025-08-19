cd /bioinformatics/RNA_Seq_Pipeline/

SECONDS=0

#-------------------------------------------------------------------------------#

#Step1. Fetch SRA Files from NCBI SRA Database Using SRA Toolkit.

#export PATH=$PWD/apps/sratoolkit.ver-os-structure/bin:$PATH

#which fastq-dump

#vdb-config -i

#fastq-dump --stdout -X 2 SRR390728

#prefetch -p SRR_No

#fasterq-dump -p SRR_No --split-files --outdir results/Raw_reads/

#-------------------------------------------------------------------------------#

#Step2. Quality control using FastQC.

#cd apps/FastQC

#chmod u+x fastqc

#./fastqc ../../results/Raw_reads/sample_1_1.fq.gz ../../results/Raw_reads/sample_1_2.fq.gz -t N

#cd ../../results/Raw_reads/

#multiqc .

#-------------------------------------------------------------------------------#

#Step3. Trim adapter sequence and discard unwanted fragments with Cutadapt.

#read1=results/Raw_reads/Gz/sample_1_1.fq.gz

#read2=results/Raw_reads/Gz/sample_1_2.fq.gz

#output1=results/Clean_reads/sample_1_1_clean.fq.gz

#output2=results/Clean_reads/sample_1_2_clean.fq.gz

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=30 --poly-a -m 50 --pair-filter=any -j N -o $output1 -p $output2 $read1 $read2

#cd apps/FastQC

#chmod u+x fastqc

#./fastqc ../../results/Clean_reads/sample_1_1_clean.fq.gz ../../results/Clean_reads/sample_1_2_clean.fq.gz -t N

#cd ../../results/Clean_reads/

#multiqc .

#-------------------------------------------------------------------------------#

#Step4. Use HISAT2 to align the reads to the reference genome.

#export PATH=$PWD/apps/hisat2-2.2.1:$PATH

#gffread results/Indices/Genome_1/Genome_1.gff3 -T -o results/Indices/Genome_1/Genome_1.gtf

#hisat2_extract_splice_sites.py results/Indices/Genome_1/Genome_1.gtf > results/Indices/Genome_1/Genome_1.ss

#hisat2_extract_exons.py results/Indices/Genome_1/Genome_1.gtf > results/Indices/Genome_1/Genome_1.exon

#hisat2-build -p 8 --ss results/Indices/Genome_1/Genome_1.ss --exon results/Indices/Genome_1/Genome_1.exon results/Indices/Genome_1/Genome_1.fa results/Indices/Genome_1/Genome_1_idx

#sudo ln -s /Library/Frameworks/Python.framework/Versions/3.x.y/bin/python3 /usr/local/bin/python

#hisat2 -q --rna-strandness FR -p N -x results/Indices/Genome_1/Genome_1_idx -1 results/Clean_reads/sample_1_1_clean.fq -2 results/Clean_reads/sample_1_2_clean.fq | samtools sort -o results/Aligned/sample_1_FR.bam

#-------------------------------------------------------------------------------#

#Step5. Quantify the reads with featureCounts.

#export PATH=$PWD/apps/subread-2.1.1-Linux-x86_64/bin:$PATH

#cd results/Aligned/

#featureCounts -T N -p -a ../Indices/Genome_1/Genome_1.gtf -o ../Quantification/Genome_1/Genome_1.txt sample_1_FR.bam

#cd ../Quantification/Genome_1/

#multiqc .

#cat Genome_1.txt | cut -f1,7 | tr '\t' ',' > Genome_1_count_matrix.csv

#cat Sample_metadata.txt | tr '\t' ',' > Sample_metadata.csv

#-------------------------------------------------------------------------------#

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."