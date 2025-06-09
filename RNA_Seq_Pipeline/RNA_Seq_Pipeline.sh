#####################################ATTENTION###################################
#                                                                               #
#     Please check and change the file name and path to the files that will     #
#     be processed later by this file.                                          #
#                                                                               #
#####################################ATTENTION###################################


#Change the working directory to the directory for RNA-Seq analysis pipeline.
#cd /Users/yeecheng/Bioinformatics/RNA_Seq_Pipeline/

#Set the timer of the terminal to 0 so that the timer can record the time length for running the whole script.
SECONDS=0

#Step1. Fetch SRA Files from NCBI SRA Database Using SRA Toolkit.

#First, we need to set up the PATH each time we are going to fetch SRA files from NCBI SRA Database using SRA Toolkit. This could be set up once and for all, yet it could involve the change of some base file in the computer, which should be carefully handled.
#export PATH=$PWD/apps/sratoolkit.3.1.1-mac-arm64/bin:$PATH

#Verify whether SRA Toolkit can be executed successfully. It should send the path the file fastq-dump was located back.
#which fastq-dump

#Configure SRA Toolkit.
#vdb-config -i

#Now, a blue background should appear with some texts on it. Press tab and enter to move between different options and confirm choices, respectively. Make sure the box in the front of "Enable Remote Access" on the main page contained a X. Go to the "Cache" tab, and make sure the box in the front of "local file-caching" contained a X. Set the "Location of user-repository" to an empty directory, which will store the sra files fetched by the toolkit. At this step, we may have to type the route manually. At last, go to exit option and press enter to save and leave the configuration mode.

#After that, we can type the script in the next row to test whether SRA Toolkit can function normally. If the outpus are the first two result of the sequencing result with quality indices, it means the SRA Toolkit can function normally.
#fastq-dump --stdout -X 2 SRR390728

#Once the verifying procedure were all passed, SRA Toolkit can be used multiple times to fetch the SRA files we need without further setting.

#Rrefetch the SRA file. We can download SRA files without prefetching step, but the downloading process would be very long and cannot be rescued once the download process was disrupted. SRA files contained the information SRA Toolkit need to construct the whole fastq file, only require less storgae than the fastq files. If the SRA files weren't deleted, the corresponding fastq files can be construct multiple times. Also, if the prefetching step was interrupted, we can still resume the process by typing the same command again. With -p, we can see the progress bar so that we can take a grasp of how much files have already been downloaded.
#prefetch -p SRR7959046

#Once the prefetch files have already been downloaded, just type the script in the following row to construct the fastq file we need. The description of --split-files can be found in the documentation of SRA Toolkit on its pages in GitHub. Also, with the argument --outdir, we can specify the path we want the fastq files to be located in.
#fasterq-dump -p SRR7959046 --split-files --outdir results/Raw_reads/

#-------------------------------------------------------------------------------#

#Step2. Quality control using FastQC

#FastQC can be performed in two ways, in GUI, or in CLI. Only the CLI can be incorporated into the pipeline within single run. Once the analysis was completed, FastQC will generate a HTML file report. The interpretation of the results can be found on the FastQC website.

#To run FastQC in a CLI, put the main directory of FastQC to the same folder where the main folder of SRA Toolkit was located. 

#Go into the directory of FastQC.
#cd apps/FastQC

#Type the script in the next row to authorize the user to run fastqc program.
#chmod u+x fastqc

#Next, we can perform FastQC on the fastq files downloaded and stored in the "results" directory. Type the script in the next row to run FastQC. Please change the SRR_accession1 and SRR_accession2 in the script to specify the SRA files. More files can be added in the script to perform batch analysis. The command "./fastqc" was used to wake up the program.
#./fastqc ../../results/Raw_reads/gz_files/CK21_L1_2.fq -t 8

#Since FastQC generate the report the each fastq file, we can use MultiQC to gather all the report for all the fastq files in the Rae_reads folder.
#cd ../../results/Raw_reads/

#multiqc .

#-------------------------------------------------------------------------------#

#cd ../../

#Step3. Trim adapter sequence and poly-A tail (poly-T tail on the reverse strand) and discard fragments with lengths shorter than 50 bp with Cutadapt.

#Since Illumina has announced the adapter sequence publicly, we can manually input the adapter sequence following the -a and -A argument, respectively. Also, we can assign some strings with the path for the files for the input and output to shorten the main usage command.

#--pair-filter=any was used to discard the fragments either one or both of the reads meet the criteria. Argument -j can be used to specify the threads we can allocate for the program to use. 

#read1=results/Raw_reads/gz_files/CK21_L1_1.fq
#read2=results/Raw_reads/gz_files/CK21_L1_2.fq

#output1=results/Clean_reads/CK21_L1_1_clean.fastq
#output2=results/Clean_reads/CK21_L1_2_clean.fastq

#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --poly-a -m 50 --pair-filter=any -j 8 -o $output1 -p $output2 $read1 $read2

#After the trimming process was complete, we can perform FastQC again to confirm the status of the trimmed files meet our standard (or become the clean reads).

#cd apps/FastQC

#chmod u+x fastqc

#./fastqc ../../results/Clean_reads/CK21_L1_1_clean.fastq ../../results/Clean_reads/CK21_L1_2_clean.fastq -t 8

#cd ../../results/Clean_reads/

#multiqc .

#cd ../../

#-------------------------------------------------------------------------------#

#Step4. Use HISAT2 to align the reads to the reference genome.

#export PATH=$PWD/apps/hisat2-2.2.1:$PATH

#To build the index for the species of interest, use the following commands to generate the splice sites info and exons info with the gtf file. Sometimes the database won't provide the gtf file, but provide the gff file. We can convert the gff file with gffread program in the Cufflinks.
#gffread results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.gff3 -T -o results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.gtf
#hisat2_extract_splice_sites.py results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.gtf > results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.ss
#hisat2_extract_exons.py results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.gtf > results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.exon
#hisat2-build -p 8 --ss results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.ss --exon results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_v4.exon results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_genome_v4.fa results/Indices/Cucumis_melo_DHL92_genome_v4/DHL92_genome_v4_idx

#Then, the following command is the main usage of HISAT2. The "-q" means our input files come in fastq format. The "-x" means the basename of the indices of the reference genome. Usually, they start with the name of the directory that stored the index files, followed by the title of the indices that each file should contain in the beginnning in their file name. samtools is employed in the pipeline to convert the sam files generated by HISAT2 to bam files simultaneously.

#If hisat2-build can't run and show "env: python no such file or directroy" message, it's probably because the default path of python in HISAT2 is different than that of python that was installed on this computer. Entering the follwoing command to make a symlink, which is kind of like a temporary link of the PATH, so that the computer can temporarily add the python to the path we assigned, which is /usr/local/bin/python in this case.

#sudo ln -s /Library/Frameworks/Python.framework/Versions/3.13/bin/python3 /usr/local/bin/python

#hisat2 -q --rna-strandness FR -p 8 -x results/Indices/Cucumis_melo_ASM554921v1/GCA_005549215.1_ASM554921v1_idx -1 results/Clean_reads/CK21_L1_1_clean.fastq -2 results/Clean_reads/CK21_L1_2_clean.fastq | samtools sort -o results/Aligned/CK21_FR.bam

#-------------------------------------------------------------------------------#

#Step5. Quantify the reads with featureCounts.

#export PATH=$PWD/apps/subread-2.0.8-macOS-arm64/bin:$PATH

#cd results/Aligned/

#The usage of featureCounts require the annotation file comes in gtf format. It's recommended to move to the directory containing the bam files to make the output count matrix looks more tidy.
#featureCounts -T 8 -p -a ../Indices/Cucumis_melo_ASM554921v1/GCA_005549215.1_ASM554921v1_genomic.gtf -o ../Quantification/Melon_ASM554921v1/Melon_ASM554921v1.txt CK21_FR.bam

#cd results/Quantification/Melon_ASM554921v1/

#Multiqc can be used to check the status of the quantification.
#multiqc .

#Use the following pipeline to extract the info we need to construct the count matrix.
#cat Melon_ASM554921v1.txt | cut -f1,7 | tr '\t' ',' > Melon_ASM554921v1.txt_count_matrix.csv
cat Melon_sample_data.txt | tr '\t' ',' > Melon_sample_data.csv

#-------------------------------------------------------------------------------#

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."