#!/bin/bash

###################################################
# miRNA analysis pipeline                         #
# May 2014                                        #
# By Adrià Closa                                  #
#                                                 #
#                                                 #
# Adapted by Ainhoa Garcia in 2020/21             #
#                                                 #
#                                                 #
# Dependencies (version used):                    #
# - miRBase (22.1)                                #                           
# - Python (2.7.15)                  	          #
# - FASTX-toolkit (0.0.13)          	          #
# - bowtie (1.2.2)                 	              #
# - SOLiD preprocess filter v2   	  	          #
# - cutadapt (1.0)                                #
# - samtools (1.8)                                #
###################################################


# -------------------------------------------------------------------
# STEP 0: Paths and useful pointers as environment variables
# -------------------------------------------------------------------
export MY_PATH="/shared/projects/COLONOMICS/mirna/userAG"
export AC_PATH="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics"
module load apps/python-2.7.15                                                  
export EXTRACT_HSA="${AC_PATH}/scripts/extract_hsa/extract_hsa.py"
export BOWTIE_PATH="/share/apps/bowtie/1.2.2"
  

# ------------------------------------------------------
# STEP 1: Downloading and processing reference sequences
# ------------------------------------------------------
# We download and uncompress the mature.fa.gz file from miRBase release 22.1
# at ftp://mirbse.org/pub/mirbase/22.1/mature.fa.gz, we rename the file mature.fa to mature_v22_1.fa

# STEP 1.1: Using a custom Python script, we keep only human miRNA sequences
# --------------------------------------------------------------------------
python $EXTRACT_HSA $MY_PATH/data/ref/mature_v22_1.fa $MY_PATH/data/ref/mature_hsa_v22_1.fa
# 2042 miRNA sequences left in release 20
# 2656 miRNA sequences left in release 22.1 (including 918 3p and 910 5p)


# STEP 1.2: Converting Uracils (U) to Timines (T) in reference sequence file
# --------------------------------------------------------------------------
/shared/projects/colab/LValle__miRNA_sang/scripts/fasta_nucleotide_changer -i $MY_PATH/data/ref/mature_hsa_v22_1.fa -o $MY_PATH/data/ref/mature_hsa_dna_v22_1.fa -d


# STEP 1.3: Creating bowtie index file for the miRNA mature reference in CS
# -------------------------------------------------------------------------
# Our index will be called mature_hsa_v22_1_index
$BOWTIE_PATH/bowtie-build -C $MY_PATH/data/ref/mature_hsa_dna_v22_1.fa $MY_PATH/data/ref/mature_hsa_v22_1


# STEP 1.4: Indexing genome for IGV and samtools mpileup
# ------------------------------------------------------
/share/apps/samtools/1.8/bin/samtools faidx $MY_PATH/data/ref/mature_hsa_dna_v22_1.fa



# ---------------------------------------
# STEP 2: Reads filtering and processing
# ---------------------------------------

# SOLID PREPROCESS FILTER 1
# *** This step was not executed since there is no need to repeat, we will keep Adrià's QVassess files ***

# Quality assessment for all the samples to assess
# the number of reads that do not pass different
# quality thresholds. This assessment is done with SOLiD preprocessing
# filter (https://github.com/fls-bioinformatics-core/HTS_SOLiD_Preprocessing)

# Output for each sample: matrices with rows for positions and columns for quality scores, each value is the count for that score in that position.
# Results in /shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/solid_preprocess_filter/QV_assess/SAMPLE.QVassess.txt

# START R CODE ---------------------------------------------------------------
# Get quality files names (254 quality files)
qual.files.short <- gsub(".qual", "", list.files("/shared/projects/COLONOMICS/data/mirna/ICO/Fasta_Files/CLX", pattern=".qual", full.names=FALSE));
qual.files.full <- list.files("/shared/projects/COLONOMICS/data/mirna/ICO/Fasta_Files/CLX", pattern=".qual", full.names=TRUE);

#Output folder
output <- paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/solid_preprocess_filter/QV_assess/", qual.files.short, ".QVassess.txt", sep="");

#Creating the command to call the perl script
commands <- paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/scripts/solid_preprocess_filter/QV_run_assessment.pl -i ", qual.files.full," -m y -o ", output, sep="");

#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_run.sh");
for (i in (1:length(commands)))
{
  cat("#!/bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_",i,".sh",sep=""),append=TRUE);
  cat("#$ -N preproc_",i,"\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_",i,".sh",sep=""),sep="",append=TRUE);
  cat(commands[i],"\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_",i,".sh",sep=""), append=TRUE);
  cat("qsub",paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_",i,".sh",sep=""),"&\n",
      file="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/QV_assess/QV_assess_run.sh", append=TRUE);
}
# END R CODE ---------------------------------------------------------------
./QV_assess_run.sh


# SOLID PREPROCESS FILTER 2
# *** This step was not executed since there is no need to repeat, we will keep Adrià's QVassess files ***

# Load QV assessment results and remove reads with miscalls
# also with SOLiD preprocessing filter (https://github.com/fls-bioinformatics-core/HTS_SOLiD_Preprocessing)

# Obtain the total reads, sampled reads, reads once miscalls removed, >=Q20 reads, and different pXeX filters for each sample and visualize results.

# START R CODE ---------------------------------------------------------------
#Loading QV files results for all samples
QV.assess.files.full <- list.files("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/solid_preprocess_filter/QV_assess", pattern=".txt", full.names=TRUE);
samples <- gsub(".Q", "", substr(list.files("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/solid_preprocess_filter/QV_assess", pattern=".txt"), 1, 9));
total.reads <- NULL;
sampled.reads <- NULL;
miscalls <- NULL;
QV.20 <- NULL;

# p is for minimum polyclonal counts: entire read poor quality or missequenced, two templates amplified in a single bead resulting in a hybrid sequence without match in the ref. genome
# e is for maximum errors infentifyed: single color call errors, QV of 10 or less had higher probability of being erroneous
p1.eoff <- NULL;
p3.eoff <- NULL;
p5.eoff <- NULL;
p1.e5 <- NULL;
p1.e3 <- NULL;
p5.e0 <- NULL;

# Numbers 2, 65, 124, 183, 242... are hardcoded, if you rerun this script 
# please check in the future if they are still correct!!
for (i in (1:length(QV.assess.files.full)))
{
  tmp <- readLines(QV.assess.files.full[i]);
  
  aux <- unlist(strsplit(tmp[1]," - "))[2];
  total.reads <- c(total.reads, as.numeric(aux));
  
  aux <- unlist(strsplit(tmp[2]," - "))[2];
  sampled.reads <- c(sampled.reads, as.numeric(aux));
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[65]," "))));
  miscalls <- c(miscalls, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[124]," "))));
  QV.20 <- c(QV.20, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[183]," "))));
  p1.eoff <- c(p1.eoff, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[242]," "))));
  p3.eoff <- c(p3.eoff, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[301]," "))));
  p5.eoff <- c(p5.eoff, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[360]," "))));
  p1.e5 <- c(p1.e5, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[419]," "))));
  p1.e3 <- c(p1.e3, aux);
  
  aux <- sum(as.numeric(unlist(strsplit(tmp[478]," "))));
  p5.e0 <- c(p5.e0, aux);
}

names(total.reads) <- samples;
names(sampled.reads) <- samples;
names(miscalls) <- samples;
names(QV.20) <- samples;
names(p1.eoff) <- samples;
names(p3.eoff) <- samples;
names(p5.eoff) <- samples;
names(p1.e5) <- samples;
names(p1.e3) <- samples;
names(p5.e0) <- samples;

miscalls <- miscalls/sampled.reads;
QV.20 <- QV.20/sampled.reads;
p1.eoff <- p1.eoff/sampled.reads;
p3.eoff <- p3.eoff/sampled.reads;
p5.eoff <- p5.eoff/sampled.reads;
p1.e5 <- p1.e5/sampled.reads;
p1.e3 <- p1.e3/sampled.reads;
p5.e0 <- p5.e0/sampled.reads;

# Determining which samples fall outside the big majority for the different parameters
# Number of total reads
# WE DO NOT APPLY THIS FILTER. ALL SAMPLES HAVE AT LEAST ALMOST 2M READS.
totalreads.fail <- names(total.reads[total.reads < 2000000]); # "A2096_N" "D2125_M"

# Percentage of missings without calls
misc.fail <- names(miscalls[miscalls < 0.94]); # "E2023_N" "E2023_T" "E2092_N" "E2092_T" "G2120_M"

# Percentage of reads with an average QV > 20
QV20.fail <- names(QV.20[QV.20 < 0.6]); # "E2023_T"   "Q2063_T"   "T2070_T_1"

# p=1 => minimum 1 call with QV>=25 in ONLY the first 10 bases of the read
# e=5 => maximum 5 calls with QV<=10 in ALL the read

# Percentage of reads kept after applying filter p=1 + e=off
p1eoff.fail <- names(p1.eoff[p1.eoff < 0.93]); # "E2023_N" "E2023_T" "E2092_N" "E2092_T" "G2120_M" "Q2040_T" "X2011_T" "Z2061_T"

# Percentage of reads kept after applying filter p=3 + e=off
p3eoff.fail <- names(p3.eoff[p3.eoff < 0.84]); # "X2011_T"

# Percentage of reads kept after applying filter p=5 + e=off
p5eoff.fail <- names(p5.eoff[p5.eoff < 0.78]); #"X2011_T"

# Percentage of reads kept after applying filter p=1 + e=5
p1e5.fail <- names(p1.e5[p1.e5 < 0.4]); # "B2081_T" "K2068_N"

# Percentage of reads kept after applying filter p=1 + e=3
p1e3.fail <- names(p1.e3[p1.e3 < 0.3]); #"B2081_T" "C2090_T" "D2056_N" "G2005_N" "K2068_N" "M2029_N" "M2052_T" "X2080_N"

# Percentage of reads kept after applying filter p=5 + e=0
p5e0.fail <- names(p5.e0[p5.e0 < 0.1]); #"0""B2081_T" "C2090_T" "D2056_N" "G2005_N" "K2068_N" "L2089_N" "M2029_N""M2052_T" "P2032_T" "S2062_T" "X2080_N"


# Removal of reads containing miscalls and all reads containing negative quality scores.
# The output will result in 4 files: 2 csfasta and 2 QV files for each sample: T for those who passed the filter and U for those that were not (misscalls).

#Getting cs fasta and quality files names
sample.names <- gsub(".csfasta", "", list.files("/shared/projects/COLONOMICS/data/mirna/ICO/Fasta_Files/CLX", pattern=".csfasta", full.names=FALSE));

#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_run.sh");
for (i in (1:length(sample.names)))
{
  cat("#!/bin/bash\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_",i,".sh",sep=""),
      append=TRUE);
  cat("#$ -N remove_",i,"\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_",i,".sh",sep=""),sep="",
      append=TRUE);
  cat(paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/scripts/solid_preprocess_filter/SOLiD_preprocess_filter_v2.pl ",
            "-f /shared/projects/COLONOMICS/data/mirna/ICO/Fasta_Files/CLX/",
            sample.names[i],
            ".csfasta ",
            "-g /shared/projects/COLONOMICS/data/mirna/ICO/Fasta_Files/CLX/",
            sample.names[i],
            ".qual ",
            "-o /shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/solid_preprocess_filter/remove_miscalls/",
            sample.names[i],
            "_nomis ",
            "-v -x no -p no -y no -n yes\n",
            sep=""),
      file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_",i,".sh",sep=""),
      append=TRUE);
  
  cat("qsub", paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_",i,".sh",sep=""), "&\n",
      file="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/solid_preprocess_filter/remove_miscalls/remove_miscalls_run.sh", append=TRUE);
}

# END R CODE ---------------------------------------------------------------

./remove_miscalls_run.sh

# ------------------------
# STEP 3: Adaptor trimming
# ------------------------

# Options chosen for cutadapt are:
# -- allow 1 mismatch on the adapter sequence (=> -e 0.1 is OK)
# -- min and max lengths = 16, 28
# -- do not keep sequences with NO MATCH with the adapter

# START R CODE ---------------------------------------------------------------
#Getting quality files names
path <- "/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/solid_preprocess_filter/remove_miscalls";
output.path <- "/shared/projects/COLONOMICS/mirna/userAG/cutadapt";

csfasta.files.short <- gsub(".csfasta", "", list.files(path, pattern="_T_F3.csfasta", full.names=FALSE));
csfasta.files.full <- list.files(path, pattern="_T_F3.csfasta", full.names=TRUE);

# Test cutadapt with same parameters than Adrià: -e (maximum error rate), -z (compression level), -m and -M (minimum and maixmum length).
#Creating the command to call the perl script
commands <- paste("/share/apps/2_old/cutadapt-1.0/cutadapt -c -e 0.1 -z -m 16 -M 28 ",
                  "--untrimmed-output=",output.path,"/untrimmed/",csfasta.files.short,"_untrimmed.csfastq ",
                  "-a 330201030313112312 ",
                  path,"/",csfasta.files.short,".csfasta ",
                  path,"/",csfasta.files.short,"_QV.qual",
                  " > ",output.path,"/trimmed/",csfasta.files.short,"_trimmed.csfastq ",
                  sep="");

#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_run.sh");
for (i in (1:length(commands)))
{
  cat("#!/bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_",i,".sh",sep=""),append=TRUE);
  cat("#$ -N cut_",i,"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_",i,".sh",sep=""),sep="",append=TRUE);	
  cat("#$ -e ",output.path,"/summary/",csfasta.files.short[i],"_summary.txt\n",
      file=paste("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_",i,".sh",sep=""),sep="",append=TRUE);
  cat(commands[i],"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_",i,".sh",sep=""), append=TRUE);
  cat("qsub", paste("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_",i,".sh",sep=""), "&\n",
      file="/shared/projects/COLONOMICS/mirna/userAG/cutadapt/scripts/cutadapt_run.sh", append=TRUE);
}
# END R CODE ---------------------------------------------------------------

./cutadapt_run.sh

# Summary of cutadapt parameters (as returned by cutadapt after its execution - it is the same for all files)
# Type of Input - F
# Polyclonal Analysis - n
# Min count for Polyclonal Analysis - 0
# Min QV for Polyclonal Analysis - 25
# Error Analysis - n
# Max count permitted errors - 3
# Max QV to consider an error - 10
# Removal of reads with negative QV score - y
# Truncation is off - n
# Quality Value analysis is off - n
# Output QV files - y

# START R CODE ---------------------------------------------------------------
# Now we process the output of cutadapt for all samples
path <- "/shared/projects/COLONOMICS/mirna/userAG/cutadapt/";
csfasta.files.short <- gsub(".csfasta", "", list.files(paste(path,"/trimmed",sep=""), pattern="_T_F3_trimmed.csfastq", full.names=FALSE));
sum.files <- gsub("_trimmed.csfastq", "", paste(path,"summary/",csfasta.files.short,"_summary.txt",sep=""));
cut.sum <- vector("list",length(sum.files));
i <- 1;
for (f in sum.files)
{
  cut.sum[[i]] <- readLines(f);
  cat(i,"\n");
  i <- i + 1;
}

samples <- gsub("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/summary/","",sum.files);
samples <- gsub("_nomis_T_F3_summary.txt","",samples);

trimmed.reads <- as.numeric(unlist(lapply(strsplit(unlist(lapply(cut.sum,FUN=function(x){x[[5]]}))," "),FUN=function(x){gsub("%)","",x[length(x)])})));
jpeg("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/graphs/trimmed_reads.jpg",width=1024,height=1024,quality=100);
hist(trimmed.reads,main="Trimmed reads",xlab="Percentage of trimmed reads",prob=TRUE,ylim=c(0,0.055));
lines(density(trimmed.reads),col="blue",lwd=3);
dev.off();

aux <- data.frame(samples,trimmed.reads);
aux[aux[,2]<50,];

adapter.lengths <- numeric(16);
for (i in (1:16))
{
  adapter.lengths[i] <- sum(as.numeric(gsub(paste(i+2,"\t",sep=""),"",unlist(lapply(cut.sum,FUN=function(x){x[[i+16]]})))));
}

jpeg("/shared/projects/COLONOMICS/mirna/userAG/cutadapt/graphs/adapter_lengths.jpg",width=1024,height=1024,quality=100);
barplot(adapter.lengths,main="Adapter length",xlab="Length (nucleotides)",names.arg=3:18,ylab="Frequency");
dev.off();

# Determining what samples fall outside the big majority regarding the proportion of trimmed reads
trimmed.fail <- samples[trimmed.reads < 50]; #"B2081_T" "C2090_T" "D2056_N" "G2005_N" "K2068_N" "M2029_N" "M2052_T" "S2062_T" "T2070_T_1" "X2080_N"



# -----------------------------------------------------
# STEP 4A: Alignment against mature miRNAs (1 mismatch)
# -----------------------------------------------------
# Parameters for bowtie: -t (time execution), --sam (sam output), -C (colorspace, for SOLiD reads), -q (input fastq)
# -n (number of max missmatches allowed), --all (report all valid alignments per read or pair), --un (output unaligned reads)

# Bowtie version with colorspace (v.1.2.2). Bowtie2 does not align colorspace reads.


#Getting csfastq files names
path <- "/shared/projects/COLONOMICS/mirna/userAG/";

csfasta.files.short <- gsub(".csfastq", "", list.files(paste(path,"cutadapt/trimmed",sep=""), pattern=".csfastq", full.names=FALSE));
csfasta.files.full <- list.files(paste(path,"cutadapt/trimmed",sep=""), pattern=".csfastq", full.names=TRUE);

commands <- paste("/share/apps/bowtie/1.2.2/bowtie -t --sam -C -q -n 1 --all ",
                  "--un ",path,"bowtie/n_1/unaligned/",csfasta.files.short,"_unaligned.sam ",
                  path,"data/ref/mature_hsa_v22_1 ",
                  path,"data/raw/cutadapt/trimmed/",csfasta.files.short,".csfastq ",
                  path,"bowtie/n_1/aligned/",csfasta.files.short,".sam",
                  sep="");

#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_run.sh");
for (i in (1:length(commands)))
{
  # We need to 
  cat("#!/bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_",i,".sh",sep=""),append=TRUE);
  cat("#$ -N bowtie_",i,"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_",i,".sh",sep=""),sep="",append=TRUE);	
  cat("#$ -e ",path,"bowtie/n_1/summary/",csfasta.files.short[i],"_bowtie_summary.txt\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_",i,".sh",sep=""),sep="",append=TRUE);
  cat(commands[i],"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_",i,".sh",sep=""), append=TRUE);
  cat("qsub", paste("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_",i,".sh",sep=""), "&\n", file="/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/scripts/bowtie_run.sh", append=TRUE);
}

# END R CODE ---------------------------------------------------------------
./bowtie_run.sh

# START R CODE ---------------------------------------------------------------
path <- "/shared/projects/COLONOMICS/mirna/userAG/";
prop.mapped.reads <- numeric(0);
aligned.reads <- numeric(0);
total.reads <- numeric(0);
f <- list.files(paste(path,"bowtie/n_1/summary",sep=""),full.names=TRUE);
samples <- gsub("_nomis_T_F3_trimmed_bowtie_summary.txt","",f);
samples <- gsub("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_1/summary/","",samples);
for (i in f)
{
  aux <- readLines(i)[6];
  aux2 <- unlist(strsplit(unlist(strsplit(aux,": "))[2]," "));
  
  aux4 <- as.numeric(aux2[1]);
  
  aux2 <- as.numeric(substr(aux2[2],2,6));
  aux2 <- as.numeric(aux2);
  
  aux3 <- readLines(i)[5];
  aux3 <- as.numeric(unlist(strsplit(aux3,": "))[2]);
  
  prop.mapped.reads <- c(prop.mapped.reads, aux2);
  total.reads <- c(total.reads, aux3);
  aligned.reads <- c(aligned.reads, aux4);
}

names(prop.mapped.reads) <- samples;
names(total.reads) <- samples;
names(aligned.reads) <- samples;


# Filter by the absolute number of aligned reads
aligned.fail <- names(aligned.reads[aligned.reads < 500000]); # "B2081_T" "C2067_T_1" "C2090_T" "D2056_N" "K2068_N" "M2029_N" "M2052_T" "N2036_T" "R2025_T" "X2080_N" "Z2061_T" 


# -----------------------------------------------------
# STEP 4B: Alignment against mature miRNAs (2 mismatches)
# -----------------------------------------------------
# Same than in step 4A but allowing 2 mismatches

#Getting csfastq files names
path <- "/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/";

csfasta.files.short <- gsub(".csfastq", "", list.files(paste(path,"cutadapt/trimmed",sep=""), pattern=".csfastq", full.names=FALSE));
csfasta.files.full <- list.files(paste(path,"cutadapt/trimmed",sep=""), pattern=".csfastq", full.names=TRUE);

#Creating the command to call the perl script
commands <- paste("/share/apps/bowtie/1.2.2/bowtie -t --sam -C -q -n 2 --all ",
                  "--un ",path,"bowtie/n_1/unaligned/",csfasta.files.short,"_unaligned.sam ",
                  path,"data/ref/mature_hsa_v22_1 ",
                  path,"data/raw/cutadapt/trimmed/",csfasta.files.short,".csfastq ",
                  path,"bowtie/n_1/aligned/",csfasta.files.short,".sam",
                  sep="");


#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_run.sh");
for (i in (1:length(commands)))
{
  # We need to 
  cat("#!/bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_",i,".sh",sep=""),append=TRUE);
  cat("#$ -N bowtie_",i,"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_",i,".sh",sep=""),sep="",append=TRUE);	
  cat("#$ -e ",path,"bowtie/n_2/summary/",csfasta.files.short[i],"_bowtie_summary.txt\n",file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_",i,".sh",sep=""),sep="",append=TRUE);
  cat(commands[i],"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_",i,".sh",sep=""), append=TRUE);
  cat("qsub", paste("/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_",i,".sh",sep=""), "&\n", file="/shared/projects/COLONOMICS/mirna/userAC/hairpin_colonomics/sh_hairpin/bowtie/n_2/bowtie_run.sh", append=TRUE);
}
# END R CODE ---------------------------------------------------------------

./bowtie_run.sh

# START R CODE -------------------------------------------------------------
path <- "/shared/projects/COLONOMICS/mirna/userAG/";
prop.mapped.reads <- numeric(0);
aligned.reads <- numeric(0);
total.reads <- numeric(0);
f <- list.files(paste(path,"bowtie/n_2/summary",sep=""),full.names=TRUE);
samples <- gsub("_nomis_T_F3_trimmed_bowtie_summary.txt","",f);
samples <- gsub("/shared/projects/COLONOMICS/mirna/userAG/bowtie/n_2/summary/","",samples);
for (i in f)
{
  aux <- readLines(i)[6];
  aux2 <- unlist(strsplit(unlist(strsplit(aux,": "))[2]," "));
  
  aux4 <- as.numeric(aux2[1]);
  
  aux2 <- as.numeric(substr(aux2[2],2,6));
  aux2 <- as.numeric(aux2);
  
  aux3 <- readLines(i)[5];
  aux3 <- as.numeric(unlist(strsplit(aux3,": "))[2]);
  
  prop.mapped.reads <- c(prop.mapped.reads, aux2);
  total.reads <- c(total.reads, aux3);
  aligned.reads <- c(aligned.reads, aux4);
}

names(prop.mapped.reads) <- samples;
names(total.reads) <- samples;
names(aligned.reads) <- samples;

# Filter by the absolute number of aligned reads
aligned.fail <- names(aligned.reads[aligned.reads < 500000]); # "B2081_T" "C2090_T" "D2056_N" "K2068_N" "M2029_N" "M2052_T" "X2080_N"

# ------------------------------------------------------------------------------------
# STEP 5: #SAM2BAM conversion AND filtering of unmapped reads and reverse strand reads
# ------------------------------------------------------------------------------------
# samtools parameters: remove unmapped reads (F4) and reverse strand reads (F16)

# START R CODE ---------------------------------------------------------------
#Getting quality files names
path <- "/shared/projects/COLONOMICS/mirna/userAG/";

sam.files.short <- gsub(".sam", "", list.files(paste(path,"bowtie/n_2/aligned",sep=""), pattern=".sam", full.names=FALSE));
sam.files.full <- list.files(paste(path,"bowtie/n_2/aligned",sep=""), pattern=".sam", full.names=TRUE);

#Creating the command to call the perl script
commands <- paste("/share/apps/samtools/1.8/bin/samtools view -uS ",
                  sam.files.full,"  ",
                  "/share/apps/samtools/1.8/bin/samtools view -u -F4 -F16 - | /share/apps/samtools/1.8/bin/samtools sort -o ",path,"samtools/sam2bam/n_2/bam_files/",sam.files.short,".bam\n",
                  "/share/apps/samtools/1.8/bin/samtools index /shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/bam_files/",sam.files.short,".bam",
                  sep="");

#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/samtools_run.sh");
for (i in (1:length(commands)))
{
  cat("#!/bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""),append=TRUE);
  cat("#$ -N sam_",i,"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""),sep="",append=TRUE);	
  cat("#$ -j y\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""),sep="",append=TRUE);	
  cat(paste("#$ -o /shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/logs/samtools_",i,"\n",sep=""),file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""),sep="",append=TRUE);
  cat(commands[i],"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""), append=TRUE);
  cat("qsub", paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/scripts/samtools_",i,".sh",sep=""), "&\n", file="/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/samtools_run.sh", append=TRUE);
}

# END R CODE ---------------------------------------------------------------

./samtools_run.sh

# -----------------------------------------------------------------------
# STEP 6: Computing coverage with samtools mpileup and processing results
# -----------------------------------------------------------------------

# START R CODE ---------------------------------------------------------------
path <- "/shared/projects/COLONOMICS/mirna/userAG/samtools/sam2bam/n_2/bam_files/";

bam.files.short <- gsub(".bam", "", list.files(path, pattern=".bam", full.names=FALSE));
bam.files.full <- list.files(path, pattern=".bam", full.names=TRUE);

bam.files.short <- bam.files.short[-grep(".bai",bam.files.short)];
bam.files.full <- bam.files.full[-grep(".bai",bam.files.full)];

#Creating the command to call the perl script
# At a position, read maximally 1000000 reads per input file.
commands <- paste("/share/apps/samtools/1.8/bin/samtools mpileup -d 1000000 -f /shared/projects/COLONOMICS/mirna/userAG/data/ref/mature_hsa_dna_v22_1.fa ",
                  bam.files.full," | cut -f1-4",
                  sep="");

#Creating the .sh files
cat("#!/bin/bash\n", file="/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_run.sh");
for (i in (1:length(commands)))
{
  cat("#!/bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""));
  cat("#$ -S /bin/bash\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""),append=TRUE);
  cat("#$ -N mpil_",i,"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""),sep="",append=TRUE);
  cat("#$ -o /shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/counts/",bam.files.short[i],"_counts.txt\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""),sep="",append=TRUE);
  cat("#$ -cwd\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""),append=TRUE);
  cat(commands[i],"\n",file=paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""), append=TRUE);
  
  cat("qsub", paste("/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_",i,".sh",sep=""), "&\n", file="/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/scripts/samtools_mpileup_run.sh", append=TRUE);
}

# END R CODE -----------------------------------------------------------------
./samtools_mpileup_run.sh

# START R CODE ---------------------------------------------------------------

path <- "/shared/projects/COLONOMICS/mirna/userAG/samtools/mpileup/n_2/counts";

count.files.full <- list.files(path, pattern="_counts.txt", full.names=TRUE);
count.files.short <- gsub("_nomis_T_F3_trimmed", "", gsub("_counts.txt", "", list.files(path, pattern="_counts.txt", full.names=FALSE)));

counts <- list();
for (i in (1:length(count.files.full)))
{
  aux <- read.table(count.files.full[i], sep="\t", header=FALSE, as.is=TRUE, strip.white=TRUE);
  counts[[i]] <- aux;
  names(counts)[i] <- count.files.short[i];
}

counts.med <- lapply(counts,FUN=function(x){tapply(x$V4,x$V1,FUN=median)});

mirnas <- sort(unique(unlist(lapply(counts.med,FUN=function(x){names(x)}))));

counts.raw <- matrix(0,nrow=length(mirnas),ncol=length(count.files.short));
rownames(counts.raw) <- mirnas;
colnames(counts.raw) <- count.files.short;

i <- 1;
while(i <= length(counts.med))
{
  j <- 1;
  while (j <= length(counts.med[[i]]))
  {
    counts.raw[names(counts.med[[i]])[j], names(counts.med)[i]] <- counts.med[[i]][j];
    
    j <- j + 1;
  }
  
  i <- i + 1;
}
# Finally we got the raw counts matrix!
counts.raw <- trunc(counts.raw);
write.table(counts.raw,file="/shared/projects/COLONOMICS/mirna/userAG/raw_counts_matrix/n_2/mirna_raw_counts_n2.txt",quote=FALSE,sep="\t");

# END R CODE ---------------------------------------------------------------