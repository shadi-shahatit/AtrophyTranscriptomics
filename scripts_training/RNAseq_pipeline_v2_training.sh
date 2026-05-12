################################# Transcriptomics Analysis Training - Al-Zghoul Lab, JUST
################################# Shadi Shahatit - RA, 2026

#### Table of Contents

#### 1. RNA-seq Quantification via Salmon 
####	1.0	Setup Environment and Directories
####	1.1	Run FastQC on Raw Reads
####	1.2	Run Trimmomatic
####	1.3	Run FastQC on Trimmed Reads
####	1.4	Reference Transcriptome Indexing 
####	1.5	Salmon Quantification
####	1.6	Appendix

#########################################################################################################################################################################################################
## RNA-seq Quantification via Salmon 
#########################################################################################################################################################################################################

#################################
## Step 0: Setup Environment and Directories
#################################

## Conda setup

# Create and name conda environment
conda create -y -n salmon
# Activate conda env 
source activate salmon
# install fastqc, trimmomatic, and salmon from the links in the slides

## Directories setup

mkdir -p raw_fastq fastqc_reports trimmed_fastq quant_r115

# RUN samples one by one instead of using a loop

## Raw data setup

cd raw_fastq

# Fetch the files for raw data fastq files from the links in fastq2down_liver_training.txt in https://github.com/shadi-shahatit/AtrophyTranscriptomics/tree/main/scripts_training

# sample names
# 1st run (1 sample; 2 files)
D22_con_Liver2_1.fastq.gz
D22_con_Liver2_2.fastq.gz
# 2nd run (1 sample; 2 files)
D22_tm_Liver1_1.fastq.gz
D22_tm_Liver1_2.fastq.gz

# List all files with sizes
ls -lah

# Check the integrity of fastq files after download with md5sum
# note to install md5sum, use conda: https://anaconda.org/channels/conda-forge/packages/cms-md5/overview
# note to download md5sum code file, use checksums_liver_training.txt in https://github.com/shadi-shahatit/AtrophyTranscriptomics/tree/main/scripts_training

md5sum -c checksums_liver.txt

cd ../

#################################
## Step 1: Run FastQC on Raw Reads
#################################

fastqc raw_fastq/D22_con_Liver2_1.fastq.gz -o fastqc_reports/
fastqc raw_fastq/D22_con_Liver2_2.fastq.gz -o fastqc_reports/

echo "FastQC finished running!"

# duration = 241 minutes and 7 seconds

#################################
## Step 2: Run Trimmomatic
#################################

# make sure to fix the path of the adapters_truseq_v2.fa

trimmomatic PE -threads 8 -phred33 \
            raw_fastq/D22_con_Liver2_1.fastq.gz raw_fastq/D22_con_Liver2_2.fastq.gz \
            trimmed_fastq/D22_con_Liver2_trim_1.fastq.gz trimmed_fastq/D22_con_Liver2_trim_1.fastq.gz_unpair_1.fastq.gz \
            trimmed_fastq/D22_con_Liver2_trim_2.fastq.gz trimmed_fastq/D22_con_Liver2_trim_2.fastq.gz_unpair_2.fastq.gz \
            ILLUMINACLIP:adapters_truseq_v2.fa:2:30:10:5:true \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Trimmomatic finished running!"

# duration = ~465 minutes and 50 seconds

#################################
## Step 3: Run FastQC on Trimmed Reads
#################################

fastqc trimmed_fastq/D22_con_Liver2_trim_1.fastq.gz -o fastqc_reports/
fastqc trimmed_fastq/D22_con_Liver2_trim_2.fastq.gz -o fastqc_reports/

echo "FastQC finished running!"

# duration = 229 minutes and 58 seconds

#################################
## Step 4: Reference Transcriptome Indexing 
#################################

# Download the ref files from Ensembl - most recent release (115)
# wget https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/cdna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz 	# ref trans - cdna fa release 115

index_dir_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115"
ref_trans_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"

salmon index -t $ref_trans_r155 -i $index_dir_r155

echo "Indexing finished running!"

# duration = 

#################################
## Step 5: Salmon Quantification
#################################

salmon quant -i $index_dir_r155 -l A \
        -1 trimmed_fastq/D22_con_Liver2_trim_1.fastq.gz \
        -2 trimmed_fastq/D22_con_Liver2_trim_2.fastq.gz \
        -p 8 --validateMappings \
        -o quant_r115/D22_con_Liver2_quant

echo "Salmon finished running!"

# duration = 121 minutes and 8 seconds

#################################
## 1.6 Appendix
#################################

# RUN all samples via loops
# BUT make sure to fix the paths to the directories

## Directories setup for the loop

# Define directories
work_dir="/TMBroilers_Transcriptomics/liver" ## change this to your sys path
raw_dir="NewTranscriptiomeProject/liver"
fastqc_dir="${work_dir}/fastqc_reports"
trim_dir="${work_dir}/trimmed_fastq"
quant_dir_r155="${work_dir}/quant_r115"
index_dir_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115"
ref_trans_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"
adapters="/TMBroilers_Transcriptomics/adapters_truseq_v2.fa"

# Create directories
mkdir -p "$fastqc_dir" "$trim_dir" "$quant_dir_r155"

## Step 1: Run FastQC on Raw Reads

fastqc "${raw_dir}"/*.fastq.gz -o "$fastqc_dir"

echo "FastQC finished running!"

## Step 2: Run Trimmomatic

for file1 in ${raw_dir}/*_1.fastq.gz; do
    # Extract sample name prefix
    base=$(basename "$file1" _1.fastq.gz)
    file2="${raw_dir}/${base}_2.fastq.gz"

    if [[ -f "$file2" ]]; then
        echo "Trimming sample: $base ..."
        trimmomatic PE -threads 8 -phred33 \
            "$file1" "$file2" \
            "${trim_dir}/${base}_trim_1.fastq.gz" "${trim_dir}/${base}_unpair_1.fastq.gz" \
            "${trim_dir}/${base}_trim_2.fastq.gz" "${trim_dir}/${base}_unpair_2.fastq.gz" \
            ILLUMINACLIP:${adapters}:2:30:10:5:true \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        echo "Finished sample: $base"
    else
        echo "WARNING: Paired file for $file1 not found. Skipping."
    fi
done

echo "Trimmomatic finished running!"

## Step 3: Run FastQC on Trimmed Reads

fastqc "${trim_dir}"/*_trim_*.fastq.gz -o "$fastqc_dir"

echo "FastQC finished running!"

## Step 5: Salmon Quantification

for fn in ${trim_dir}/*_trim_1.fastq.gz; do
    samp=$(basename "$fn" "_trim_1.fastq.gz")
    echo "Processing sample ${samp}"
    salmon quant -i $index_dir_r155 -l A \
        -1 "${trim_dir}/${samp}_trim_1.fastq.gz" \
        -2 "${trim_dir}/${samp}_trim_2.fastq.gz" \
        -p 8 --validateMappings \
        -o "${quant_dir_r155}/${samp}_quant"
    echo "Finished quantification for ${samp}"
done

echo "Salmon finished running!"


