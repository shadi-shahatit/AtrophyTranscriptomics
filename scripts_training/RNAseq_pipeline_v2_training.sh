################################# Transcriptomics Analysis Training - Al-Zghoul Lab, JUST
################################# Shadi Shahatit - RA, 2026

#### Table of Contents

#### 1. RNA-seq Quantification via Salmon 
####	1.0	Setup Environment and Prepare Directories
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
## Step 0: Setup Environment and Prepare Directories
#################################

SECONDS=0

# Create and name conda environment
conda create -y -n salmon
# Activate conda env 
source activate salmon
# install fastqc, trimmomatic, and salmon from the links in the slides

# Define directories
work_dir="/TMBroilers_Transcriptomics/muscle"
adapters="/TMBroilers_Transcriptomics/adapters_truseq_v2.fa"

raw_dir="NewTranscriptiomeProject/muscle"
fastqc_dir="${work_dir}/fastqc_reports"
trim_dir="${work_dir}/trimmed_fastq"

quant_dir_r155="${work_dir}/quant_r115"
index_dir_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115"
ref_trans_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"

mkdir -p "$fastqc_dir" "$trim_dir" "$quant_dir_r155"

#################################
## Step 1: Run FastQC on Raw Reads
#################################

fastqc "${raw_dir}"/*.fastq.gz -o "$fastqc_dir"

echo "FastQC finished running!"

# duration = 241 minutes and 7 seconds

#################################
## Step 2: Run Trimmomatic
#################################

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

# duration = ~465 minutes and 50 seconds

# If confused, run samples one by one instead of using a loop in Appendix

#################################
## Step 3: Run FastQC on Trimmed Reads
#################################

fastqc "${trim_dir}"/*_trim_*.fastq.gz -o "$fastqc_dir"

echo "FastQC finished running!"

# duration = 229 minutes and 58 seconds

#################################
## Step 4: Reference Transcriptome Indexing 
#################################

# Download the ref files from Ensembl - most recent release (115)
# wget https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/cdna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz 	# ref trans - cdna fa release 115

SECONDS=0

work_dir="/TMBroilers_Transcriptomics/muscle"
trim_dir="${work_dir}/trimmed_fastq"
quant_dir_r155="${work_dir}/quant_r115"
index_dir_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115"
ref_trans_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"

mkdir -p "$quant_dir_r155"

salmon index -t $ref_trans_r155 -i $index_dir_r155

echo "Indexing finished running!"

# Pipeline duration
duration=$SECONDS
echo "Pipeline completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."

# duration = 

#################################
## Step 5: Salmon Quantification
#################################

SECONDS=0

work_dir="/TMBroilers_Transcriptomics/muscle"
trim_dir="${work_dir}/trimmed_fastq"
quant_dir_r155="${work_dir}/quant_r115"
index_dir_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115"
ref_trans_r155="/TMBroilers_Transcriptomics/Ggallus_Refs/cdna_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz"

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

# Pipeline duration
duration=$SECONDS
echo "Pipeline completed in $(($duration / 60)) minutes and $(($duration % 60)) seconds."

# duration = 121 minutes and 8 seconds

#################################
## 1.6 Appendix
#################################

# run samples one by one instead of using a loop:
# BUT make sure to fix the paths to the dir

# sample names examples

sample1_1.fastq.gz
sample1_2.fastq.gz

# trimmomatic

trimmomatic PE -threads 8 -phred33 \
            sample1_1.fastq.gz sample1_2.fastq.gz \
            trimmed_fastq/sample1_trim_1.fastq.gz trimmed_fastq/sample1_trim_1.fastq.gz_unpair_1.fastq.gz \
            trimmed_fastq/sample1_trim_2.fastq.gz trimmed_fastq/sample1_trim_2.fastq.gz_unpair_2.fastq.gz \
            ILLUMINACLIP:adapters_truseq_v2.fa:2:30:10:5:true \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# fastqc

fastqc sample1_1.fastq.gz -o fastqc_reports/
fastqc sample1_2.fastq.gz -o fastqc_reports/
fastqc sample1_trim_1.fastq.gz -o fastqc_reports/
fastqc sample1_trim_2.fastq.gz -o fastqc_reports/

# salmon

salmon quant -i $index_dir -l A \
        -1 sample1_trim_1.fastq.gz \
        -2 sample1_trim_2.fastq.gz \
        -p 8 --validateMappings \
        -o output_files


