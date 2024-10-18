#!/bin/bash -l
#SBATCH -A [project_code]
#SBATCH -J kallisto
#SBATCH -t 20:00:00
#SBATCH -p core
#SBATCH -n 6
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alikarl@kth.se

module purge
module load bioinfo-tools
module load kallisto/0.48.0
module load TrimGalore/0.6.1

# Specify directories
rawdata=data/fastq
trimmed=data/trimmed
kallisto_output=data/kallisto_out
mkdir ${trimmed}
mkdir ${kallisto_output}

# Inputs - ref used from https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/
kallisto_ref="Homo_sapiens.GRCh38.cdna.all.fa.gz"
kallisto_idx="Homo_sapiens.GRCh38.cdna.all_kallisto_index.idx"
meta_file="metadata.txt"
log="kallisto.log"


# 1. CREATE KALLISTO IDX (unless it already exists)
if [ ! -f "$kallisto_idx" ]; then
    echo "Kallisto idx not found -> creating idx"
    kallisto index -i $kallisto_idx $kallisto_ref
else
    echo "Kallisto idx already exists"
fi
wait

# 2. TRIM GALORE
echo "Initialize trim galore: $(date)" >> $log

 Paired files: *R1_001.fastq.gz and *R2_001.fastq.gz
 https://www.biostars.org/p/359116/
ls ${rawdata}/*R1_001.fastq.gz | xargs -P6 -I@ bash -c 'trim_galore --fastqc -o data/trimmed --paired "$1" ${1%R1_001.fastq.gz}R2_001.fastq.gz' _ @
wait

echo "Trim galore finished: $(date)" >> $log &&
echo "Initialize kallisto: $(date)" >> $log &&

# 3. KALLISTO
for sample in $(cut -f1 ${meta_file} | tail -n +2); do
 echo "Kallisto processing sample: $sample" >> $log

 kallisto quant -i ${kallisto_idx} \
        -o ${kallisto_output}/${sample} \
        ${trimmed}/${sample}*R1_001_val_1.fq.gz \
        ${trimmed}/${sample}*R2_001_val_2.fq.gz
done
wait

echo "Kallisto is finished: $(date)" >> $log


