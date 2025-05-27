#!/bin/bash

ADAPTERS="/home/bandana/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
TRIMMOMATIC_JAR="/home/bandana/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar"
THREADS=8

OUTPUT_DIR="/media/bandana/DATA/Bandana/trans/DEGs/control"
REF_GENOME_DIR="/media/bandana/DATA/Bandana/trans/DEGs/control"
REF_FASTA="/media/bandana/DATA/Bandana/trans/DEGs/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

#STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $REF_GENOME_DIR --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.109.gtf --sjdbOverhang 100 --genomeSAindexNbases 12 --genomeChrBinNbits 18

for R1 in *1.fastq; do
    
    SAMPLE=$(basename "$R1" _1.fastq)
    
    # Define input/output files
    R2="${SAMPLE}_2.fastq"
    
    TRIMMED_R1="${SAMPLE}_1_trimmed.fastq"
    UNPAIRED_R1="${SAMPLE}_1_unpaired.fastq"
    TRIMMED_R2="${SAMPLE}_2_trimmed.fastq"
    UNPAIRED_R2="${SAMPLE}_2_unpaired.fastq"
    
    # Run Trimmomatic
    java -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" -phred33 "$R1" "$R2" "$TRIMMED_R1" "$UNPAIRED_R1" "$TRIMMED_R2" "$UNPAIRED_R2" ILLUMINACLIP:"$ADAPTERS":2:30:10:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    echo "Trimming completed for sample: $SAMPLE"			

    echo "Processing sample: $SAMPLE"
    
    STAR --genomeDir "$REF_GENOME_DIR" --readFilesIn "$TRIMMED_R1" "$TRIMMED_R2" --readFilesCommand cat  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_" --quantMode GeneCounts
    
    grep -E "Uniquely mapped|%|Mapping efficiency" "$OUTPUT_DIR/${SAMPLE}_Log.final.out"
      
    samtools view -T "$REF_FASTA" -C -o "${OUTPUT_DIR}/${SAMPLE}.cram" "${OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    
    echo "Alignment with STAR completed for sample: $SAMPLE"
    
    /home/bandana/Downloads/subread-2.1.1-Linux-x86_64/bin/featureCounts -p -O -T 8 -a Homo_sapiens.GRCh38.109.gtf -o "${OUTPUT_DIR}/${SAMPLE}_featureCounts.txt" "${OUTPUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    
    echo "feature counts completed for sample: $SAMPLE"
 
done
