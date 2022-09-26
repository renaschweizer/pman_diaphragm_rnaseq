# pman_diaphragm_rnaseq
This repository is for data processing and analysis related to Schweizer et al. ([submitted](https://doi.org/10.1101/2022.09.24.509328)), "Gene regulatory changes underlie developmental plasticity in respiration and aerobic performance in highland deer mice" with Catherine M. Ivy, Chandrasekhar Natarajan, Graham R. Scott, Jay F. Storz, and Zachary A. Cheviron. 

Below, users will find example commands for analyzing raw RNAseq fastq data to generate count data. In the .Rmarkdown file, users will find all the analyses to perform gene expression analyses in R, including differential gene expression analysis with EdgeR, identifying regulatory modules in WGCNA, testing for correlations of WGCNA modules with phenotypes, testing for effects of treatment on module expression, and identifying overap with candidate genes for high-altitude adaptation. 

# Processing raw RNAseq fastq reads

1. Run BBduk to remove adapters from paired end RNAseq data. 

```{bash}
ls *_R1_001.fastq.gz | while read file; do
        name=$(echo "${file}" | cut -d "/" -f 8 | cut -d "." -f 1 | cut -d "_" -f 1)
        echo "${name}"                                                                  
        if [ -f "fastq_trimmed/${name}.1.fq.gz" ]; then
                echo "exists; skipping..."                              #
        else
               bbmap/bbduk.sh -Xmx1g \
                in1="${name}_R1_001.fastq.gz" \
                in2="${name}_R2_001.fastq.gz" \
                out1="fastq_trimmed/${name}.1.bb.fq.gz" \
                out2="fastq_trimmed/${name}.2.bb.fq.gz" \
                ref=adapters.fa \
                ktrim=r k=23 mink=11 hdist=1 tpe tbo
        fi
done         
echo "Done!"

```
2. Run trimmomatic to trim junk off of reads. 

```{unix}
ls *.1.bb.fq.gz | while read file; do 
        name=$(echo "${file}" | cut -d "." -f 1) 		
        echo "${name}" 
        if [ -f "${name}.fastp.trimmed_1P.fq.gz" ]; then 		
                echo "exists; skipping..."				
        else

                trimmomatic PE -threads 8 -phred33 -trimlog \
                	"${name}.trimmomatic.log" "${file}" "${name}.2.bb.fq.gz" \
                	-baseout "${name}.fastp.trimmed.fq.gz" \
                	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:36 
        fi
done
echo "Done!"
```
3. Align reads to deer mouse genome with HiSat2. 

```{unix}

ls *.trimmed_1P.fq.gz | while read file; do 
        name=$(echo "${file}" | cut -d "." -f 1) 		
        echo "${name}" 									
	~/hisat2-2.2.1/hisat2 -p 48 --mp 2,0 -q -x ~/pman2_ref/GCA_003704035.1_HU_Pman_2.1_genomic_wMt.fa -1 "${name}.fastp.trimmed_1P.fq.gz" -2 "${name}.fastp.trimmed_2P.fq.gz" -U "${name}.fastp.trimmed_1U.fq.gz" --summary-file "${name}_align_mp20.summ.txt" --un "${name}.mp20.unmapped.fq" | samtools view -Sbo "../bam/${name}_Halign_wMt_mp20.bam" 
	samtools sort ../bam/${name}_Halign_wMt_mp20.bam -o ../bam/${name}_Halign_wMt_mp20.s.bam
done
echo "Done!" 

```

4. Generate gene count data with featureCounts. 

```{unix}
featureCounts -p -O -F GTF -a ~/pman2_ref/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.104.gtf -o pman_diaphragm_wMt_mp20_readcounts.txt *_Halign_wMt_mp20.bam
```

# Analyzing gene expression count data in R

See diaphragm_RNAseq_analysis.Rmd
