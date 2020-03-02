#!/bin/bash
#SBATCH --job-name=gogogo    # Job name
#SBATCH --mail-user=chongmin@andrew.cmu.edu     # Where to send mail
#SBATCH --time=10:00:00
#SBATCH -N 1                 # Run on a single CPU
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --ntasks-per-node 20
#SBATCH --output=gogogo.log   # Standard output and error log

. $HOME/.bashrc
conda activate eclip

# star build index step
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir index --genomeFastaFiles genome/GRCh38.p13.genome.fa --sjdbGTFfile genome/gencode.v33.annotation.gtf

# fastqc step
fastqc -o raw_data/rep2/fastqc/before raw_data/rep2/read1.fastq raw_data/rep2/read2.fastq

# cutadpat step
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18     -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA     -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA     -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG     -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG     -A GGAAGAGCGTCGTGT -o raw_data/rep2/cutadapt/round1/read_1_adapterTrim.fastq  -p raw_data/rep2/cutadapt/round1/read_2_adapterTrim.fastq raw_data/rep2/read1.fastq raw_data/rep2/read2.fastq
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18     -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG     -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC     -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT     -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o raw_data/rep2/cutadapt/round2/read_1_adapterTrim.fastq -p raw_data/rep2/cutadapt/round2/read_2_adapterTrim.fastq  raw_data/rep2/cutadapt/round1/read_1_adapterTrim.fastq raw_data/rep2/cutadapt/round1/read_2_adapterTrim.fastq

# fastqc step
fastqc -o raw_data/rep2/fastqc/after raw_data/rep2/cutadapt/round2/read_1_adapterTrim.fastq raw_data/rep2/cutadapt/round2/read_2_adapterTrim.fastq

# prinseq step
prinseq-lite.pl -verbose  -fastq raw_data/rep2/cutadapt/round2/read_1_adapterTrim.fastq -graph_data raw_data/rep2/prinseq/read1.gd -out_good null -out_bad null -graph_stats ld,ns,ts,dn
prinseq-lite.pl -verbose  -fastq raw_data/rep2/cutadapt/round2/read_2_adapterTrim.fastq -graph_data raw_data/rep2/prinseq/read2.gd -out_good null -out_bad null -graph_stats ld,ns,ts,dn

# star adding headers step
awk -v l=10 'BEGIN{OFS=FS= } substr(, 1, 1) == @ {print @ substr(, (l+3), 500) _ substr(, 2, l)    };     substr(, 1, 1) != @ {print}; ' raw_data/rep2/cutadapt/round2/read_1_adapterTrim.fastq > raw_data/rep2/cutadapt/round2/read_1_adapterTrim.bc.fastq
awk -v l=10 'BEGIN{OFS=FS= } substr(, 1, 1) == @ {print @ substr(, (l+3), 500) _ substr(, 2, l)    };     substr(, 1, 1) != @ {print}; ' raw_data/rep2/cutadapt/round2/read_2_adapterTrim.fastq > raw_data/rep2/cutadapt/round2/read_2_adapterTrim.bc.fastq

# star alignment step
STAR --outFilterType BySJout --genomeDir genome --readFilesIn raw_data/rep2/cutadapt/round2/read_1_adapterTrim.bc.fastq,raw_data/rep2/cutadapt/round2/read_2_adapterTrim.bc.fastq --outFileNamePrefix raw_data/rep2/star/rep_     --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999     --outFilterMismatchNoverLmax 0.04 --scoreDelOpen -1 --alignIntronMin 20 --alignIntronMax 1000000     --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --runThreadN 10 --alignEndsType EndToEnd

# umi step
samtools index raw_data/rep2/star/rep_Aligned.sortedByCoord.out.bam
umi_tools dedup -I raw_data/rep2/star/rep_Aligned.sortedByCoord.out.bam --paired -S raw_data/rep2/pre_peak/rep_aligned.f.duplRm.bam

# fastqc step
fastqc -o raw_data/rep1/fastqc/before raw_data/rep1/read1.fastq raw_data/rep1/read2.fastq

# cutadpat step
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18     -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA     -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA     -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG     -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG     -A GGAAGAGCGTCGTGT -o raw_data/rep1/cutadapt/round1/read_1_adapterTrim.fastq  -p raw_data/rep1/cutadapt/round1/read_2_adapterTrim.fastq raw_data/rep1/read1.fastq raw_data/rep1/read2.fastq
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18     -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG     -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC     -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT     -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o raw_data/rep1/cutadapt/round2/read_1_adapterTrim.fastq -p raw_data/rep1/cutadapt/round2/read_2_adapterTrim.fastq  raw_data/rep1/cutadapt/round1/read_1_adapterTrim.fastq raw_data/rep1/cutadapt/round1/read_2_adapterTrim.fastq

# fastqc step
fastqc -o raw_data/rep1/fastqc/after raw_data/rep1/cutadapt/round2/read_1_adapterTrim.fastq raw_data/rep1/cutadapt/round2/read_2_adapterTrim.fastq

# prinseq step
prinseq-lite.pl -verbose  -fastq raw_data/rep1/cutadapt/round2/read_1_adapterTrim.fastq -graph_data raw_data/rep1/prinseq/read1.gd -out_good null -out_bad null -graph_stats ld,ns,ts,dn
prinseq-lite.pl -verbose  -fastq raw_data/rep1/cutadapt/round2/read_2_adapterTrim.fastq -graph_data raw_data/rep1/prinseq/read2.gd -out_good null -out_bad null -graph_stats ld,ns,ts,dn

# star adding headers step
awk -v l=10 'BEGIN{OFS=FS= } substr(, 1, 1) == @ {print @ substr(, (l+3), 500) _ substr(, 2, l)    };     substr(, 1, 1) != @ {print}; ' raw_data/rep1/cutadapt/round2/read_1_adapterTrim.fastq > raw_data/rep1/cutadapt/round2/read_1_adapterTrim.bc.fastq
awk -v l=10 'BEGIN{OFS=FS= } substr(, 1, 1) == @ {print @ substr(, (l+3), 500) _ substr(, 2, l)    };     substr(, 1, 1) != @ {print}; ' raw_data/rep1/cutadapt/round2/read_2_adapterTrim.fastq > raw_data/rep1/cutadapt/round2/read_2_adapterTrim.bc.fastq

# star alignment step
STAR --outFilterType BySJout --genomeDir genome --readFilesIn raw_data/rep1/cutadapt/round2/read_1_adapterTrim.bc.fastq,raw_data/rep1/cutadapt/round2/read_2_adapterTrim.bc.fastq --outFileNamePrefix raw_data/rep1/star/rep_     --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999     --outFilterMismatchNoverLmax 0.04 --scoreDelOpen -1 --alignIntronMin 20 --alignIntronMax 1000000     --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --runThreadN 10 --alignEndsType EndToEnd

# umi step
samtools index raw_data/rep1/star/rep_Aligned.sortedByCoord.out.bam
umi_tools dedup -I raw_data/rep1/star/rep_Aligned.sortedByCoord.out.bam --paired -S raw_data/rep1/pre_peak/rep_aligned.f.duplRm.bam

# fastqc step
fastqc -o raw_data/control/fastqc/before raw_data/control/read1.fastq raw_data/control/read2.fastq

# cutadpat step
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18     -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA     -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA     -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG     -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG     -A GGAAGAGCGTCGTGT -o raw_data/control/cutadapt/round1/read_1_adapterTrim.fastq  -p raw_data/control/cutadapt/round1/read_2_adapterTrim.fastq raw_data/control/read1.fastq raw_data/control/read2.fastq
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18     -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG     -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC     -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT     -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o raw_data/control/cutadapt/round2/read_1_adapterTrim.fastq -p raw_data/control/cutadapt/round2/read_2_adapterTrim.fastq  raw_data/control/cutadapt/round1/read_1_adapterTrim.fastq raw_data/control/cutadapt/round1/read_2_adapterTrim.fastq

# fastqc step
fastqc -o raw_data/control/fastqc/after raw_data/control/cutadapt/round2/read_1_adapterTrim.fastq raw_data/control/cutadapt/round2/read_2_adapterTrim.fastq

# prinseq step
prinseq-lite.pl -verbose  -fastq raw_data/control/cutadapt/round2/read_1_adapterTrim.fastq -graph_data raw_data/control/prinseq/read1.gd -out_good null -out_bad null -graph_stats ld,ns,ts,dn
prinseq-lite.pl -verbose  -fastq raw_data/control/cutadapt/round2/read_2_adapterTrim.fastq -graph_data raw_data/control/prinseq/read2.gd -out_good null -out_bad null -graph_stats ld,ns,ts,dn

# star adding headers step
awk -v l=10 'BEGIN{OFS=FS= } substr(, 1, 1) == @ {print @ substr(, (l+3), 500) _ substr(, 2, l)    };     substr(, 1, 1) != @ {print}; ' raw_data/control/cutadapt/round2/read_1_adapterTrim.fastq > raw_data/control/cutadapt/round2/read_1_adapterTrim.bc.fastq
awk -v l=10 'BEGIN{OFS=FS= } substr(, 1, 1) == @ {print @ substr(, (l+3), 500) _ substr(, 2, l)    };     substr(, 1, 1) != @ {print}; ' raw_data/control/cutadapt/round2/read_2_adapterTrim.fastq > raw_data/control/cutadapt/round2/read_2_adapterTrim.bc.fastq

# star alignment step
STAR --outFilterType BySJout --genomeDir genome --readFilesIn raw_data/control/cutadapt/round2/read_1_adapterTrim.bc.fastq,raw_data/control/cutadapt/round2/read_2_adapterTrim.bc.fastq --outFileNamePrefix raw_data/control/star/rep_     --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999     --outFilterMismatchNoverLmax 0.04 --scoreDelOpen -1 --alignIntronMin 20 --alignIntronMax 1000000     --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --runThreadN 10 --alignEndsType EndToEnd

# umi step
samtools index raw_data/control/star/rep_Aligned.sortedByCoord.out.bam
umi_tools dedup -I raw_data/control/star/rep_Aligned.sortedByCoord.out.bam --paired -S raw_data/control/pre_peak/rep_aligned.f.duplRm.bam

# samtools merge step
samtools merge -f samtools/rep1_2_sortedByCoord_merged.bam raw_data/rep2/pre_peak/rep_Aligned.sortedByCoord.out.bam raw_data/rep1/pre_peak/rep_Aligned.sortedByCoord.out.bam
samtools index samtools/rep1_2_sortedByCoord_merged.bam
mv raw_data/control/pre_peak/rep_aligned.f.duplRm.bam samtools/control_sortedByCoord.bam
mv raw_data/control/pre_peak/rep_aligned.f.duplRm.bam.bai samtools/control_sortedByCoord.bam.bai

# pureclip step
pureclip -i samtools/rep1_2_sortedByCoord_merged.bam -bai samtools/rep1_2_sortedByCoord_merged.bam.bai -g genome/GRCh38.p13.genome.fa -iv 'chr1;chr2;chr3;' -nt 10 -o peak/rep1_2_pureclip_crosslink_sites.bed

# pureclip step
pureclip -i samtools/control_sortedByCoord.bam -bai samtools/control_sortedByCoord.bam.bai -g genome/GRCh38.p13.genome.fa -iv 'chr1;chr2;chr3;' -nt 10 -o peak/control_pureclip_crosslink_sites.bed
