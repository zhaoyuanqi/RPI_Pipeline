from glob import glob
from os.path import join, basename
from os import makedirs
import os

def fastqc(read1_path, read2_path, output_path):
    prompt = '\n# fastqc step\n'
    cmd_str = 'fastqc -o {} {} {}'.format(output_path, read1_path, read2_path)
    return prompt + cmd_str


def cutadapt(round1_cut1_path, round1_cut2_path, round2_cut1_path, round2_cut2_path, read1_path, read2_path):
    prompt = '\n# cutadpat step\n'
    cmd_round1_str = 'cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
    -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA \
    -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA \
    -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG \
    -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG \
    -A GGAAGAGCGTCGTGT -o {}  -p {} {} {}'.format(round1_cut1_path, round1_cut2_path, read1_path, read2_path)

    cmd_round2_str = 'cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 \
    -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG \
    -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC \
    -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT \
    -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT -o {} -p {}  {} {}'.format(round2_cut1_path, round2_cut2_path, round1_cut1_path, round1_cut2_path)

    cmd_str = cmd_round1_str+'\n'+cmd_round2_str
    return prompt + cmd_str

def prinseq(round2_cut1_path, round2_cut2_path, out_graph1_path, out_graph2_path):
    prompt = '\n# prinseq step\n'
    base_str = 'prinseq-lite.pl -verbose  -fastq {} -graph_data {} -out_good null -out_bad null -graph_stats ld,ns,ts,dn'
    cmd_str1 = base_str.format(round2_cut1_path, out_graph1_path)
    cmd_str2 = base_str.format(round2_cut2_path, out_graph2_path)
    cmd_str = cmd_str1+'\n'+cmd_str2
    return prompt + cmd_str


def star_build_index(index_dir, genome_fasta_path, gtf_path):
    prompt = '\n# star build index step\n'
    cmd_str = 'STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {}'.format(index_dir, genome_fasta_path, gtf_path)
    return prompt + cmd_str

def star(round2_cut1_path, round2_cut2_path, bc1_path, bc2_path, output_prefix, index_dir):
    prompt = '\n# star adding headers step\n'
    ## awk step
    bc_base_str = 'awk -v l=10 \'BEGIN{{OFS=FS=" "}} substr($1, 1, 1) == "@" {{print "@" substr($1, (l+3), 500) "_" substr($1, 2, l) " " $2 }}; \
    substr($1, 1, 1) != "@" {{print}}; \' {} > {}'
    bc_cmd_str1 = bc_base_str.format(round2_cut1_path, bc1_path)
    bc_cmd_str2 = bc_base_str.format(round2_cut2_path, bc2_path)
    bc_cmd_str = prompt + bc_cmd_str1+'\n'+bc_cmd_str2

    prompt = '\n# star alignment step\n'
    ## star step
    star_base_str = 'STAR --outFilterType BySJout --genomeDir {} --readFilesIn {},{} --outFileNamePrefix {} \
    --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 --scoreDelOpen -1 --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --runThreadN 10 --alignEndsType EndToEnd'
    star_cmd_str = prompt + star_base_str.format(index_dir, bc1_path, bc2_path, output_prefix)

    cmd_str = bc_cmd_str+'\n'+star_cmd_str
    return cmd_str


def umi(bam_path, output_bam_path):
    prompt = '\n# umi step\n'
    ## samtools step
    samtools_cmd_str = 'samtools index {}'.format(bam_path)

    ## umi step
    umi_cmd_str = 'umi_tools dedup -I {} --paired -S {}'.format(bam_path, output_bam_path)

    cmd_str = samtools_cmd_str+'\n'+umi_cmd_str
    return prompt + cmd_str

def samtools_index(bam_path):
    prompt = '\n# samtools index step\n'
    cmd_str = 'samtools index {}'.format(bam_path)
    return cmd_str

def samtools_merge(input1_path, input2_path, output_path):
    prompt = '\n# samtools merge step\n'
    cmd_str1 = 'samtools merge -f {} {} {}'.format(output_path, input1_path, input2_path)
    cmd_str2 = samtools_index(output_path)
    cmd_str = cmd_str1+'\n'+cmd_str2
    return prompt + cmd_str


def pureclip(bam_path, bai_path, genome_path, bed_path):
    prompt = '\n# pureclip step\n'
    cmd_str = "pureclip -i {} -bai {} -g {} -iv 'chr1;chr2;chr3;' -nt 10 -o {}".format(bam_path, bai_path, genome_path, bed_path)
    return prompt + cmd_str


def sbatch_generator(cmd, sbatch_template='sbatch_template.sh'):

    os.system('cp {} {}'.format(sbatch_template, 'sbatch.sh'))
    os.system('echo "{}" >> sbatch.sh'.format(cmd))
    os.system('sbatch sbatch.sh')

def read_process(config):

    ## star build index

    cmd = []

    index_dir = config.genome.index_dir
    genome_dir = config.genome.genome_dir
    makedirs(index_dir, exist_ok=True)
    makedirs(genome_dir, exist_ok=True)
    genome_fasta_path = join(genome_dir, config.genome.fasta_path)
    if config.genome.build:
        gtf_path = join(genome_dir, config.genome.gtf_path)
        cmd_tmp = star_build_index(index_dir, genome_fasta_path, gtf_path)
        cmd.append(cmd_tmp)
        # sbatch_generator(cmd_tmp)

    rep_dirs = glob(join(config.raw.raw_dir, 'rep*'))
    control_dir = glob(join(config.raw.raw_dir, 'control'))
    dirs = rep_dirs + control_dir
    for dir in dirs:

        read1_path = join(dir, 'read1.fastq')
        read2_path = join(dir, 'read2.fastq')

        ## fastqc before
        fastqc_before_dir = join(dir, 'fastqc', 'before')
        makedirs(fastqc_before_dir, exist_ok=True)
        cmd_tmp = fastqc(read1_path, read2_path, fastqc_before_dir)
        # sbatch_generator(cmd_tmp)
        cmd.append(cmd_tmp)

        ## cutadpat
        cutadapt_round1_dir = join(dir, 'cutadapt', 'round1')
        cutadapt_round2_dir = join(dir, 'cutadapt', 'round2')
        makedirs(cutadapt_round1_dir, exist_ok=True)
        makedirs(cutadapt_round2_dir, exist_ok=True)
        round1_cut1_path = join(cutadapt_round1_dir, 'read_1_adapterTrim.fastq')
        round1_cut2_path = join(cutadapt_round1_dir, 'read_2_adapterTrim.fastq')
        round2_cut1_path = join(cutadapt_round2_dir, 'read_1_adapterTrim.fastq')
        round2_cut2_path = join(cutadapt_round2_dir, 'read_2_adapterTrim.fastq')
        cmd_tmp = cutadapt(round1_cut1_path, round1_cut2_path, round2_cut1_path, round2_cut2_path, read1_path, read2_path)
        # sbatch_generator(cmd_tmp)
        cmd.append(cmd_tmp)

        ## fastqc after
        fastqc_after_dir = join(dir, 'fastqc', 'after')
        makedirs(fastqc_after_dir, exist_ok=True)
        cmd_tmp = fastqc(round2_cut1_path, round2_cut2_path, fastqc_after_dir)
        # sbatch_generator(cmd_tmp)
        cmd.append(cmd_tmp)

        ## prinseq
        prinseq_dir = join(dir, 'prinseq')
        makedirs(prinseq_dir, exist_ok=True)
        out_graph1_path = join(prinseq_dir, 'read1.gd')
        out_graph2_path = join(prinseq_dir, 'read2.gd')
        cmd_tmp = prinseq(round2_cut1_path, round2_cut2_path, out_graph1_path, out_graph2_path)
        # sbatch_generator(cmd_tmp)
        cmd.append(cmd_tmp)

        ## star
        star_dir = join(dir, 'star')
        makedirs(star_dir, exist_ok=True)
        bc1_path = join(cutadapt_round2_dir, 'read_1_adapterTrim.bc.fastq')
        bc2_path = join(cutadapt_round2_dir, 'read_2_adapterTrim.bc.fastq')
        output_prefix = join(star_dir, 'rep_')
        cmd_tmp = star(round2_cut1_path, round2_cut2_path, bc1_path, bc2_path, output_prefix, genome_dir)
        # sbatch_generator(cmd_tmp)
        cmd.append(cmd_tmp)

        ## umi
        pre_preak_dir = join(dir, 'pre_peak')
        makedirs(pre_preak_dir, exist_ok=True)
        bam_path = output_prefix+'Aligned.sortedByCoord.out.bam'
        output_bam_path = join(pre_preak_dir, 'rep_aligned.f.duplRm.bam')
        cmd_tmp = umi(bam_path, output_bam_path)
        # sbatch_generator(cmd_tmp)
        cmd.append(cmd_tmp)

    ## samtools merge rep1 and rep2
    samtools_dir = config.samtools.samtools_dir
    makedirs(samtools_dir, exist_ok=True)
    merged_bam_path = join(samtools_dir, 'rep1_2_sortedByCoord_merged.bam')
    bam_paths = [join(rep_dir, 'pre_peak', 'rep_Aligned.sortedByCoord.out.bam') for rep_dir in rep_dirs]
    cmd_tmp = samtools_merge(bam_paths[0], bam_paths[1], merged_bam_path)
    # sbatch_generator(cmd_tmp)
    cmd.append(cmd_tmp)

    ## samtools index control
    control_bam_path = join(control_dir[0], 'pre_peak', 'rep_aligned.f.duplRm.bam')
    control_bai_path = control_bam_path+'.bai'
    samtools_index(bam_path)
    cmd_tmp1 = 'mv {} {}'.format(control_bam_path, join(samtools_dir, 'control_sortedByCoord.bam'))
    cmd_tmp2 = 'mv {} {}'.format(control_bai_path, join(samtools_dir, 'control_sortedByCoord.bam.bai'))
    cmd_tmp = cmd_tmp1+'\n'+cmd_tmp2
    cmd.append(cmd_tmp)


    peak_dir = config.peak.peak_dir
    makedirs(peak_dir, exist_ok=True)
    ## pureclip for merged rep_1_2
    peak_path = join(peak_dir, 'rep1_2_pureclip_crosslink_sites.bed')
    bai_path = join(samtools_dir, 'rep1_2_sortedByCoord_merged.bam.bai')
    cmd_tmp = pureclip(merged_bam_path, bai_path, genome_fasta_path, peak_path)
    # sbatch_generator(cmd_tmp)
    cmd.append(cmd_tmp)

    ## pureclip for control
    peak_path = join(peak_dir, 'control_pureclip_crosslink_sites.bed')
    control_bam_path = join(samtools_dir, 'control_sortedByCoord.bam')
    control_bai_path = join(samtools_dir, 'control_sortedByCoord.bam.bai')
    cmd_tmp = pureclip(control_bam_path, control_bai_path, genome_fasta_path, peak_path)
    # sbatch_generator(cmd_tmp)
    cmd.append(cmd_tmp)

    cmd = '\n'.join(cmd)
    # print(cmd)
    sbatch_generator(cmd)
    print('Uploaded reads processing job to server...')