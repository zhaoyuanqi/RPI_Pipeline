import requests
import pandas as pd
import numpy as np
from tqdm import tqdm

from os.path import join, exists, basename
from os import makedirs, remove
import os

def download_single_file(url, save_path):

    if exists(save_path):
        print("{} has already been downloaded".format(save_path))
        return True

    print('Downloading file from {}'.format(url))
    r = requests.get(url, allow_redirects=True, stream=True)
    total_size = int(r.headers.get("content-length", 0))
    block_size = 1024
    t = tqdm(total=total_size, unit="iB", unit_scale=True)
    with open(save_path, "wb") as f:
        for data in r.iter_content(block_size):
            t.update(len(data))
            f.write(data)
    t.close()
    if total_size != 0 and t.n != total_size:
        remove(save_path)
        return False
    else:
        return True


def download(config):

    ## download meta file
    base_url = "https://www.encodeproject.org/metadata/?type=Experiment&searchTerm={}"
    makedirs(config.raw.raw_dir, exist_ok=True)
    metafile_path = join(config.raw.raw_dir, '{}.tsv'.format(config.gene))

    metafile_url = base_url.format(config.gene)
    if not download_single_file(metafile_url, metafile_path):
        raise Exception("Failed to download metafile...")

    meta_tb = pd.read_csv(metafile_path, sep='\t')
    if meta_tb.shape[0] == "0":
        raise Exception("Empty metafile!!!")
    else:
        print("Successfully downloaded metafile".format(config.gene))

    # download eclip files
    eclip_fastq_tb = meta_tb[(meta_tb['Assay'] == 'eCLIP') & (meta_tb['File format'] == 'fastq')]
    eclip_cell = set(eclip_fastq_tb['Biosample term name'].values)
    if config.cell not in eclip_cell:
        print('Available cell line with eCLIP experiments: {}'.format(eclip_cell))
        raise Exception('No eCLIP experiment with {}'.format(config.cell))
    else:
        eclip_fastq_tb = eclip_fastq_tb[eclip_fastq_tb['Biosample term name'] == config.cell]

    eclip_rep = set(eclip_fastq_tb['Biological replicate(s)'].values)
    print('{} has {} eCLIP biological replicate(s)'.format(config.gene, len(eclip_rep)))

    eclip_control_fastq_tb = eclip_fastq_tb[eclip_fastq_tb['Controlled by'].isna()]
    eclip_fastq_tb = eclip_fastq_tb[eclip_fastq_tb['Controlled by'].notna()]


    for r in eclip_rep:
        r_dir = join(config.raw.raw_dir, 'rep{}'.format(r))
        makedirs(r_dir, exist_ok=True)
        r_reads = eclip_fastq_tb[eclip_fastq_tb['Biological replicate(s)'] == r][['File download URL', 'Paired end']].values
        for r_url, r_read_num in r_reads:
            r_fastq_path = join(r_dir, 'read{}.fastq.gz'.format(int(r_read_num)))
            if not download_single_file(r_url, r_fastq_path):
                raise Exception("Failed to download {}", basename(r_fastq_path))
            print('unzip {}...'.format(r_fastq_path))
            os.system('gunzip {}'.format(r_fastq_path))

    eclip_control_dir = join(config.raw.raw_dir, 'control')
    makedirs(eclip_control_dir, exist_ok=True)
    eclip_control_reads = eclip_control_fastq_tb[['File download URL', 'Paired end']].values
    for r_url, r_read_num in eclip_control_reads:
        r_fastq_path = join(eclip_control_dir, 'read{}.fastq.gz'.format(int(r_read_num)))
        if not download_single_file(r_url, r_fastq_path):
            raise Exception("Failed to download {}", basename(r_fastq_path))
        print('unzip {}...'.format(r_fastq_path))
        os.system('gunzip {}'.format(r_fastq_path))

    print('Finish downloading eCLIP data...')

    ## download shRNA files
    shrna_tsv_tb = meta_tb[(meta_tb['Assay'] == 'shRNA knockdown followed by RNA-seq') & (meta_tb['File format'] == 'tsv') &\
                           (meta_tb['Output type'] == 'gene quantifications') & (meta_tb['Biosample term name'] == config.cell) &\
                           (meta_tb['Assembly'] == config.assembly)]
    shrna_rep = set(shrna_tsv_tb['Biological replicate(s)'].values)
    shrna_dir = join(config.raw.raw_dir, 'shrna')
    makedirs(shrna_dir, exist_ok=True)
    for r in shrna_rep:
        r_tsv_path = join(shrna_dir, 'rep{}.tsv'.format(r))
        r_url = shrna_tsv_tb[shrna_tsv_tb['Biological replicate(s)'] == r]['File download URL'].values[0]
        if not download_single_file(r_url, r_tsv_path):
            raise Exception("Failed to download {}", basename(r_tsv_path))

    ## download shRNA control files
    shrna_fastq_tb = meta_tb[
        (meta_tb['Assay'] == 'shRNA knockdown followed by RNA-seq') & (meta_tb['File format'] == 'fastq') & \
        (meta_tb['Biosample term name'] == config.cell)]
    shrna_control = basename(shrna_fastq_tb['Controlled by'].values[0].split(',')[0][:-1])
    shrna_control_metafile_url = base_url.format(shrna_control)
    shrna_control_path = join(shrna_dir, '{}_shrna_control.tsv'.format(config.gene))
    if not download_single_file(shrna_control_metafile_url, shrna_control_path):
        raise Exception("Failed to download {}", basename(shrna_control_path))

    shrna_control_tb = pd.read_csv(shrna_control_path, sep='\t')
    shrna_control_tsv_tb = shrna_control_tb[(shrna_control_tb['Assay'] == 'shRNA knockdown followed by RNA-seq') & (shrna_control_tb['File format'] == 'tsv') &\
                                            (shrna_control_tb['Output type'] == 'gene quantifications') & (shrna_control_tb['Biosample term name'] == config.cell) &\
                                            (shrna_control_tb['Assembly'] == config.assembly)]

    shrna_control_rep = set(shrna_control_tsv_tb['Biological replicate(s)'].values)
    for r in shrna_control_rep:
        r_tsv_path = join(shrna_dir, 'control_rep{}.tsv'.format(r))
        r_url = shrna_control_tsv_tb[shrna_control_tsv_tb['Biological replicate(s)'] == r]['File download URL'].values[0]
        if not download_single_file(r_url, r_tsv_path):
            raise Exception("Failed to download {}", basename(r_tsv_path))



    print('Finish downloading shRNA data...')

    genome_dir = config.genome.genome_dir
    genome_fasta_path = join(genome_dir, config.genome.fasta_path)
    gtf_path = join(genome_dir, config.genome.gtf_path)
    os.system('wget -O {} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.p13.genome.fa.gz'.format(genome_fasta_path))
    os.system('wget -O {} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz'.format(gtf_path))
    os.system('gunzip {}'.format(genome_fasta_path))
    os.system('gunzip {}'.format(gtf_path))
    print('Finish downloading genome...')
    return