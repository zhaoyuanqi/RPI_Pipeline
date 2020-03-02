from glob import glob
from os.path import join
from os import makedirs

import os

def de_analysis(config, control_tsv, knock_tsv, output_de_path, output_de_gene_list_path):

    cmd = []
    input_tsv = control_tsv + knock_tsv
    input_str = ' '.join(input_tsv)
    tmp_dir = 'tmp'
    makedirs(tmp_dir, exist_ok=True)
    matrix_path = join(tmp_dir, 'result.matrix')
    os.system('rsem-generate-data-matrix {} > {}'.format(input_str, matrix_path))
    ebseq_path = join(tmp_dir, 'result.ebseq')
    cmd.append('rsem-run-ebseq {} {},{} {}'.format(matrix_path, len(control_tsv), len(knock_tsv), ebseq_path))
    cmd.append('rsem-control-fdr {} {} {}'.format(ebseq_path, config.de.fdr, output_de_path))
    cmd.append("awk 'NR!=1{{ print $1}}' {} > {}".format(output_de_path, output_de_gene_list_path))
    cmd = '\n'.join(cmd)
    os.system('echo "{}" > de.sh'.format(cmd))
    os.system('sh de.sh')
    return

def de(config):

    if not config.de.build:
        print('differential expression analysis was turned off...')
        print('set config.de.build to be True in your config file if you want to do differential expression analysis')
        return

    shrna_dir = join(config.raw.raw_dir, 'shrna')
    control_tsv_list = glob(join(shrna_dir, 'control_rep*.tsv'))
    knock_tsv_list = glob(join(shrna_dir, 'rep*.tsv'))
    fdr = config.de.fdr
    de_dir = config.de.de_dir
    makedirs(de_dir, exist_ok=True)
    output_de_path = join(de_dir, 'shrna_de_gene_fdr_{}.txt'.format(fdr))
    output_de_gene_list_path = join(de_dir, 'shrna_de_gene_list.txt')
    de_analysis(config, control_tsv_list, knock_tsv_list, output_de_path, output_de_gene_list_path)
    print('Finished Differential expression analysis...')
    return


