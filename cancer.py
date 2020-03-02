from os.path import join
from os import makedirs
import os


def cancer(config):

    if not config.cancer.build:
        print('cancer analysis was turned off...')
        print('set config.cancer.build to be True in your config file if you want to do cancer analysis')
        return

    cmd = []

    cancer_dir = config.cancer.cancer_dir
    makedirs(cancer_dir, exist_ok=True)
    shrna_de_gene_path = join(config.de.de_dir, 'shrna_de_gene_list.txt')
    control_bed_path = join(config.peak.peak_dir, 'rep1_2_pureclip_crosslink_sites.bed')
    rep1_2_bed_path = join(config.peak.peak_dir, 'control_pureclip_crosslink_sites.bed')
    bed_output_path = join(cancer_dir, 'eclip_crosslink_gene_name.txt')
    shrna_output_path = join(cancer_dir, 'shrna_de_gene_name.txt')
    cmd.append('Rscript convertBedToGeneNames.R {} {} {} '.format(control_bed_path, rep1_2_bed_path, bed_output_path))
    cmd.append('Rscript ensemblToHGNC.R {} {} '.format(shrna_de_gene_path, shrna_output_path))
    cmd.append('Rscript hypergeometricPval.R {} {} {}'.format(bed_output_path, shrna_output_path, config.cancer.sample_path))

    bed_output_path_2 = join(cancer_dir, 'eclip_cancer.csv')
    shrna_output_path_2 = join(cancer_dir, 'shrna_cancer.csv')
    bed_info_output_path = join(cancer_dir, 'eclip_cancer.txt')
    shrna_info_output_path = join(cancer_dir, 'shrna_cancer.txt')
    cmd.append('python -W ignore Predict_Diff_Expression.py --input-file={} --output-file={} > {}'.format(bed_output_path, bed_output_path_2, bed_info_output_path))
    cmd.append('python -W ignore Predict_Diff_Expression.py --input-file={} --output-file={} > {}'.format(shrna_output_path, shrna_output_path_2, shrna_info_output_path))
    cmd = '\n'.join(cmd)
    os.system('echo "{}" > cancer.sh'.format(cmd))
    os.system('sh cancer.sh')
    print('Finished cancer analysis...')
    return