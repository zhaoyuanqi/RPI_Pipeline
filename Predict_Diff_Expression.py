#!/usr/bin/env python
# coding: utf-8

import collections
import csv
import numpy as np
import os
from scipy.stats import ttest_ind
from scipy.stats import pearsonr
from scipy.stats import shapiro
from scipy.stats import normaltest
from scipy.stats import kendalltau
from scipy.stats import spearmanr
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from scipy.stats import f_oneway
from sklearn.cluster import KMeans
import sys


# Student t distribution
def test_t(x,y,p_threshold=0.01):
    stat, p = ttest_ind(x, y)
    return (p <= p_threshold)

# Anova
def test_anova(x,y,p_threshold=0.01):
    stat, p = f_oneway(x, y)
    return (p <= p_threshold)

# Mann-Whitney U
def test_mannwhitneyu(x,y,p_threshold=0.01):
    stat, p = mannwhitneyu(x, y)
    return (p <= p_threshold)

# Kruskal-Wallis H 
def test_kruskal(x,y,p_threshold=0.01):
    stat, p = kruskal(x, y)
    return (p <= p_threshold)

# Pearson
def test_pearson(x,y,p_threshold=0.01):
    stat, p = pearsonr(x, y)
    return (p <= p_threshold)

# Spearman
def test_spearmanr(x,y,p_threshold=0.01):
    stat, p = spearmanr(x, y)
    return (p <= p_threshold)

# Kendall
def test_kendalltau(x,y,p_threshold=0.01):
    stat, p = kendalltau(x, y)
    return (p <= p_threshold)

# Shapiro-Wilk
def test_normality_shapiro(x,y,p_threshold=0.01):
    stat, p1 = shapiro(x)
    stat, p2 = shapiro(y)
    return (p1 > p_threshold) and (p2 > p_threshold)

# Agostino
def test_normality_agostino(x,y,p_threshold=0.01):
    stat, p1 = normaltest(x)
    stat, p2 = shapiro(y)
    return (p1 > p_threshold) and (p2 > p_threshold)

def get_raw_data(name):
    # This function outputs six lists: raw, genes, genes_filtered, cancer, normal, normal_repeat
    # raw contains unfiltered gene expression data
    # genes contain names of the genes
    # genes_filtered contains names of the genes with significant samples
    # cancer contains gene expressions from cancer samples
    # normal contains gene expressions from non-cancer samples
    # normal_repeat scales normal list up to match size of corresponding cancer list
    genes = []
    raw = []
    genes_filtered = []
    cancer = []
    normal = []
    normal_repeat = []
    sample_ids = []
    with open(name) as f:
        content = f.readlines()
        print('Processing ', len(content), ' rows from', name )
        n = len(content) - 1
        for i in range(n):
            r = content[i]
            if i==0 :
                l = r.split('\t')
                # Get the sample ids
                sample_ids = [int(item.split('-')[3][0:2]) for item in l if item.startswith('TCGA')]
                h = len(l)
                c_count = h

            l = r.split('\t')
            assert(c_count == len(l))
            if i > 1:
                raw_temp = [float(i) for i in l[1:]]
                assert len(raw_temp) == len(sample_ids)
                raw.append(raw_temp)
                genes.append(l[0].split('|')[0])
                cancer_temp = []
                normal_temp = []
                for (code, value) in zip(sample_ids, raw_temp):
                    # We score cancer samples based on severity of cancer
                    if code in [1, 3, 5]:
                        # These codes are for primary cancer samples -> Score: 1
                        cancer_temp.append(value)
                    elif code in [2, 4]:
                        # These codes are for recurrent cancer samples -> Score: 2
                        cancer_temp += [value] * 2
                    elif code in [6, 7]:
                        # These codes are for metastatic samples -> Score: 4
                        cancer_temp += [value] * 4
                    elif code in [10, 11, 12, 13, 14]:
                        # These codes are for non-cancer samples
                        normal_temp.append(value)                
                    else:
                        print('Unexpected id ', code, ' occured for file ', name)

                if len(normal_temp) > 5 and len(cancer_temp) > 5: #only include if data is not too small
                    mult = int(len(cancer_temp)/len(normal_temp))
                    remainder = len(cancer_temp)%len(normal_temp)
                    #we scale non-cancer values up without disturbing the distribution
                    normal_rep_temp = normal_temp*mult
                    normal_rep_temp += normal_temp[0:remainder]

                    cancer.append(cancer_temp)
                    normal.append(normal_temp)
                    normal_repeat.append(normal_rep_temp)
                    genes_filtered.append(l[0].split('|')[0])

            if i%10000 == 0:
                print('Fetched ',i,' genes')
    return raw, genes, genes_filtered, cancer, normal, normal_repeat

def get_raw_data_from_files(input_path):
    #Data should be located in input_path subdirectory
    files = [item for item in os.listdir(input_path) if 'RSEM_genes_normalized' in item]
    if files == []:
        raise Exception("No data files found. Please download from firebrowser.")

    cancer_types = [file.split('.')[0] for file in files]
    data = []
    for file in files:
        data += [get_raw_data(os.path.join(input_path, file))]
    return data, cancer_types


def get_clusters(input_path, cancer_types, data):
    print('We are calculating clusters based on rbf kernel.')
    print('Please note that generating clusters will take long time and huge memory because of large number of genes involved')
    assert len(cancer_types) == len(data)
    genes = list(data[0][1])
    # At first, we will concatenate data for all cancer types
    data_concatenated = [[] for i in range(len(data[0][0]))]
    n = len(data_concatenated)
    for i in range(len(data)):
        assert genes == data[i][1]
        for j in range(len(data_concatenated)):
            data_concatenated[j].extend(list(data[i][0][j]))
        print('Finished appending ',len(data[i][0][0]),' rows for', cancer_types[i], ' data')
    
    
    # Given the concatenated data for all genes, we apply the rbf kernel(gaussian similarity matrix) on it
    # We take the laplacian matrix of the similarity matrix and then calculate eigenvalues and eigenvectors
    # Then we fit the eigenvectors for clustering
    # Because of the large size of the concatenated dataset, we do not calculate rbf kernel applied matrix seperately
    # We do as much processing as we can in the loops itself
    # The following two loops will take significant time and memory for large number of genes (10,000 onwards)
    data_np = np.asarray(data_concatenated)
    D = np.zeros((n,n))
    kernel = 1.0/len(data_concatenated[0]) #1/number of samples
    for i in range(n):
        if i%300 == 0 or i == n - 1:
            print('Done with calculating distance for ', i, ' genes out of ', n,' genes')
        for j in range(i):
            # We calculate kernel applied distance similarity at once
            D[j,i] = D[i,j] = -np.exp(-np.sum(np.square(data_np[i]-data_np[j]))*kernel)

    # We convert the matrix into Laplacian
    row_sum = np.sum(D, axis=1)
    for i in range(n):
        D[i,i] -= row_sum[i]

    # Now we calculate eigenvalues and eigenvectors from the Laplacian to cluster
    n_clusters = int(len(data[0][1])/100) # Dynamic cluster size based on number of genes
    clustering_successful = False
    try:
        e_val, e_vec = np.linalg.eig(D)
        ind = np.argsort(e_val)[:n_clusters]
        eigen_vectors = e_vec[:, ind]
    
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(eigen_vectors)
        labels = kmeans.labels_
    
        with open(os.path.join(input_path, 'gene_clustering_labels.csv'), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['gene', 'label'])
            for i in range(len(data[0][1])):
                writer.writerow({'gene': data[0][1][i], 'label': labels[i]})
        clustering_successful = True
    except:
        print('Some error occured during clustering. We will not be using clustering.')
        print('Please check the values using pdb.')
        labels = []
    return labels, clustering_successful

def detect_differential_expression(input_path, cancer_type, genes, cancer_data, normal_data, normal_repeat_data):
    print('Started detecting differential expression for cancer type: ', cancer_type)
    if(len(cancer_data) == 0 or len(normal_data) == 0):
        return []
    # We need the following to be of same length, requirement for some of the tests
    assert len(cancer_data) == len(normal_repeat_data)
    cancer_genes = [] # List of genes with differential expression
    for i in range(len(cancer_data)):
        try:
            is_gaussian = test_normality_shapiro(cancer_data[i],normal_data[i]) or test_normality_agostino(cancer_data[i],normal_data[i]) 
        except:
            is_gaussian = False
        if is_gaussian:
            # This is the gaussian case
            try:
                res_t = test_t(cancer_data[i],normal_data[i],p_threshold=0.1)
            except:
                res_t = False
            try:
                res_pearson = test_pearson(cancer_data[i],normal_repeat_data[i],p_threshold=0.1)
            except:
                res_pearson = False
            try:
                res_anova = test_anova(cancer_data[i],normal_data[i],p_threshold=0.1)
            except:
                res_anova = False
            if res_t or res_pearson or res_anova:
                cancer_genes.append(genes[i])
        else:
            # This is non-gaussian case
            try:
                res_spearmann = test_spearmanr(cancer_data[i],normal_repeat_data[i])
            except:
                res_spearmann = False
            try:
                res_mannwhitneyu = test_mannwhitneyu(cancer_data[i],normal_data[i])
            except:
                res_mannwhitneyu = False
            try:
                res_kendalltau = test_kendalltau(cancer_data[i],normal_repeat_data[i])
            except:
                res_kendalltau = False
            try:
                res_kruskal = test_kruskal(cancer_data[i],normal_data[i])
            except:
                res_kruskal = False
            if res_spearmann or res_mannwhitneyu or res_kendalltau or res_kruskal:
                cancer_genes.append(genes[i])

    print('Category ', cancer_type, ': Total number of genes ', len(cancer_data), ' differential expression detected in ', len(cancer_genes), ' genes')
    f = open(os.path.join(input_path, cancer_type + '_differential_gene_list.txt'), 'w')
    f.writelines(",".join(cancer_genes))
    f.close()
    return cancer_genes

def generate_differential_expression_gene_list(input_path, cancer_types, data):
    assert len(cancer_types) == len(data)
    gene_list = []
    for i in range(len(cancer_types)):
        # Remember: each data entry should have six entries: raw, genes, genes_filtered, cancer, normal, normal_repeat
        gene_list.append(detect_differential_expression(input_path, cancer_types[i], data[i][2], data[i][3], data[i][4], data[i][5]))
    return gene_list

def manual():
    print('Manual')
    print('=============================================================================================================')
    print('(1)\tEnsure input_file subdirectory is valid')
    print('If there are precalculated differential gene list, please put them under input path')
    print('File name format: <cancer category>_differential_gene_list.txt')
    print('(2)\tPut raw normlalized gene expression files in input path')
    print('File name format: <cancer category>.<arbitrary string>RSEM_genes_normalized<arbitrary string>')
    print('Ensure that input file exists')
    print('(3)\tGene names in relation with cancers will be saved as output path/cancer_gene_prediction.csv')
    print('Rows in cancer_gene_prediction.csv denote query genes')
    print('Columns in cancer_gene_prediction.csv denote cancer category')
    print('')
    print('Manadatory arguments:')
    print('(1)\tUse --input-file=x, where ''x'' denotes input gene lists. It should have one gene name in one line')
    print('(2)\tUse --output-file=x.csv, where ''x'' denotes output file name. It has to be a csv')
    print('')
    print('Optional arguments:')
    print('(1)\tUse --use-clustering argument to specify if KMean clustering is required')
    print('Note that corresponding label containing file should be saved as gene_clustering_labels.csv in input directory')
    print('If there are no labels, clusters will be recalculated')
    print('')
    print('(2)\tUse --force-analyze to reanalyze gene expressions even though the enriched data is present in subdirectory')
    print('')
    print('(3)\tUse --help to print this manual')


if __name__ == '__main__':
    # Check if it is --help
    call_for_manual = len(sys.argv) > 1 and '--help' in sys.argv[1:]
    if call_for_manual:
        manual()
        exit()

    if len(sys.argv) < 3 or len(sys.argv) > 5 or (len(sys.argv) > 1 and len([item for item in sys.argv[1:] \
           if not (item in ['--help', '--force-analyze', '--use-clustering'] or \
                   item.startswith('--input-file') or item.startswith('--output-file'))]) > 0):
        print('Incorrect option provided. Please use --help option for the manual.')
        exit()

    force_calculate = len(sys.argv) > 1 and '--force-analyze' in sys.argv[1:]

    # Check input directory sanctity
    input_file = [item for item in sys.argv[1:] if item.startswith('--input-file')][0].split('--input-file=', 1).pop()
    input_file_path = os.path.join(os.getcwd(), input_file)
    if not os.path.exists(input_file_path):
        print('Provide valid input file')
        exit()
    input_path = os.path.dirname(os.path.abspath(input_file_path))
    print('Input directory specified ', input_path)

    # Check output directory sanctity
    output_file = [item for item in sys.argv[1:] if item.startswith('--output-file')][0].split('--output-file=', 1).pop()
    output_file_path = os.path.join(os.getcwd(), output_file)
    output_path = os.path.dirname(os.path.abspath(output_file_path))
    if not os.path.isdir(output_path):
        print('Provide valid output directory')
        exit()
    elif not output_file.endswith('.csv'):
        print('Provide valid output file name. It should be ending with .csv')
        exit()
    print('Output directory specified ', output_path)

    # Now, we check if we need to cluster genes
    use_clustering = len(sys.argv) > 1 and '--use-clustering' in sys.argv[1:]
    data_fetched = False

    if use_clustering:
        try:
            print('We will be using clustering information')
            # Check if the labels are already present in enriched_data directory
            label_exists = len([item for item in os.listdir(input_path) if 'gene_clustering_labels.csv' in item]) > 0
            if label_exists:
                print('Found cluster labels')
                with open(os.path.join(input_path, 'gene_clustering_labels.csv'), mode='r') as g:
                    reader = csv.DictReader(g)
                    line_count = 0
                    labels = []
                    for row in reader:
                        labels.append(int(row['label']))
            else:
                data, cancer_types = get_raw_data_from_files(input_path)
                data_fetched = True
                labels, clustering_successful = get_clusters(input_path, cancer_types, data)
                use_clustering = clustering_successful
        except:
            use_clustering = False

    # Now we check if some clusters do not have excessive number of genes
    # Our number of clusters is number of genes divided by 100
    # So we check if number of genes is less than 3*100 = 300.
    # If not, we have too many genes in the cluster, making clusters redundant
    if use_clustering:
        counter=collections.Counter(labels)
        for frequency in counter.values():
            if frequency > 300:
                use_clustering = False

    # We create a dictionary from gene names and labels
    if use_clustering:
        if not data_fetched:
            data, cancer_types = get_raw_data_from_files(input_path)
        if(len(data[0][1]) == len(labels)):
            dictionary = dict(zip(data[0][1], labels))
        else: # Length did not match, something is wrong. Do not use clustering.
            use_clustering = False

    #At first, we check if we already have precomputed genelists with differential expression
    differential_gene_files = [item for item in os.listdir(input_path) if 'differential_gene_list' in item]
    if len(differential_gene_files) > 0 and not force_calculate:
        print('Differential gene files found in ',input_path ,' directory. These will be used.')
        differential_gene_list = []
        cancer_types = []
        for i in range(len(differential_gene_files)):
            f = open(os.path.join(input_path, differential_gene_files[i]), 'r')
            content = f.readlines()
            if(len(content) > 0):
                differential_gene_list.append(content[0].split(','))
                cancer_types.append(differential_gene_files[i].split('_')[0])
            f.close()
    else:
        if not data_fetched:
            data, cancer_types = get_raw_data_from_files(input_path)
        differential_gene_list = generate_differential_expression_gene_list(input_path, cancer_types, data)
        data_fetched = True

    # We scan through the input genelist
    with open(input_file_path) as f:
        genelist = [ item.strip('\n ') for item in f.readlines()]
        f.close()
    
    with open(os.path.join(output_file_path), 'w', newline='') as f:
        # This block will output into a csv 
        # Top header entries contain cancer categories, e.g. BRCA, GBMLGG, KIPAN and so on.
        # Cells that are at intersection of one query gene (row) and one cancer category (column)
        # denote the genes that have differential gene expressions in cancer vs non-cancer samples
        # in the corresponding cancer category and also have names starting with corresponding query gene name.
        writer = csv.writer(f)
        writer.writerow(['gene_names'] + cancer_types)
        flattened_differential_gene_list = list(set([item for sublist in differential_gene_list for item in sublist]))
        print('')
        for gene in genelist:
            cancerous = False
            related_categories = []
            related_genes = []
            for i in range(len(cancer_types)):
                # Filter the list to find the genes starting with query
                filtered_list = [item for item in differential_gene_list[i] if item.startswith(gene)]
                if(len(filtered_list) > 0):
                    cancerous += 1
                    related_categories.append(cancer_types[i])
                    related_genes.append((";".join(filtered_list))) # Append gene names so that it fits in one cell in csv
                else:
                    related_genes.append('') # An empty entry in csv
            writer.writerow([gene] + related_genes)
            if cancerous:
                print('Gene ',gene, ' is found to be cancerous in ', cancerous, ' categories of cancer')
            else:
                print('Gene ',gene, ' is not found to be cancerous')
                if use_clustering: # If we use clustering, we check for similar genes which have differential expressions
                    label = dictionary.get(gene, default=-1)
                    genes_in_same_cluster = list(dictionary.keys())[list(dictionary.values()).index(label)]
                    # We check if any of the genes in the same cluster appears in differential_gene_list
                    gene_intersection = list(set(genes_in_same_cluster) and set(flattened_differential_gene_list))
                if use_clustering and len(gene_intersection) > 0:
                    print('We found ', len(gene_intersection), ' genes in the same cluster having differential expression')
                    print('It might be useful to check the functions of gene ', gene)
        print('')
        print('Further details about gene names in relation with cancers is given in enriched_data/cancer_gene_prediction.csv')
        print('Rows indicate individual genes, columns indicate cancer category')
        print('')
        print('Use --help option if needed (e.g. python Predict_Diff_Expression.py --help)')
        exit()
