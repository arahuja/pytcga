import os
import tarfile
import pandas as pd

from pytcga.tcga_requests import tcga_request
from pytcga.tcga_clinical import load_clinical_data

def prefetch_rnaseq_data(disease_code,
                        wait_time=30,
                        cache=True,):

    archive_path = tcga_request(disease=disease_code,
                               level='3',
                               center='7',
                               platformType='RNASeqV2',
                               platform='IlluminaHiSeq_RNASeqV2',
                               wait_time=30)

    return archive_path


def load_rnaseq_data(disease_code,
                     with_clinical=False,
                     wait_time=30):
    # Fetch RNA data
    archive_path = prefetch_rnaseq_data(disease_code)

    # Unpack tar file
    archive = tarfile.open(archive_path)
    result_dir = os.path.join(os.path.dirname(archive_path), disease_code, 'gene_expression')

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    gene_quantification_files = set(el for el in archive if 'genes.normalized_results' in el.name)
    if len(gene_quantification_files.difference(os.listdir(result_dir))) > 0:
        archive.extractall(members=gene_quantification_files, path=result_dir)
        archive.extract('FILE_SAMPLE_MAP.txt', path=result_dir)

    gene_quantification_files = [el.name for el in gene_quantification_files]

    # Load map from samples to RNA files
    rna_file_sample_map = pd.read_csv(os.path.join(result_dir, 'FILE_SAMPLE_MAP.txt'), sep='\t')
    rna_file_sample_map_id_split = rna_file_sample_map['barcode(s)'].str.rsplit('-', n=4, expand=True)
    rna_file_sample_map_id_split.columns = ['TCGA_ID', 'SampleID', 'PortionID', 'PlateID', 'CenterID']

    rna_file_sample_map = rna_file_sample_map.join(rna_file_sample_map_id_split)

    gene_filter = rna_file_sample_map['filename'].str.contains('genes.normalized_results')
    gene_rna_file_sample_map = rna_file_sample_map[gene_filter]

    rna_dfs = []
    for (f, sample) in list(zip(gene_rna_file_sample_map['filename'], gene_rna_file_sample_map['TCGA_ID'])):
        sample_rna_df = pd.read_csv(os.path.join(result_dir, f), sep='\t')
        sample_rna_df['TCGA_ID'] = sample

        sample_rna_df['gene_name'] = sample_rna_df.gene_id.str.split('|').str.get(0)
        rna_dfs.append(sample_rna_df)

    rna_df = pd.concat(rna_dfs)

    if with_clinical:
        patient_data_df = load_clinical_data(disease_code)
        merged = rna_df.merge(patient_data_df,
                              how='outer',
                              left_on='TCGA_ID', 
                              right_on='bcr_patient_barcode')

        return merged
    else: 
        return rna_df