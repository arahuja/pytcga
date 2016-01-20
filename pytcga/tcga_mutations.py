import os
import logging
import tarfile
import pandas as pd 

from pytcga.tcga_requests import tcga_request
from pytcga.tcga_clinical import load_clinical_data

def prefetch_mutation_data(disease_code,
                          wait_time=30,
                          cache=True,):

    archive_path = tcga_request(disease=disease_code,
                               level='2',
                               center='BI',
                               platformType='Somatic Mutations',
                               platform='Automated Mutation Calling',
                               wait_time=wait_time,
                               cache=cache)

    return archive_path

def load_mutation_data(disease_code,
                       with_clinical=False,
                       variant_type='all',
                       wait_time=30):

    archive_path = prefetch_mutation_data(disease_code,
                                      wait_time=wait_time,
                                      cache=True)

    # Unpack tar file
    archive = tarfile.open(archive_path)

    result_dir = os.path.join(os.path.dirname(archive_path), disease_code, 'mutations')
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)

    archive.extractall(path=result_dir)

    # Filter to MAF files
    maf_files = [f 
                  for f in os.listdir(result_dir) 
                  if f.endswith('.maf')]

    mutation_df = pd.concat([pd.read_csv(os.path.join(result_dir, maf_file),
                                    sep='\t', 
                                    na_values='[Not Available]') 
                    for maf_file in maf_files], copy=False)

    # Expand out the TCGA barcode to retrieve the TCGA ID
    tcga_info = mutation_df['Tumor_Sample_Barcode'].str.rsplit('-', n=4, expand=True)
    tcga_info.columns = ['TCGA_ID', 'SampleID', 'PortionID', 'PlateID', 'CenterID']

    mutations = mutation_df.join(tcga_info, how='left')

    if variant_type != 'all':
        if variant_type == 'indel':
            mutations = mutations[
                            (mutations['Variant_Type'] == 'INS') |
                            (mutations['Variant_Type'] == 'DEL')
                        ]          
        else:     
            mutations = mutations[mutations['Variant_Type'] == variant_type]

    logging.info("Loaded {} mutations for {} tumors from {} patients".format(
                    len(mutations),
                    mutations['Tumor_Sample_Barcode'].nunique(),
                    mutations['TCGA_ID'].nunique()
                )
    )

    if with_clinical:
        patient_data_df = load_clinical_data(disease_code)
        merged = mutations.merge(patient_data_df,
                            how='outer',
                            left_on='TCGA_ID', 
                            right_on='bcr_patient_barcode')

        logging.info("Patients: {}, Tumor Samples: {}, Mutations {}".format(
                    merged['bcr_patient_barcode'].nunique(),
                    merged['Tumor_Sample_Barcode'].nunique(),
                    len(merged[~merged['TCGA_ID'].isnull()])
                )
        )
        return merged
    else: 
        return mutations
