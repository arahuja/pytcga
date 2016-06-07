import os
import logging
import requests
from bs4 import BeautifulSoup
import pandas as pd

from .tcga_requests import cache_data_dir
from .tcga_utils import load_tcga_tabfile
from .clinical_data_dictionary import clinical_data_dictionary

TCGA_CLINICAL_URL = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/{}/bcr/biotab/clin/"

PATIENT_DATA_FILE_CODE = 'clinical_patient'
def request_clinical_data(disease_code,
                  cache=True,
                  block_size=1024):
    """Downloads TCGA public clinical data from the TCGA FTP site

    Parameters
    ----------
    disease_code : str
        TCGA disease type, i.e. 'LUAD', 'BLCA', 'BRCA' etc.
    cache : bool, optional
        Whether to cache the results of the request
    block_size : int, optional
        Block size for file downloads

    Returns
    -------
    patient_data_path : str
        Path to TCGA patient data file after downloading
    """
    # Create directory to save clinical data
    disease_code_dir = os.path.join(cache_data_dir(), disease_code)

    if cache and os.path.exists(disease_code_dir):
        patient_data_file = [f for f in os.listdir(disease_code_dir) if PATIENT_DATA_FILE_CODE in f]

        if len(patient_data_file) == 1:
            return os.path.join(disease_code_dir, patient_data_file[0])

    if not os.path.exists(disease_code_dir):
        os.makedirs(disease_code_dir)

    clinical_data_directory = TCGA_CLINICAL_URL.format(disease_code.lower())
    r = requests.get(clinical_data_directory)
    soup = BeautifulSoup(r.content, "html.parser")

    # Retrieve list of files and filter to txt files
    file_links = [link.get('href')
                    for link in soup.find_all('a')]
    clinical_files = [link for link in file_links if link.endswith('.txt')]

    # Download all clinical data files
    patient_data_path = None
    for clinical_file in clinical_files:
        output_file = os.path.join(disease_code_dir, clinical_file)
        logging.debug('Saving {} clinical data request to {}'.format(clinical_file, output_file))

        if PATIENT_DATA_FILE_CODE in output_file:
            patient_data_path = output_file

        with open(output_file, 'wb') as archive:
            archive_response = requests.get(clinical_data_directory + '/' + clinical_file, stream=True)

            for block in archive_response.iter_content(block_size):
                archive.write(block)

    return patient_data_path

def load_patient_data(disease_code, recode_columns=True):
    return load_clinical_data(disease_code, recode_columns)

def load_clinical_data(disease_code, recode_columns=True):
    """Downloads and loads the TCGA clinical data into a Pandas dataframe

    Parameters
    ----------
    disease_code : str
        TCGA disease type, i.e. 'LUAD', 'BLCA', 'BRCA' etc.

    Returns
    -------
    patient_data_df : Dataframe
        Returns a Pandas dataframe with the patient data
    """
    patient_data_path = request_clinical_data(disease_code, cache=True)

    patient_data_df = load_tcga_tabfile(patient_data_path, skiprows=1)

    if recode_columns:
        for column in clinical_data_dictionary:
            if column in patient_data_df.columns:
                patient_data_df[column] = patient_data_df[column].map(
                        lambda val: clinical_data_dictionary[column].get(val, val)
                    )

    logging.info("Loaded {} rows of clinical data from {} patients".format(
            len(patient_data_df),
            patient_data_df['bcr_patient_barcode'].nunique()
        )
    )

    return patient_data_df

def find_clinical_files(search_tag, disease_code_dir):
    files = [f for f in os.listdir(disease_code_dir) if search_tag in f]
    return files

def _load_samples(disease_code, filter_vial=None):
    disease_code_dir = os.path.join(cache_data_dir(), disease_code)
    sample_files = find_clinical_files('_biospecimen_sample_', disease_code_dir)
    sample_df = pd.concat(
        [load_tcga_tabfile(os.path.join(disease_code_dir, f))
            for f in sample_files],
        copy=False)
    if filter_vial:
        sample_df = sample_df[sample_df.vial_number == filter_vial]
    return sample_df

def _load_analytes(disease_code):
    disease_code_dir = os.path.join(cache_data_dir(), disease_code)
    analyte_files = find_clinical_files('_biospecimen_analyte_', disease_code_dir)
    analyte_df = pd.concat(
        [load_tcga_tabfile(os.path.join(disease_code_dir, f))
            for f in analyte_files],
        copy=False)
    return analyte_df

def load_treatments(disease_code):
    """Load the treatment entries for each patient

    Parameters
    ----------
    disease_code : str
        TCGA disease code

    Returns
    -------
    treatment_df : Pandas dataframe
        Dataframe of treatment entries for each patient
    """
    disease_code_dir = os.path.join(cache_data_dir(), disease_code)
    treatment_files = find_clinical_files('_clinical_drug', disease_code_dir)

    treatment_df = pd.concat(
        [load_tcga_tabfile(os.path.join(disease_code_dir, f), skiprows=1)
            for f in treatment_files],
        copy=False)
    return treatment_df

def load_patient_samples(disease_code, recode_columns=True, filter_vial=None):
    """Load the samples taken per patient"""
    patient_data = load_patient_data(disease_code, recode_columns)
    samples = _load_samples(disease_code, filter_vial)

    return patient_data.merge(samples, how='left')

def load_patient_analytes(disease_code, recode_columns=True):
    """Load the analytes per sample. Possible analytes include RNA or DNA"""
    patient_data = load_patient_data(disease_code, recode_columns)
    analytes = _load_analytes(disease_code)

    return patient_data.merge(analytes, how='left')

def load_sample_and_analytes(disease_code, filter_vial=None):

    samples = _load_samples(disease_code, filter_vial=filter_vial)
    analytes = _load_analytes(disease_code)

    return samples.merge(analytes)

def load_aliquots(disease_code, recode_columns=True):
    """Load the aliqouts taken per patient"""
    disease_code_dir = os.path.join(PYTCGA_BASE_DIRECTORY, disease_code)
    aliquot_files = find_clinical_files('_biospecimen_aliquot_', disease_code_dir)
    aliquot_df = pd.concat(
        [load_tcga_tabfile(os.path.join(disease_code_dir, f))
            for f in aliquot_files],
        copy=False)
    return aliquot_df
