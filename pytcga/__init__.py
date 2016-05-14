from pytcga.tcga_requests import tcga_request, RequestError
from pytcga.tcga_clinical import load_clinical_data, load_patient_data, load_patient_samples, load_patient_analytes, load_treatments, load_sample_and_analytes, load_aliquots
from pytcga.tcga_mutations import load_mutation_data
from pytcga.tcga_rna import load_rnaseq_data
from pytcga.tcga_utils import load_studies
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
