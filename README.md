## pytcga

Python library for accessing and processing public TCGA data


#### Example query

```python

import pytcga

# Downloading LUAD mutations
luad_mutations = \
    pytcga.request(
        disease='LUAD',
        level='2',
        center='BI',
        platformType='Somatic Mutations'
        platform='Automated Mutation Calling')

# Downloading LUAD RNASeq
luad_rna = \
    pytcga.request(
        disease='LUAD',
        level='3',
        center='7',
        platformType='RNASeqV2',
        platform='IlluminaHiSeq_RNASeqV2',
        wait_time=30)

```
