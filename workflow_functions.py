import pandas as pd

#### Python functions for snakemake pipeline ####

''' Read sample info table into pandas dataframe and confirm column names that are needed for the pipeline'''
def get_sample_info_df(input_tsv):
    sample_info = pd.read_table(input_tsv, dtype=str)
    # Verify column names
    if not {'ref', 'patient_id', 'srx_id', 'accession', 'layout', 'group', 'age', 'sex', 'ethnicity'}.issubset(sample_info.columns.values):
        raise KeyError("The sample info file must contain the following named columns: 'ref', 'patient_id', 'srx_id', 'accession', 'layout', 'group', 'age', 'sex' 'ethnicity'")
    sample_info = sample_info.set_index(["srx_id"], drop=False)
    return sample_info
