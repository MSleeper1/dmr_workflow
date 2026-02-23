import marimo

__generated_with = "0.2.13"
app = marimo.App()


@app.cell
def _():
    # notes

    # # divide index 0 by index 1 for first row to find proportion of sequences methylated
    # (cancer_beta[0:1, 0] / cancer_beta[0:1, 1])[0]

    # # create a dict of outputs from prop_meth function
    # proportion_methylated = {
    #     'cancer': prop_meth(cancer_beta_df),
    #     'control': prop_meth(control_beta_df),
    #     'cancerM': prop_meth(cancerM_beta_df)
    #     }

    # # create a column in cancer_beta_df that is the proportion of methylated sequences out of total sequences
    # cancer_beta_df['prop_meth'] = cancer_beta_df['num_meth'].div(cancer_beta_df['total_reads'])
    return


@app.cell
def _():
    import numpy as np
    import pandas as pd

    cancer = "/Users/meghansleeper/Desktop/farm-files/data/tissue-samples/271/merged/wgbstools-out/cancer-SRX381569-merged.beta"
    cancerM = "/Users/meghansleeper/Desktop/farm-files/data/tissue-samples/271/merged/wgbstools-out/cancerM-SRX381585-merged.beta"
    control = "/Users/meghansleeper/Desktop/farm-files/data/tissue-samples/271/merged/wgbstools-out/control-SRX381553-merged.beta"

    cancer_beta = np.fromfile(cancer, dtype=np.uint8).reshape((-1, 2))
    cancerM_beta = np.fromfile(cancerM, dtype=np.uint8).reshape((-1, 2))
    control_beta = np.fromfile(control, dtype=np.uint8).reshape((-1, 2))

    cancer_beta_df = pd.DataFrame(cancer_beta, columns=['num_meth', 'total_reads'])
    cancerM_beta_df = pd.DataFrame(cancerM_beta, columns=['num_meth', 'total_reads'])
    control_beta_df = pd.DataFrame(control_beta, columns=['num_meth', 'total_reads'])

    return (
        cancer,
        cancerM,
        cancerM_beta,
        cancerM_beta_df,
        cancer_beta,
        cancer_beta_df,
        control,
        control_beta,
        control_beta_df,
        np,
        pd,
    )


@app.cell
def _(np):
    def calc_beta_info(beta):
        length = len(beta)
        mean_coverage = np.mean(beta[:, 1])
        median_coverage = np.median(beta[:, 1])
        std_coverage = np.std(beta[:, 1])
        quantile_1 = np.quantile(beta[:, 1], 0.25)
        quantile_2 = np.quantile(beta[:, 1], 0.50)
        quantile_3 = np.quantile(beta[:, 1], 0.75)
        max_coverage = np.max(beta[:, 1])
        min_coverage = np.min(beta[:, 1])

        print("BETA SUMMARY: \n Rows: {0} (each is a CpG site)".format(length))
        print(" Columns: [# of methylated sequences, # of sequences total] \n", beta)
        print(" \n STATS FOR SEQ COVERAGE BY CPG SITE: \n   Mean:", mean_coverage)
        print("   Median:", median_coverage)
        print("   Standard deviation:", std_coverage)
        print("   25th, 50th, and 75th percentiles: {0}, {1}, {2}.".format(quantile_1, quantile_2, quantile_3))
        print("   Maximum: {} (limited to 255 by unit8 format of beta file)".format(max_coverage))
        print("   Minimum:", min_coverage)
        print("\n")

        return length, mean_coverage, median_coverage, std_coverage, quantile_1, quantile_2, quantile_3, max_coverage, min_coverage
        
    return calc_beta_info,


@app.cell
def _(calc_beta_info, cancerM_beta, cancer_beta, control_beta, pd):
    file_info = {
        'cancer': calc_beta_info(cancer_beta),
        'control': calc_beta_info(control_beta),
        'cancerM': calc_beta_info(cancerM_beta)
        }

    beta_info_df = pd.DataFrame.from_dict(file_info, orient='index', 
                                columns=['length', 'mean_cov', 'median_cov', 'std_cov',
                                        'quant_1', 'quant_2', 'quant_3', 'max_cov', 'min_cov'])
    return beta_info_df, file_info


@app.cell
def _():
    # divide number of sequences methylated by total sequences to find proportion of sequences methylated
    def prop_meth(beta_df):
        beta_df['prop_meth'] = beta_df['num_meth'].div(beta_df['total_reads'])
        return beta_df

    return prop_meth,


@app.cell
def _(cancerM_beta_df, cancer_beta_df, control_beta_df, prop_meth):
    cancer_beta_df = prop_meth(cancer_beta_df)
    control_beta_df = prop_meth(control_beta_df)
    cancerM_beta_df = prop_meth(cancerM_beta_df)
    return cancerM_beta_df, cancer_beta_df, control_beta_df


@app.cell
def _(cancerM_beta_df, cancer_beta_df, control_beta_df, pd):
    # combine the three dataframes into one dataframe
    beta_df = pd.concat([cancer_beta_df, control_beta_df, cancerM_beta_df], axis=1, keys=['cancer', 'control', 'cancerM'])

    beta_df.head()

    return beta_df,


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

