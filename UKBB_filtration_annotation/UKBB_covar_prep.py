import pandas as pd
import os


def with_proxy(PD_covar, proxy_covar, control_covar):
    # STEP 1: Prepare covar dataframe
    covar_df = pd.concat([PD_covar, control_covar.iloc[1:], proxy_covar.iloc[1:]], ignore_index=True)
    covar_df.insert(1, "FID", covar_df.iloc[:, 0])
    covar_df.insert(2, "IID", covar_df.iloc[:, 0])
    covar_df.drop(covar_df.columns[0], axis=1, inplace=True)
    # Females encoded as 0 and males as 1 in the original file from UKBB
    # https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=9
    covar_df['Sex'] = covar_df['Sex'].map({0: 2, 1: 1}).fillna(0)

    # STEP 2: Create the samples dataframe with FID and IID columns
    samples_df = covar_df.iloc[:, [0, 1]].copy()
    # STEP 3: Create sex dataframe with FID, IID, and Sex columns
    sex_df = covar_df.iloc[:, [0, 1, 6]].copy()
    # STEP 4: Create the phenotype dataframe with FID, IID and Status columns
    pheno_df = covar_df.iloc[:, [0, 1]].copy()
    pheno_df['Status'] = pheno_df['FID'].apply(
        lambda x: 2 if x in PD_covar['ID'].values or x in proxy_covar['ID'].values else 1)

    covar_df['Status'] = pheno_df['Status']

    # check for NA values
    samples_df.dropna(axis=0, inplace=True)
    sex_df.dropna(axis=0, inplace=True)
    pheno_df.dropna(axis=0, inplace=True)
    covar_df.dropna(axis=0, inplace=True)
    # write the files
    samples_df.to_csv(os.path.join(output_path_proxy, 'UKBB.samples.list'),
                      sep='\t', index=False)
    sex_df.to_csv(os.path.join(output_path_proxy, 'sex_' + 'UKBB' + '.txt'),
                  sep='\t', index=False)
    pheno_df.to_csv(os.path.join(output_path_proxy, 'pheno_' + 'UKBB' + '.txt'),
                    sep='\t', index=False)
    covar_df.to_csv(os.path.join(output_path_proxy, 'covar_' + 'UKBB' + '.txt'),
                    sep='\t', index=False)


def without_proxy(PD_covar, control_covar):
    # STEP 1: Prepare covar dataframe
    covar_df = pd.concat([PD_covar, control_covar.iloc[1:]], ignore_index=True)
    covar_df.insert(1, "FID", covar_df.iloc[:, 0])
    covar_df.insert(2, "IID", covar_df.iloc[:, 0])
    covar_df["FID"] = covar_df.iloc[:, 0].astype(int)
    covar_df["IID"] = covar_df.iloc[:, 0].astype(int)
    covar_df.drop(covar_df.columns[0], axis=1, inplace=True)
    # Females encoded as 0 and males as 1 in the original file from UKBB
    # https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=9
    covar_df['Sex'] = covar_df['Sex'].map({0: 2, 1: 1}).fillna(0)

    # STEP 2: Create the samples dataframe with FID and IID columns
    samples_df = covar_df.iloc[:, [0, 1]].copy()
    # STEP 3: Create sex dataframe with FID, IID, and Sex columns
    sex_df = covar_df.iloc[:, [0, 1, 6]].copy()
    # STEP 4: Create the phenotype dataframe with FID, IID and Status columns
    pheno_df = covar_df.iloc[:, [0, 1]].copy()
    pheno_df['Status'] = pheno_df['FID'].apply(
        lambda x: 2 if x in PD_covar['ID'].values else 1)

    covar_df['Status'] = pheno_df['Status']

    # check for NA values
    samples_df.dropna(axis=0, inplace=True)
    sex_df.dropna(axis=0, inplace=True)
    pheno_df.dropna(axis=0, inplace=True)
    covar_df.dropna(axis=0, inplace=True)

    # write the files
    samples_df.to_csv(os.path.join(output_path_no_proxy, 'UKBB.samples.list'),
                      sep='\t', index=False)
    sex_df.to_csv(os.path.join(output_path_no_proxy, 'sex_' + 'UKBB' + '.txt'),
                  sep='\t', index=False)
    pheno_df.to_csv(os.path.join(output_path_no_proxy, 'pheno_' + 'UKBB' + '.txt'),
                    sep='\t', index=False)
    covar_df.to_csv(os.path.join(output_path_no_proxy, 'covar_' + 'UKBB' + '.txt'),
                    sep='\t', index=False)


file_path = "/Users/cemparlar/Desktop/Gan-Or/*OCT_LRRK2_UKBB/new covar/"
output_path_proxy = "/Users/cemparlar/Desktop/Gan-Or/*OCT_LRRK2_UKBB/covar_UKBB_proxy/"
output_path_no_proxy = "/Users/cemparlar/Desktop/Gan-Or/*OCT_LRRK2_UKBB/covar_UKBB_no_proxy/"

# create the directories for the output files if the folders don't exist

if not os.path.exists(output_path_proxy):
    os.mkdir(output_path_proxy)

if not os.path.exists(output_path_no_proxy):
    os.mkdir(output_path_no_proxy)

# Load the initial covar files
PD_covar = pd.read_csv(file_path + "ukbb_PD_covar.txt", delimiter=',')

# With Proxy
control_covar = pd.read_csv(file_path + "ukbb_PD_control_covar.txt", delimiter=',')
proxy_covar = pd.read_csv(file_path + "ukbb_PD_proxy_covar.txt", delimiter=',')

# Without Proxy
num_samples = 65000
header = control_covar.iloc[0]
control_covar_no_proxy = control_covar.sample(n=num_samples, random_state=42)  # Set random_state for reproducibility
control_covar_no_proxy = pd.concat([header.to_frame().T, control_covar_no_proxy])

print(PD_covar.shape, control_covar.shape, proxy_covar.shape, control_covar_no_proxy.shape)

# WRITE WITH PROXY
with_proxy(PD_covar, proxy_covar, control_covar)

# WRITE WITHOUT PROXY
without_proxy(PD_covar, control_covar_no_proxy)







