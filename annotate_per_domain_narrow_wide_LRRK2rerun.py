import pandas as pd

# Required files
# {cohort_name}_FILTERED_IMPORTANT_GENES.bim
# {genes_list} or {domains_list}, make sure to have the correct headers
# {cohort_name}_CADD.txt, {cohort_name}_ENCODE.txt, {cohort_name}_LOF.txt, {cohort_name}_nonsyn.txt

# Cohort name
base_dir = '/Users/cemparlar/Desktop/Gan-Or/*JULY2024_yesnoG2019S_LRRK2/SETID_PREP'
cohort_name = 'ISRAEL'

# Read the input files into pandas dataframes
domains_df = pd.read_excel(f'{base_dir}/LRRK2_GRCh38.xlsx',
                           sheet_name="Sheet2")  # headers are: Chromosome, Gene, Domain, start, end
domains_df['Chromosome'] = domains_df['Chromosome'].astype(str)
domains_df['Chromosome'] = domains_df['Chromosome'].str.replace('chr', '')
# drop nan values in start and end
domains_df['narrow_start'] = domains_df['narrow_start'].astype(int)
domains_df['narrow_end'] = domains_df['narrow_end'].astype(int)
domains_df['wide_start'] = domains_df['wide_start'].astype(int)
domains_df['wide_end'] = domains_df['wide_end'].astype(int)
# remove all spaces from gene names
domains_df['Gene'] = domains_df['Gene'].str.strip().str.replace(' ', '')
domains_df['Domain'] = domains_df['Domain'].str.strip().str.replace(' ', '')

# Separate wide and narrow into different dataframes
domains_df_narrow = domains_df[['Chromosome', 'Gene', 'Domain', 'narrow_start', 'narrow_end']].rename(
    columns={'narrow_start': 'start', 'narrow_end': 'end'})
domains_df_wide = domains_df[['Chromosome', 'Gene', 'Domain', 'wide_start', 'wide_end']].rename(
    columns={'wide_start': 'start', 'wide_end': 'end'})

domains_df_narrow = domains_df_narrow.dropna(subset=['start', 'end'])
domains_df_wide = domains_df_wide.dropna(subset=['start', 'end'])

bim_variants_df = pd.read_csv(f'{base_dir}/{cohort_name}.bim',
                              sep='\t', header=None,
                              names=['Chr', 'variantID', 'morgan_position', 'position', 'A1', 'A2'])
bim_variants_df['Chr'] = bim_variants_df['Chr'].astype(str)


# functional snps list retrieved from ANNOVAR
# functional_highCADD_df = pd.read_csv(f'{base_dir}/{cohort_name}_functional_hjghCADD.txt', header=None, names=['variantID'])
# functional_snp_df = pd.read_csv(f'{base_dir}/{cohort_name}_functional_snp.txt', header=None, names=['variantID'])
ALL_df = bim_variants_df[['variantID']]
CADD_df = pd.read_csv(f'{base_dir}/{cohort_name}_CADD.txt', header=None, names=['variantID'])
encode_df = pd.read_csv(f'{base_dir}/{cohort_name}_ENCODE.txt', header=None, names=['variantID'])
LOF_df = pd.read_csv(f'{base_dir}/{cohort_name}_LOF.txt', header=None, names=['variantID'])
nonsyn_df = pd.read_csv(f'{base_dir}/{cohort_name}_nonsyn.txt', header=None, names=['variantID'])


list_of_dfs = {'CADD': CADD_df, 'EXONIC': encode_df, 'LOF': LOF_df, 'nonsyn': nonsyn_df}

# Iterate through each row in the DNM3 dataframe


def narrow_wide(domains_df, string_narrow_or_wide):

    setid_filename = f'{base_dir}/{cohort_name}.SETID'

    with open(setid_filename, 'a') as output_file:
        for row in domains_df.itertuples(index=False):
            # Extract relevant information from the DNM3 dataframe
            start, end = sorted([row.start, row.end])
            gene = row.Gene
            chromosome = row.Chromosome
            domain = row.Domain
            # Filter variants based on the interval defined by the current gene/domain row

            filtered_variants = bim_variants_df[
                (bim_variants_df['Chr'] == chromosome) & (bim_variants_df['position'] >= start) & (
                            bim_variants_df['position'] <= end)]

            # Use set operations for faster checks
            if cohort_name in ['AMP_PD', 'UKBB']:
                matching_variants = set(filtered_variants['variantID'])
                for df_key, df_value in list_of_dfs.items():
                    common_variants = matching_variants.intersection(df_value.variantID)
                    for variantID in common_variants:
                        output_file.write(
                            f'{gene}~{domain}~{string_narrow_or_wide}~{df_key}\t{variantID}\t"{gene}"\t"{domain}"\n')
            else:
                filtered_variants.loc[:, 'Custom_Column_1'] = filtered_variants.apply(lambda row: f"{row['Chr']}_{row['position']}_{row['A2']}_{row['A1']}", axis=1)
                matching_variants_df = filtered_variants[['Custom_Column_1', 'variantID']]
                for df_key, df_value in list_of_dfs.items():
                    common_variants = matching_variants_df[matching_variants_df['Custom_Column_1'].isin(df_value['variantID'])]
                    for index, row in common_variants.iterrows():
                        output_file.write(
                            f'{gene}~{domain}~{string_narrow_or_wide}~{df_key}\t{row["variantID"]}\t"{gene}"\t"{domain}"\n')

    output_file.close()


with open(f'{base_dir}/{cohort_name}.SETID', 'w'):
    pass  # This opens the file and immediately closes it, effectively clearing its contents

narrow_wide(domains_df_narrow, 'narrow')
narrow_wide(domains_df_wide, 'wide')

# Get stats on how many variants per SET were found in SETID

# Read the SETID file into a dataframe
setid_df = pd.read_csv(f'{base_dir}/{cohort_name}.SETID', sep='\t', header=None,
                       names=['SETID', 'variantID'])

print(setid_df['SETID'].value_counts())
