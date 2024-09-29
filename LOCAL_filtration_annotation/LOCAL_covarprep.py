import pandas as pd
import os
import pprint

df = pd.read_excel('/Users/cemparlar/Desktop/Gan-Or/*McGill Data/2023-04-25-all-info-by-ROUnum.xlsx')


output_path = '/Users/cemparlar/Desktop/Gan-Or/cohort_files'

groups = df.groupby('Samples_GOnum')


# Step 2-5: Define function to determine preferred sample
def get_preferred_sample(group):
    # Step 3: Check for duplicates with non-empty cells in 'PD_MIPs_submission' and 'Lys_MIPs_submission'
    non_empty_duplicates = group[(group['PD_MIPs_submission'].notnull()) & (group['Lys_MIPs_submission'].notnull())]
    if not non_empty_duplicates.empty:
        return non_empty_duplicates.iloc[0]  # Return the first duplicate with non-empty cells in both columns

    # Step 4: Check for duplicates with at least one non-empty cell
    non_empty_cells = group[(group['PD_MIPs_submission'].notnull()) | (group['Lys_MIPs_submission'].notnull())]
    if not non_empty_cells.empty:
        return non_empty_cells.iloc[0]  # Return the first duplicate with at least one non-empty cell

    # Step 5: Return any sample if none of the duplicates have any information in either column cell
    return group.iloc[0]


# Step 6: Apply the function to each group
preferred_samples = groups.apply(get_preferred_sample)

# Step 7: Concatenate preferred samples into a new DataFrame
df = preferred_samples.reset_index(drop=True)

# Remove unwanted New York cohort samples in the main database
excluded_values = ['GPD9{:02d}'.format(i) for i in range(0, 1000)] + ['GPD9{:03d}'.format(i) for i in range(1000, 10000)] + ['GPD998C', 'GPD999C']
df = df[~df['Sender_number'].isin(excluded_values)]

# Filter for PD or control in 'Disease' column for the whole database but not case-sensitive
df = df[df['Disease'].str.contains('PD|control', case=False)]
# Set up the cohort names
initial_cohort_names = df['Cohort'].unique()
cohort_names = {}

for cohort in initial_cohort_names:
    info = str(cohort).split('_')
    location = ''
    country = ''
    if len(info) == 1: # CORRECTING FROM OLD NAMING FORMAT (Quebec, France, CONTROLS), they all belong to FRENCH
        if info[0].lower() == 'controls':
            location = 'QUEBEC'
            country = 'CANADA'
        elif info[0].lower() == 'quebec':
            location = 'QUEBEC'
            country = 'CANADA'
        elif info[0].lower() == 'france':
            country = 'FRANCE'
        else:
            continue
    elif len(info) == 2:  # SENDER_COUNTRY
        country = info[1]
    elif len(info) == 3:  # SENDER_LOCATION_COUNTRY
        location = info[1]
        country = info[2]
    else:
        continue

    cohort_name = country.upper()
    if cohort_name == 'CANADA' and location.lower() in ['montreal', 'quebec city', 'quebec']:
        cohort_name = 'FRENCH'
    if cohort_name == 'FRANCE':
        cohort_name = 'FRENCH'
    if cohort_name in cohort_names:
        cohort_names[cohort_name].append(cohort)
    else:
        cohort_names[cohort_name] = [cohort]

pprint.pprint(cohort_names)
# DO QUALITY CONTROL NOW
for cohort_name in cohort_names:
    if cohort_name == "FRENCH":
        invalid_values = ['Gender mismatch', 'Heterozygocity outlier', 'M - non-euro ancestry - duplicate',
                          'Non-European', 'M - gender mismatch - duplicate', 'M - gender mismatch - non-euro ancestry',
                          'Related', 'yes - gender mismatch']
    else:
        invalid_values = ['Gender mismatch', 'Heterozygocity outlier', 'M - gender mismatch - duplicate',
                          'M - gender mismatch - non-euro ancestry', 'Related', 'yes - gender mismatch']
    df = df[~df['Patients_Exclude_from_analysis'].isin(invalid_values)]

# create the sample lists
sample_dict = {}

for key, value in cohort_names.items():
    sample_list = []
    for cohort in value:
        # from the main database selects the samples that belong to the cohort and have a sample number
        samples = df[(df['Cohort'] == cohort) & df['Snum'].notnull()]['Snum'].tolist()
        sample_list += samples
    # creates a dictionary with cohort names as keys and sample lists as values
    sample_dict[key] = sample_list
    # prints the cohort and the cohort sample size for each cohort
    print(key, len(sample_list))

# create the cohort folders and output the sample lists
try:
    os.mkdir(output_path)
except FileExistsError:
    pass
empty_cohorts = []
full_cohorts = []
total_valid_samples = []
for cohort in cohort_names:

    samples = pd.DataFrame(sample_dict[cohort], columns=['Sample ID'])
    # Create sex_$cohort file (FID, IID, Sex)
    sex_df = pd.DataFrame({'FID': samples['Sample ID'], 'IID': samples['Sample ID'], 'Sex': None})
    for i, row in samples.iterrows():
        sample_id = row['Sample ID']
        gender_series = df.loc[(df['Snum'] == sample_id), 'Gender']
        gender = gender_series.values[0]
        if not isinstance(gender, str):
            sex = 0
        else:
            sex: int = 1 if gender == 'M' else 2 if gender == 'F' else 0
        sex_df.loc[i, 'Sex'] = sex

    # Create pheno_$cohort file (FID, IID, Status)
    pheno_df = pd.DataFrame({'FID': samples['Sample ID'], 'IID': samples['Sample ID'], 'Status': None})
    for i, row in samples.iterrows():
        sample_id = row['Sample ID']
        status_series = df.loc[(df['Snum'] == sample_id), 'Disease']
        if status_series.empty:
            pheno_df.drop(i, inplace=True)
            continue
        status = status_series.values[0]
        status: int = 1 if status.upper() == 'CONTROL' else 2 if status.upper() == 'PD' else -9
        if status == -9:
            pheno_df.drop(i, inplace=True)
            continue
        pheno_df.loc[i, 'Status'] = status

    # Create covar_$cohort file (FID, IID, Status, Sex, Age, Ethn)
    covar_df = pd.DataFrame({'FID': samples['Sample ID'], 'IID': samples['Sample ID'],
                             'Status': None, 'Sex': None, 'Age': None, 'Ethn': '1'})
    for i, row in samples.iterrows():
        sample_id = row['Sample ID']
        # Add status
        status_series = df.loc[(df['Snum'] == sample_id), 'Disease']
        if status_series.empty:
            covar_df.drop(i, inplace=True)
            continue
        status = status_series.values[0]
        status: int = 1 if status.upper() == 'CONTROL' else 2 if status.upper() == 'PD' else -9
        if status == -9:
            covar_df.drop(i, inplace=True)
            continue
        covar_df.loc[i, 'Status'] = status

        # Add sex
        gender_series = df.loc[(df['Snum'] == sample_id), 'Gender']
        gender = gender_series.values[0]
        if not isinstance(gender, str):
            sex = 0
        else:
            sex: int = 1 if gender == 'M' else 2 if gender == 'F' else 0
        covar_df.loc[i, 'Sex'] = sex

        # Add age (THIS WORKS JUST FOR PD)
        age_series = df.loc[(df['Snum'] == sample_id), 'Age_at_onset_PD']
        if age_series.isna().any():
            age_series = df.loc[(df['Snum'] == sample_id), 'Age_at_diagnosis_PD']
            if age_series.isna().any():
                # Added step for including controls in the covar file
                age_series = df.loc[(df['Snum'] == sample_id), 'Age_at_Sample']
                if age_series.isna().any():
                    covar_df.drop(i, inplace=True) # drop the sample if no age is found
                    continue
        age = age_series.values[0]
        covar_df.loc[i, 'Age'] = int(age)

        if cohort in ['ISRAEL', 'USA']:
            # Add ethnicity
            eth_series = df.loc[(df['Snum'] == sample_id), 'Ethnicity_self-report']
            eth = eth_series.values[0]
            if cohort == 'ISRAEL':
                if not isinstance(eth, str):
                    eth = 3
                else:
                    eth: int = 1 if eth.upper() == 'ASHKENAZI' else 2 if eth.upper() == 'MIZRAHI' else 3
            elif cohort == 'USA':
                if not isinstance(eth, str):
                    eth = 2
                else:
                    eth: int = 1 if eth.upper() == 'ASHKENAZI' else 2
            covar_df.loc[i, 'Ethn'] = eth

    valid_samples = covar_df['FID'].tolist()
    total_valid_samples.extend(valid_samples)
    # Make a list of the files that have len(valid_samples) == 0
    if len(valid_samples) == 0:
        empty_cohorts.append(cohort)
        continue
    else:
        full_cohorts.append(cohort)

    # Remove duplicates
    # samples = samples.drop_duplicates(subset=samples.columns[0])
    # sex_df = sex_df.drop_duplicates(subset=sex_df.columns[0])
    # pheno_df = pheno_df.drop_duplicates(subset=pheno_df.columns[0])
    # covar_df = covar_df.drop_duplicates(subset=covar_df.columns[0])

    # try:
    #     print("COHORT: " + cohort)
    #     print("Covar file sample size:" + str(len(covar_df)))
    # except TypeError:
    #     print("COHORT: " + cohort + " has no samples")
    # print("-----------------------------------------------")

    patients_df = covar_df[covar_df['Status'] == 2]
    controls_df = covar_df[covar_df['Status'] == 1]

    # Check if the subsets have data
    if not patients_df.empty and not controls_df.empty:
        # Calculate counts
        number_patients = pd.value_counts(patients_df['Status'])
        number_controls = pd.value_counts(controls_df['Status'])
        # Print the counts
        print(cohort + ': n(PD)=' + str(number_patients[2]) + '\tn(CONTROL)=' + str(number_controls[1]))
    else:
        print("No data found for one or both subsets.")
    # Write the files
    cohort_file_name = cohort.replace(' ', '')
    try:
        os.mkdir(os.path.join(output_path, 'cohort_' + cohort_file_name))
    except FileExistsError:
        pass

    samples.to_csv(os.path.join(output_path, 'cohort_' + cohort_file_name, cohort_file_name + '.samples.list'),
                   sep='\t', index=False, header=False)
    sex_df.to_csv(os.path.join(output_path, 'cohort_' + cohort_file_name, 'sex_' + cohort_file_name + '.txt'),
                  sep='\t', index=False)
    pheno_df.to_csv(os.path.join(output_path, 'cohort_' + cohort_file_name, 'pheno_' + cohort_file_name + '.txt'),
                    sep='\t', index=False)
    covar_df.to_csv(os.path.join(output_path, 'cohort_' + cohort_file_name, 'covar_' + cohort_file_name + '.txt'),
                    sep='\t', index=False)

print(len(df))
print(len(total_valid_samples))

# Create a new Excel file with the updated DataFrame
df = df[df['Snum'].isin(total_valid_samples)]
df.insert(0, 'key_cohort', df['Cohort'].apply(lambda x: [key for key, values in cohort_names.items() if x in values]))
df.to_excel('/Users/cemparlar/Desktop/Gan-Or/indexed_cohorts.xlsx', index=False)

print("The following cohorts have no samples in covar file; hence, they were excluded:")
print(empty_cohorts)
print("The following cohorts have samples in covar file:")
print(full_cohorts)