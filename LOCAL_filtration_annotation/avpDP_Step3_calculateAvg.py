file_path = '/Users/cemparlar/Desktop/Gan-Or/*JULY2024_yesnoG2019S_LRRK2/average depth calculation/average_depth_test2.txt'

# Initialize variables
cohort_depths = {
    'USA': [],
    'ISRAEL': [],
    'FRENCH': []
}

current_cohort = None

# Open and read the file
with open(file_path, 'r') as file:
    for line in file:
        line = line.strip()
        # Check if the line indicates the start of a cohort
        if line.startswith('Processing cohort:'):
            current_cohort = line.split(':')[-1].strip()
            continue

        # Check if the line indicates the end of a cohort's depths
        if current_cohort and line.startswith(f'Cohort {current_cohort} average depth'):
            current_cohort = None
            continue

        # Check if the line contains the depth information and add it to the current cohort
        if line.startswith('Depth for') and current_cohort in cohort_depths:
            # Split the line to extract the depth information
            parts = line.split(':')
            if len(parts) == 2:
                try:
                    depth = float(parts[1].strip())
                    cohort_depths[current_cohort].append(depth)
                except ValueError:
                    # Skip lines that do not contain valid depth values
                    continue

# Calculate the average depth for each cohort
average_depths = {}
for cohort, depths in cohort_depths.items():
    if depths:
        average_depths[cohort] = sum(depths) / len(depths)
    else:
        average_depths[cohort] = 0

average_depths
print(average_depths)