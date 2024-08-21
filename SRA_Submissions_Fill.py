#!/usr/bin/env python3
import os
import glob
import pandas as pd
import numpy as np
import xlsxwriter

#curl -L https://oregonstate.box.com/s/mwkowifdk8hda4y0k73xdt458qbq1a2z --output SC2.xlsx

# List of samples to exclude from submissions
file_patterns = [
    "*Twist*",
    "*water*",
    "*ONT-PSN*",
    "*SIL-TR*",
    "*WS-SSL*",
    "*WS-Inf*",
    "*FLO-10th*",
    "*FLO-15th*",
    "*FLO-16th*",
    "*FLO-28th*",
    "*FLO-2nd*",
    "*FLO-36th*",
    "*FLO-38th*",
    "*FLO-Boy*",
    "*FLO-Hospital*",
    "*FLO-SVPS*",
    "*SALM-Eff*",
    "*COCC-Dorm*"
]

# Remove all files matching excluded samples
for file_pattern in file_patterns:
    for file_name in glob.glob(file_pattern):
        os.remove(file_name)

# use glob to get a list of files matching the R1_001.fastq.gz pattern
r1_files = sorted(glob.glob('*R1_001.fastq.gz'))

# use glob to get a list of files matching the R2_001.fastq.gz pattern
r2_files = sorted(glob.glob('*R2_001.fastq.gz'))

# create a list of tuples, where each tuple contains the R1 and R2 file names
file_list = list(zip(r1_files, r2_files))

# create a tab-separated file with the file names in two columns, sorted by both columns
with open('file_list_sorted.tsv', 'w') as f:
    for r1_file, r2_file in file_list:
        f.write(f"{r1_file}\t{r2_file}\n")

# create a new file list from the original file list
# by splitting each file name at the 6th "-" and the first "_"
new_file_list = []
for file_pair in file_list:
    new_pair = []
    for file in file_pair:
        parts = file.split("-")
        if len(parts) > 6:
            new_pair.append("-".join(parts[6:]))
    new_file_list.append(tuple(new_pair))

samples_list = []
for file_pair in new_file_list:
    samples_pair = []
    for file in file_pair:
        parts = file.split("_")
        if len(parts) > 0:
            samples_pair.append(parts[0])
    samples_list.append(samples_pair[0])

# read SuppDataTRACE.txt
trace_data = pd.read_csv('SuppDataTRACE.txt', sep='\t')

# create a dictionary for faster lookup
trace_dict = dict(zip(trace_data['MasterID'], trace_data.iloc[:, 0]))

# create a dictionary for 'Flow_LPerDay' and 'Population' using 'MasterID' as the key
flow_dict = dict(zip(trace_data['MasterID'], trace_data.iloc[:, 3]))
population_dict = dict(zip(trace_data['MasterID'], trace_data.iloc[:, 4]))

# prepend corresponding value from column 0 of SuppDataTRACE
final_samples_list = []
for sample in samples_list:
    master_id = sample.split("-")[0]
    if master_id in trace_dict:
        final_samples_list.append(str(trace_dict[master_id]) + "-" + sample)

#create dates column
dates = []
for sample in samples_list:
    parts = sample.split('-')
    if len(parts) >= 5:
        date = '-'.join(parts[2:5])
        dates.append(date)
    else:
        dates.append('')

# Create a new list called final_samples_list_with_population
final_samples_list_with_population = []
for sample in samples_list:
    master_id = sample.split("-")[0]
    population = population_dict.get(master_id, '')
    population_str = str(population) if population else ''  # Convert population to str, or use an empty string if it's None
    flow = flow_dict.get(master_id, '')
    flow_str = str(flow) if flow else ''  # Convert flow to str, or use an empty string if it's None
    trace_value = trace_dict.get(master_id, '')
    trace_value_str = str(trace_value) if trace_value else ''  # Convert trace value to str, or use an empty string if it's None
    final_samples_list_with_population.append([trace_value_str + '-' + sample, sample, population_str, flow_str])  # Prepend trace value and append sample, population, and flow as a list

# Copy the contents of the dates list to a new column in final_samples_list_with_population
for i in range(len(final_samples_list_with_population)):
    if i < len(dates):
        final_samples_list_with_population[i].append(dates[i])
    else:
        final_samples_list_with_population[i].append('')  # Add empty string if there is no corresponding date

# Read the TSV file and store the contents in a list
tsv_file = 'file_list_sorted.tsv'
tsv_data = []
with open(tsv_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            tsv_data.append(line.split('\t'))

# Add the TSV data as new columns in final_samples_list_with_population
for i in range(len(final_samples_list_with_population)):
    if i < len(tsv_data):
        final_samples_list_with_population[i].extend(tsv_data[i])
    else:
        final_samples_list_with_population[i].extend(['', ''])  
# Add empty strings if there is no corresponding TSV data

# Define the headers for the TSV file
headers = ['sample', 'ww_sample_id', 'population', 'flow', 'date', 'r1', 'r2']

# Create a list of lines for the TSV file, starting with the headers
tsv_lines = [headers] + final_samples_list_with_population

# Save the data as a TSV file
tsv_file = 'sampledata.tsv'
with open(tsv_file, 'w') as f:
    for line in tsv_lines:
        f.write('\t'.join(line) + '\n')

# Read the existing TSV file
existing_tsv_file = 'sampledata.tsv'
existing_tsv_data = pd.read_csv(existing_tsv_file, sep='\t')


# Read the new data from 'SRA_Metadata_blank.txt' and 'SARS-CoV-2.wwsurv.1.0.Blank.txt'
sra_data = pd.read_csv('SRA_Metadata_blank.txt', sep='\t', encoding='cp1252')
wwsurv_data = pd.read_csv('SARS-CoV-2.wwsurv.1.0.Blank.txt', sep='\t', encoding='cp1252')

# Extract the data from existing_tsv_data column 1 (index 0) from row 2 to last row (N)
#existing_tsv_column_data = existing_tsv_data.iloc[1:, 0].values
existing_tsv_column_data_samplenames = existing_tsv_data.iloc[0:, 0]
existing_tsv_column_data_sampleids = existing_tsv_data.iloc[0:, 1]
existing_tsv_column_data_r1 = existing_tsv_data.iloc[0:, 5]
existing_tsv_column_data_r2 = existing_tsv_data.iloc[0:, 6]
existing_tsv_column_data_population = existing_tsv_data.iloc[0:, 2]
existing_tsv_column_data_flow = existing_tsv_data.iloc[0:, 3]
existing_tsv_column_data_date = existing_tsv_data.iloc[0:, 4]



# Use xlsxwriter to write existing_tsv_column_data to sra_data
sra_data_output_file = 'SRA_Metadata_updated.xlsx'
writer = pd.ExcelWriter(sra_data_output_file, engine='xlsxwriter')
sra_data.to_excel(writer, sheet_name='Sheet1', index=False)

# Create a pandas DataFrame for existing_tsv_column_data with a dummy header to align with sra_data
existing_tsv_column_data_df = pd.DataFrame({'header': existing_tsv_column_data_samplenames})
existing_tsv_column_data_df_ids = pd.DataFrame({'header': existing_tsv_column_data_sampleids})
existing_tsv_column_data_df_r1 = pd.DataFrame({'header': existing_tsv_column_data_r1})
existing_tsv_column_data_df_r2 = pd.DataFrame({'header': existing_tsv_column_data_r2})
existing_tsv_column_data_df_population = pd.DataFrame({'header': existing_tsv_column_data_population})
existing_tsv_column_data_df_flow = pd.DataFrame({'header': existing_tsv_column_data_flow})
existing_tsv_column_data_df_date  = pd.DataFrame({'header': existing_tsv_column_data_date})


# Write existing_tsv_column_data_df to sra_data in column 2, starting from row 2
existing_tsv_column_data_df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, startcol=1, header=False)

# Write existing_tsv_column_data_df to sra_data in column 1, starting from row 2
existing_tsv_column_data_df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, startcol=0, header=False)

# Write existing_tsv_column_data_df to sra_data in column 3, starting from row 2
existing_tsv_column_data_df_ids.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, startcol=2, header=False)

# Write existing_tsv_column_data_df to sra_data in column 11, starting from row 2
existing_tsv_column_data_df_r1.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, startcol=12, header=False)

# Write existing_tsv_column_data_df to sra_data in column 12, starting from row 2
existing_tsv_column_data_df_r2.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, startcol=13, header=False)

writer.close()

# Use xlsxxwriter to write existing_tsv_column_data to sars_data
sars_data_output_file = 'SARS-CoV-2.wwsurv.1.0.updated.xlsx'
writer = pd.ExcelWriter(sars_data_output_file, engine='xlsxwriter')
wwsurv_data.to_excel(writer, sheet_name='Sheet1', index=False)

# Write existing_tsv_column_data_df to sra_data in column 1, starting from row 13
existing_tsv_column_data_df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=12, startcol=0, header=False)

# Write existing_tsv_column_data_df to sra_data in column 5, starting from row 13
existing_tsv_column_data_df_date.to_excel(writer, sheet_name='Sheet1', index=False, startrow=12, startcol=4, header=False)

# Write existing_tsv_column_data_df to sra_data in column 8, starting from row 13
existing_tsv_column_data_df_population.to_excel(writer, sheet_name='Sheet1', index=False, startrow=12, startcol=7, header=False)

# Write existing_tsv_column_data_df to sra_data in column 26, starting from row 13
existing_tsv_column_data_df_flow.to_excel(writer, sheet_name='Sheet1', index=False, startrow=12, startcol=25, header=False)

writer.close()

