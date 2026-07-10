#!/usr/bin/env python3
# long read sequencing cramino QC parser
# convert cramino reports into single summary table per cohort
import glob
import pandas as pd
import numpy as np
import argparse
import dataclasses
from dateutil.parser import isoparse
from pandas.errors import EmptyDataError

# example cramino output
# File name	Chile_404.sorted_meth.bam
# Number of alignments	7876502
# % from total reads	77.57
# Yield [Gb]	147.44
# Mean coverage	47.56
# Yield [Gb] (>25kb)	65.65
# N50	23260
# N75	16395
# Median length	16935.00
# Mean length	18718
# Median identity	99.47
# Mean identity	97.80
# Path	/data/CARD_AUX/LRS_temp/CHILE/MAPPED_BAM/Chile_404.sorted_meth.bam
# Creation time	NA

# handle command line arguments
# get cramino output file list
# function to pull fields out of cramino output dataframe
def get_fields_from_cramino(input_cramino_df,bam_type):
    # define fields_from_cramino class
    @dataclasses.dataclass
    class fields_from_cramino:
        file_name : ''
        number_of_alignments : float = 0
        percent_of_total_reads : float = 0
        yield_gb : float = 0
        mean_coverage : float = 0
        yield_gb_over_25kb : float = 0
        n50 : float = 0
        n75 : float = 0
        median_length : float = 0
        mean_length : float = 0
        median_identity : float = 0
        mean_identity : float = 0
        median_identity_q_score : float = 0
        mean_identity_q_score : float = 0
    # get values from each consecutive field in the cramino output dataframe
    # what to do if data frame is empty
    if (len(input_cramino_df[1])==1):
        fields_from_cramino.file_name = input_cramino_df[1][0]
        fields_from_cramino.number_of_alignments = 0
        fields_from_cramino.percent_of_total_reads = 0
        fields_from_cramino.yield_gb = 0
        fields_from_cramino.mean_coverage = 0
        fields_from_cramino.yield_gb_over_25kb = 0
        fields_from_cramino.n50 = 0
        fields_from_cramino.n75 = 0
        fields_from_cramino.median_length = 0
        fields_from_cramino.mean_length = 0
        fields_from_cramino.median_identity = 0
        fields_from_cramino.mean_identity = 0
        fields_from_cramino.median_identity_q_score = 0
        fields_from_cramino.mean_identity_q_score = 0
    # as shown above, these are just sequential values (row by row) in the table
    else:
        fields_from_cramino.file_name = input_cramino_df[1][0]
        fields_from_cramino.number_of_alignments = input_cramino_df[1][1]
        fields_from_cramino.percent_of_total_reads = input_cramino_df[1][2]
        fields_from_cramino.yield_gb = input_cramino_df[1][3]
        fields_from_cramino.mean_coverage = input_cramino_df[1][4]
        fields_from_cramino.yield_gb_over_25kb = input_cramino_df[1][5]
        fields_from_cramino.n50 = input_cramino_df[1][6]
        fields_from_cramino.n75 = input_cramino_df[1][7]
        fields_from_cramino.median_length = input_cramino_df[1][8]
        fields_from_cramino.mean_length = input_cramino_df[1][9]
        if (bam_type == "mapped_bam"):
            fields_from_cramino.median_identity = input_cramino_df[1][10]
            fields_from_cramino.mean_identity = input_cramino_df[1][11]
            fields_from_cramino.median_identity_q_score = round(-10*np.log10((100-float(fields_from_cramino.median_identity))/100),2)
            fields_from_cramino.mean_identity_q_score = round(-10*np.log10((100-float(fields_from_cramino.mean_identity))/100),2)
        elif (bam_type == "unmapped_bam"):
            # none of the below provided for cramino run with --ubam option (apparently identity information not provided)
            fields_from_cramino.median_identity = 0
            fields_from_cramino.mean_identity = 0
            fields_from_cramino.median_identity_q_score = 0
            fields_from_cramino.mean_identity_q_score = 0
    return fields_from_cramino
# load json file list
# user input
inparser = argparse.ArgumentParser(description = 'Extract data from long read cramino mapping QC reports into summary table')
inparser.add_argument('--bam_type', default=None, choices=['mapped_bam','unmapped_bam'], required=True, type=str, help = 'cramino report type (mapped BAM report includes identity values).')
inparser.add_argument('--cramino_dir', default=None, type=str, help = 'path to directory containing cramino files, if converting whole directory')
inparser.add_argument('--filelist', default=None, type=str, help = 'text file containing list of all cramino reports to parse')
inparser.add_argument('--output', action="store", type=str, dest="output_file", help="Output long read cramino report summary table in tab-delimited format")
args = inparser.parse_args()
# get list of files
if args.cramino_dir is not None:
    files = glob.glob(f'{args.cramino_dir}/*.txt')
elif args.filelist is not None:
    with open(args.filelist, 'r') as infile:
        files = [x.strip() for x in infile.readlines()]
else:
    quit('ERROR: No directory (--cramino_dir) or file list (--filelist) provided!')
# create output data frame
# set indices
cramino_report_df_indices = [np.arange(0,len(files))]
# set column names
cramino_report_column_names = ['Filename','Number of alignments','Percent of total reads','Yield (Gb)','Mean Coverage','Yield (Gb) [>25kb]','N50','N75','Median length','Mean length','Median identity','Mean identity','Median identity Q score','Mean identity Q score']
# initialize data frame with said column names and filenames as indexes
cramino_report_df = pd.DataFrame(index=cramino_report_df_indices,columns=cramino_report_column_names)
# main loop to process files
for idx, x in enumerate(files):
    try:
        # cramino tsv file
        # debug by printing tsv file to stdout
        # print(x)
        f = open(x, "r")
        # Reading data frame from TSV cramino output
        # read tab delimited output into pandas data frame
        data=pd.read_csv(f,sep='\t',header=None)
        # get important information
        current_data_fields = get_fields_from_cramino(data,args.bam_type)
        cramino_report_df.loc[idx] = [current_data_fields.file_name,current_data_fields.number_of_alignments,current_data_fields.percent_of_total_reads,current_data_fields.yield_gb,current_data_fields.mean_coverage,current_data_fields.yield_gb_over_25kb,current_data_fields.n50,current_data_fields.n75,current_data_fields.median_length,current_data_fields.mean_length,current_data_fields.median_identity,current_data_fields.mean_identity,current_data_fields.median_identity_q_score,current_data_fields.mean_identity_q_score]
        f.close()
    except EmptyDataError:
        print(f"No columns to parse from file {x}")
        continue
    except ValueError as e:
        print(e)
        continue
# remove empty rows
cramino_report_df.dropna(how='all',inplace=True)
# print output data frame to tab delimited tsv file
cramino_report_df.to_csv(args.output_file,sep='\t',index=False)
# end program
quit()    

