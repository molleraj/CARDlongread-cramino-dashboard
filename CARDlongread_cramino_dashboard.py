#!/usr/bin/env python3
# long read sequencing cramino QC dashboard generator
# generate violin and scatterplots from cramino QC summary table
# summary table generated by cramino report parser
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import statistics
import argparse
import dataclasses
# excel export with openpyxl
import openpyxl
# and xlsxwriter with pandas, though that could be modified to use openpyxl as well
import xlsxwriter
# for image saving
from io import BytesIO
# for date/time conversions
from datetime import datetime, timezone

# get summary statistics (min, max, range, mean, median, mode, and standard deviation for N50, sequence output, and flow cells per sample)
def get_summary_statistics(column):
    # define summary statistics class
    @dataclasses.dataclass
    class summary_statistics:
        total : int = 0
        min : float = 0
        max : float = 0
        range : float = 0
        mean : float = 0
        median : float = 0
        mode : float = 0
        stdev : float = 0
	# take single column of data frame as input
	# fill summary statistics class with each statistic
    # note that the line below should be used to give the total run count (across experiments and samples). Not useful with unique experiments, samples, etc.
    # as in flow cells per unique sample
    summary_statistics.total = len(column)
    summary_statistics.min = min(column)
    summary_statistics.max = max(column)
    summary_statistics.range = summary_statistics.max - summary_statistics.min
    summary_statistics.mean = statistics.mean(column)
    summary_statistics.median = statistics.median(column)
    summary_statistics.mode = statistics.mode(column)
    summary_statistics.stdev = statistics.stdev(column)
    # return data structure with each summary statistic as an attribute
    return summary_statistics

# make summary statistic data frame
def make_summary_statistics_data_frame(input_data_frame, property_names):
    # set column names
    column_names = ['Property', 'Total', 'Min', 'Max', 'Mean', 'Median', 'Mode', 'Standard Deviation']
    # initialize data frame
    summary_statistics_df = pd.DataFrame(index=property_names, columns=column_names)
    # iterate through each property and populate data frame with summary_statistics_set attributes
    for idx, name in enumerate(property_names):
        summary_statistics_set = get_summary_statistics(input_data_frame[name])
        summary_statistics_df.loc[name] = [property_names[idx],summary_statistics_set.total,summary_statistics_set.min,summary_statistics_set.max,summary_statistics_set.mean,summary_statistics_set.median,summary_statistics_set.mode,summary_statistics_set.stdev]
    # return populated summary statistics data frame
    return summary_statistics_df
    
# make violinplot/swarmplot figure worksheet in output workbook
def make_violinswarmplot_worksheet(data,input_variable,workbook,worksheet_name,x_axis_title=None,cutoff=None,title=None):
    # create worksheet for figure output
    worksheet=workbook.create_sheet(worksheet_name)
    # initialize raw data buffer for image
    imgdata=BytesIO()
    # initialize plot overall
    fig, ax = plt.subplots()
    # make swarm plot to show how data points overlap with distribution
    # replace color='black'
    ax = sb.swarmplot(data=data,x=input_variable,color='black')
    # add violin plot using seaborn (sb.violinplot)
    # increase transparency to improve swarmplot visibility
    ax = sb.violinplot(data=data,x=input_variable,color='white',ax=ax)
    # add x axis title if specified 
    if x_axis_title is not None:
        ax.set(xlabel=x_axis_title)
    # add title if specified
    if title is not None:
        ax.set_title(title)
    # add red line for 90GB/30X cutoff or whatever necessary for specific plots
    # cutoff line defined by x value
    if cutoff is not None:
        ax.axvline(x=cutoff,color='red')
    # put figure in variable to prep for saving into buffer
    # fig = swarmplot.get_figure()
    # save figure as 200 dpi PNG into buffer
    fig.savefig(imgdata, format='png', dpi=200)
    # close figure
    fig.clf()
    # make openpyxl image from raw data
    img = openpyxl.drawing.image.Image(imgdata)
    # set location of image in worksheet (A1)
    img.anchor = 'A1'
    # add image to worksheet
    worksheet.add_image(img)
    # close figure with matplotlib plt close
    plt.close()

# set up command line argument parser
parser = argparse.ArgumentParser(description='This program gets summary statistics from long read sequencing report data.')

# get input and output arguments
parser.add_argument('-input', action="store", dest="input_file", help="Input tab-delimited tsv file containing features extracted from long read sequencing reports.")
parser.add_argument('-output', action="store", dest="output_file", help="Output long read sequencing summary statistics XLSX")
parser.add_argument('-plot_title', action="store", default=None, dest="plot_title", help="Title for each plot in output XLSX (optional)")
# add boolean --plot_cutoff argument
parser.add_argument('--plot_cutoff', action=argparse.BooleanOptionalAction, default=True, dest="plot_cutoff", help="Include cutoff lines in violin plots (optional; default true; --no-plot_cutoff to override)")
# include failed run cutoff to exclude as well
parser.add_argument('-run_cutoff', action="store", default=1, type=float, dest="run_cutoff", help="Minimum data output per flow cell run to include (optional, 1 Gb default)")

# parse arguments
results = parser.parse_args()

# throw error if no input file provided
if results.input_file is None:
	quit('ERROR: No input file (-input) provided!')
        
# set default output filename
if results.output_file is None:
    results.output_file='output_summary_statistics.xlsx'

# read tab delimited output into pandas data frame
cramino_extract_initial=pd.read_csv(results.input_file,sep='\t')

# use functions above
# first filter out low output runs
cramino_extract = cramino_extract_initial[cramino_extract_initial['Yield (Gb)'] > results.run_cutoff]
# calculate median mapping q score
# not needed, add calculation to report parser
# cramino_extract['Median mapping Q score'] = -10*np.log10((100-cramino_extract['Median identity'])/100)
# cramino_extract['Mean mapping Q score'] = -10*np.log10((100-cramino_extract['Mean identity'])/100)
# calculate mean mapping q score
# fix indices
cramino_extract.reset_index(drop='True',inplace=True)
# summary statistics on...
cramino_summary_statistics_property_names=['Number of alignments','Percent of total reads','Yield (Gb)','Mean Coverage','Yield (Gb) [>25kb]','N50','N75','Median length','Mean length','Median identity','Mean identity','Median mapping Q score','Mean mapping Q score']
# make data frame
cramino_summary_statistics_df = make_summary_statistics_data_frame(cramino_extract, cramino_summary_statistics_property_names)
# output data frames and figures to excel spreadsheet
writer = pd.ExcelWriter(results.output_file)
# write data frames with a row between each
# write combined summary stats
start_row = 0
cramino_summary_statistics_df.to_excel(writer, startrow=start_row, index=False, sheet_name='Summary statistics report')
# close writer and save workbook
writer.close()
# then add svg figures
# use openpyxl and pipe image data into new worksheets
# append new worksheets to existing workbook
workbook = openpyxl.load_workbook(results.output_file)
# use make_violinswarmplot_worksheet function to insert figures
# plots that don't need cutoff included below to maintain order of spreadsheets (and order of figures)
# if/else depending on whether -plot_cutoff set
if results.plot_cutoff is True:
    # sequential violin/swarm plots for each property
    make_violinswarmplot_worksheet(cramino_extract,"Number of alignments",workbook,'Number of alignments plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Percent of total reads",workbook,'Percent of total reads plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Yield (Gb)",workbook,'Yield plot',None,90,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean Coverage",workbook,'Mean coverage plot',None,30,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Yield (Gb) [>25kb]",workbook,'Yield over 25 kb plot',None,90,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"N50",workbook,'N50 plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"N75",workbook,'N85 plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Median length",workbook,'Median length plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean length",workbook,'Mean length plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Median identity",workbook,'Median identity plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean identity",workbook,'Mean identity plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Median mapping Q score",workbook,'Median mapping Q score plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean mapping Q score",workbook,'Mean mapping Q score plot',None,None,results.plot_title)
else:
    # sequential violin/swarm plots for each property
    make_violinswarmplot_worksheet(cramino_extract,"Number of alignments",workbook,'Number of alignments plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Percent of total reads",workbook,'Percent of total reads plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Yield (Gb)",workbook,'Yield plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean Coverage",workbook,'Mean coverage plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Yield (Gb) [>25kb]",workbook,'Yield over 25 kb plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"N50",workbook,'N50 plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"N75",workbook,'N85 plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Median length",workbook,'Median length plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean length",workbook,'Mean length plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Median identity",workbook,'Median identity plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean identity",workbook,'Mean identity plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Median mapping Q score",workbook,'Median mapping Q score plot',None,None,results.plot_title)
    make_violinswarmplot_worksheet(cramino_extract,"Mean mapping Q score",workbook,'Mean mapping Q score plot',None,None,results.plot_title)
# save workbook when done
workbook.save(results.output_file)

# script complete
quit()
