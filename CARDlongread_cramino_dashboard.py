#!/usr/bin/env python3
# long read sequencing cramino QC dashboard generator
# generate violin and scatterplots from cramino QC summary table
# summary table generated by cramino report parser
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
def make_violinswarmplot_worksheet(data,input_variable,group_variable,legend_patches,user_palette,strip_plot_set,workbook,worksheet_name,x_axis_title=None,cutoff=None,title=None):
    # create worksheet for figure output
    worksheet=workbook.create_sheet(worksheet_name)
    # initialize raw data buffer for image
    imgdata=BytesIO()
    # initialize plot overall
    fig, ax = plt.subplots()
    # set up plots differently depending on whether group variable is set
    if (group_variable is None):
        # make swarm plot to show how data points overlap with distribution
        # replace color='black'
        if strip_plot_set is False:
            ax = sb.swarmplot(data=data,x=input_variable,color='black')
        elif strip_plot_set is True:
            ax = sb.stripplot(data=data,x=input_variable,color='black')
        # add violin plot using seaborn (sb.violinplot)
        # increase transparency to improve swarmplot visibility
        # use boxplot since only one "group" shown
        ax = sb.violinplot(data=data,x=input_variable,color='white',ax=ax)
    else:
        # make swarm plot to show how data points overlap with distribution
        # input variable on y axis and group on x axis
        # thus vertical swarm/violinplots instead of horizontal ones when no group specified
        # replace color='black' with hue set to group variable
        if strip_plot_set is False:
            # allow user set palette
            if user_palette is None:
                ax = sb.swarmplot(data=data,x=group_variable,y=input_variable,hue=group_variable,legend=False)
            else:
                ax = sb.swarmplot(data=data,x=group_variable,y=input_variable,hue=group_variable,palette=user_palette,legend=False)
        elif strip_plot_set is True:
            # allow user set palette
            if user_palette is None:
                ax = sb.stripplot(data=data,x=group_variable,y=input_variable,hue=group_variable,legend=False)
            else:
                ax = sb.stripplot(data=data,x=group_variable,y=input_variable,hue=group_variable,palette=user_palette,legend=False)
        # add violin plot using seaborn (sb.violinplot)
        # increase transparency to improve swarmplot visibility
        # include quartile lines in this context (for easily, visually comparing between groups)
        ax = sb.violinplot(data=data,x=group_variable,y=input_variable,color='white',inner="quartile",ax=ax)
    # add x axis title if specified 
    if x_axis_title is not None:
        ax.set(xlabel=x_axis_title)
    # add title if specified
    if title is not None:
        ax.set_title(title)
    # add red line for 90GB/30X cutoff or whatever necessary for specific plots
    # cutoff line defined by x value
    if cutoff is not None:
        # vertical line if non-grouped
        if (group_variable is None):
            ax.axvline(x=cutoff,color='red')
        # horizontal line if grouped (group variable on x-axis)
        else:
            ax.axhline(y=cutoff,color='red')
    # add legend with colors if requested
    if legend_patches is not None:
        plt.legend(handles=legend_patches)
    # put figure in variable to prep for saving into buffer
    # fig = swarmplot.get_figure()
    # save figure as 200 dpi PNG into buffer
    # tight layout to prevent titles from being cut off
    fig.savefig(imgdata, format='png', dpi=200, bbox_inches='tight')
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
parser.add_argument('-input', action="store", dest="input_file", nargs="+", help="Input tab-delimited tsv file containing features extracted from long read sequencing reports.")
# if multiple inputs, require input names
parser.add_argument('-names', action="store", default=None, dest="names", nargs="*", help="Names corresponding to input tsv file(s); required if more than one tsv provided.")
parser.add_argument('-output', action="store", dest="output_file", help="Output long read sequencing summary statistics XLSX")
parser.add_argument('-plot_title', action="store", default=None, dest="plot_title", help="Title for each plot in output XLSX (optional)")
# add boolean --plot_cutoff argument
parser.add_argument('--plot_cutoff', action=argparse.BooleanOptionalAction, default=True, dest="plot_cutoff", help="Include cutoff lines in violin plots (optional; default true; --no-plot_cutoff to override)")
# include failed run cutoff to exclude as well
parser.add_argument('-run_cutoff', action="store", default=1, type=float, dest="run_cutoff", help="Minimum data output per flow cell run to include (optional, 1 Gb default)")
# add option for stripplot instead of swarmplot (in case of excessive data points)
parser.add_argument('--strip_plot', action=argparse.BooleanOptionalAction, default=False, dest="strip_plot", help="Show strip plots instead of swarm plots inside violin plots (optional; default false)")
# add option for color palette
parser.add_argument('-colors', action="store", default=None, dest="colors", nargs="*", help="Color palette corresponding to sequential groups displayed (e.g., 'blue', 'red', 'blue'); optional and used only if more than one tsv provided.")
# add option for custom legend colors
parser.add_argument('-legend_colors', action="store", default=None, dest="legend_colors", nargs="*", help="Colors shown in the legend (e.g., 'blue', 'red'); optional and used only if more color palette included above. Must be palette subset.")
# add option for custom legend labels
parser.add_argument('-legend_labels', action="store", default=None, dest="legend_labels", nargs="*", help="Labels for each color in legend in order specified in -legend_colors.")
# add option to show sample size for groups in grouped violinplots
parser.add_argument('--group_count', action=argparse.BooleanOptionalAction, default=False, dest="show_group_count", help="Show group count in x-axis labels (optional; default false)")

# parse arguments
results = parser.parse_args()

# throw error if no input file provided
if results.input_file is None:
	quit('ERROR: No input file (-input) provided!')

# throw error if no names provided if multiple input files provided
if len(results.input_file)>1:
    if len(results.names)<=1:
        quit('ERROR: Multiple input files provided but not multiple names (-names).')
    elif len(results.names) != len(results.input_file):
        quit('ERROR: Number of names is different from number of input files.')
    # test if number of colors different from input files and names
    elif (results.colors is not None) and (len(results.colors) != len(results.names)):
        quit("ERROR: Color palette provided, but number of colors doesnt match number of group names.")
        
# set legend_patches to None by default
legend_patches=None
# test if legend colors and labels have proper length
if (results.legend_colors is not None) or (results.legend_labels is not None):
    if len(results.legend_colors) != len(results.legend_labels):
        quit('ERROR: Number of legend colors does not match number of legend labels.')
    else:
        # prepare legend patches
        # make list as long as legend colors (at this point same as legend_labels)
        legend_patches = [0] * len(results.legend_colors)
        for idx, i in enumerate(results.legend_colors):
            legend_patches[idx] = mpatches.Patch(color=i, label=results.legend_labels[idx])
        
# test if legend colors provided but not labels or vice versa
if ((results.legend_colors is not None) and (results.legend_labels is None)) or ((results.legend_colors is None) and (results.legend_labels is not None)):
    quit('ERROR: Either legend colors or legend labels provided but not both.')

# set default output filename
if results.output_file is None:
    results.output_file='output_summary_statistics.xlsx'
    
# read tab delimited output into pandas data frame
# case if just one input file provided
if len(results.input_file)==1:
    cramino_extract_initial=pd.read_csv(results.input_file[0],sep='\t')
    # first filter out low output runs
    cramino_extract = cramino_extract_initial[cramino_extract_initial['Yield (Gb)'] > results.run_cutoff]
    grouped=False
# what if multiple input files provided
elif len(results.input_file)>1:
    # store input tables in list as long input filename set
    cramino_extract_initial_list=[0] * len(results.input_file)
    for idx, i in enumerate(results.input_file): 
        cramino_extract_initial_list[idx]=pd.read_csv(i,sep='\t')
        # first filter out low output runs
        cramino_extract_initial_list[idx]=cramino_extract_initial_list[idx][cramino_extract_initial_list[idx]['Yield (Gb)'] > results.run_cutoff]
        # add group name to each table in list
        if results.show_group_count is True:
            # if group count specified, add group count to group name
            group_count = len(cramino_extract_initial_list[idx])
            # in this way, show n=717 or similar below group names in all plots
            cramino_extract_initial_list[idx]['Group']=results.names[idx]
            cramino_extract_initial_list[idx]['Group and count']=results.names[idx] + "\nn=" + str(group_count)
            # group_names_list[idx]=results.names[idx] + "\nn=" + str(group_count)
        else:
            cramino_extract_initial_list[idx]['Group']=results.names[idx]
    # combine groups into single concatenated data table
    cramino_extract=pd.concat(cramino_extract_initial_list[:],ignore_index=True)
    # set group variable
    grouped=True

# use functions above
# calculate median mapping q score
# not needed, add calculation to report parser
# cramino_extract['Median mapping Q score'] = -10*np.log10((100-cramino_extract['Median identity'])/100)
# cramino_extract['Mean mapping Q score'] = -10*np.log10((100-cramino_extract['Mean identity'])/100)
# calculate mean mapping q score
# fix indices
cramino_extract.reset_index(drop='True',inplace=True)
# summary statistics on...
cramino_summary_statistics_property_names=['Number of alignments','Percent of total reads','Yield (Gb)','Mean Coverage','Yield (Gb) [>25kb]','N50','N75','Median length','Mean length','Median identity','Mean identity','Median mapping Q score','Mean mapping Q score']
if grouped is False:
    # make data frame
    cramino_summary_statistics_df = make_summary_statistics_data_frame(cramino_extract, cramino_summary_statistics_property_names)
    # average median mapping Q score for plots below
    # average_median_mapping_q_score = cramino_summary_statistics_df.loc['Median mapping Q score','Mean']
    # average mean mapping Q score for plots below
    # average_mean_mapping_q_score = cramino_summary_statistics_df.loc['Mean mapping Q score','Mean']
    # output data frames and figures to excel spreadsheet
    writer = pd.ExcelWriter(results.output_file)
    # write data frames with a row between each
    # write combined summary stats
    start_row = 0
    cramino_summary_statistics_df.to_excel(writer, startrow=start_row, index=False, sheet_name='Summary statistics report')
    # close writer and save workbook
    writer.close()
elif grouped is True:
    # output data frames and figures to excel spreadsheet
    # intialize writer BEFORE for loop through results names
    writer = pd.ExcelWriter(results.output_file)
    # loop through all group names from input
    for idx, i in enumerate(results.names):
        # make data frame
        cramino_summary_statistics_df = make_summary_statistics_data_frame(cramino_extract[cramino_extract['Group'] == i], cramino_summary_statistics_property_names)
        # average median mapping Q score for plots below
        # average_median_mapping_q_score = cramino_summary_statistics_df.loc['Median mapping Q score','Mean']
        # average mean mapping Q score for plots below
        # average_mean_mapping_q_score = cramino_summary_statistics_df.loc['Mean mapping Q score','Mean']
        # write data frames with a row between each
        # write combined summary stats
        start_row = 0
        cramino_summary_statistics_df.to_excel(writer, startrow=start_row, index=False, sheet_name=i + ' statistics')
    # close writer and save workbook AFTER for loop complete
    writer.close()
# prepare legend patches
# then add svg figures
# use openpyxl and pipe image data into new worksheets
# append new worksheets to existing workbook
workbook = openpyxl.load_workbook(results.output_file)
# use make_violinswarmplot_worksheet function to insert figures
# plots that don't need cutoff included below to maintain order of spreadsheets (and order of figures)
# if/else depending on whether -plot_cutoff set
if results.plot_cutoff is True:
    # sequential violin/swarm plots for each property
    # simplify by iterating through for loop
    cramino_plot_cutoff_array=[None,None,90,30,90,None,None,None,None,None,None,None,None]
    cramino_plot_worksheet_names=['Number of alignments plot','Percent of total reads plot','Yield plot','Mean coverage plot','Yield over 25 kb plot','N50 plot','N75 plot','Median length plot','Mean length plot','Median identity plot','Mean identity plot','Median mapping Q score plot','Mean mapping Q score plot']
    for idx, i in enumerate(cramino_summary_statistics_property_names):
        # include group variable if necessary
        if grouped is False:
            make_violinswarmplot_worksheet(cramino_extract,i,None,legend_patches,results.colors,results.strip_plot,workbook,cramino_plot_worksheet_names[idx],None,cramino_plot_cutoff_array[idx],results.plot_title)
        elif grouped is True:
            if results.show_group_count is True:
                make_violinswarmplot_worksheet(cramino_extract,i,cramino_extract['Group and count'],legend_patches,results.colors,results.strip_plot,workbook,cramino_plot_worksheet_names[idx],None,cramino_plot_cutoff_array[idx],results.plot_title)
            else:
                make_violinswarmplot_worksheet(cramino_extract,i,cramino_extract['Group'],legend_patches,results.colors,results.strip_plot,workbook,cramino_plot_worksheet_names[idx],None,cramino_plot_cutoff_array[idx],results.plot_title)
else:
    # sequential violin/swarm plots for each property
    # simplify by iterating through for loop
    # no need for cutoff array this time
    cramino_plot_worksheet_names=['Number of alignments plot','Percent of total reads plot','Yield plot','Mean coverage plot','Yield over 25 kb plot','N50 plot','N75 plot','Median length plot','Mean length plot','Median identity plot','Mean identity plot','Median mapping Q score plot','Mean mapping Q score plot']
    for idx, i in enumerate(cramino_summary_statistics_property_names):
        # include group variable if necessary
        if grouped is False:
            make_violinswarmplot_worksheet(cramino_extract,i,None,legend_patches,results.colors,results.strip_plot,workbook,cramino_plot_worksheet_names[idx],None,None,results.plot_title)
        elif grouped is True:
            if results.show_group_count is True:
                make_violinswarmplot_worksheet(cramino_extract,i,cramino_extract['Group and count'],legend_patches,results.colors,results.strip_plot,workbook,cramino_plot_worksheet_names[idx],None,None,results.plot_title)
            else:
                make_violinswarmplot_worksheet(cramino_extract,i,cramino_extract['Group'],legend_patches,results.colors,results.strip_plot,workbook,cramino_plot_worksheet_names[idx],None,None,results.plot_title)
# save workbook when done
workbook.save(results.output_file)

# script complete
quit()
