# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Libraries used
##
##########################################################################################
##########################################################################################

import numpy as np
import sys
import os
import pandas as pd
import subprocess
from set_params import create_folder

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def write_merge_file_old(generic_report_file, all_report_file, write_file) :

	"""
	Function that write a new file that is the association of the two file without line
	that are present in both

	:param generic_report_file: name of the macsyfinder report for generique.xml
	:type: str
	:param all_report_file: name of the macsyfinder report for T2SS.xml, T4P.xml, Tad.xml, Archaellum.xml
	:type: str
	:param write_file: name of the file we want to write
	:type: str
	"""

	length = len(open(all_report_file, 'rt').readlines())
	list_not_generic = []

	with open(all_report_file, 'rt') as read_report :
		with open(write_file, 'w') as w_file :
			names_dataframe=['Hit_Id','Replicon_name','Position','Sequence_length','Gene','Reference_system','Predicted_system','System_Id','System_status','Gene_status','i-evalue','Score','Profile_coverage','Sequence_coverage','Begin_match','End_match']
			report_generic = pd.read_table(generic_report_file, names=names_dataframe)

			print()
			progression = 1

			for line in read_report :
				sys.stdout.write("{:.2f}% : {}/{} line of report read\r".format(progression/float(length)*100, progression,length))
				sys.stdout.flush()
				progression += 1

				split_line = line.rstrip().split("\t")

				if split_line[0] in list(report_generic.Hit_Id) :
					index = list(report_generic.Hit_Id).index(split_line[0])

					if split_line[2] == str(report_generic.iloc[index, 2]) :
						w_file.write(line)

					else :
						print()
						print("**********")
						print("ERROR::Not normal")
						print(line)
						print("\t".join(report_generic.iloc[index,:].values.tolist()))
						sys.exit(0)

					#report_generic = np.delete(report_generic, index, axis=0)
					list_not_generic.append(report_generic.iloc[index, 7])

				else:
					w_file.write(line)

			report_generic = report_generic[~report_generic.System_Id.isin(list_not_generic)]
			report_generic.to_csv(w_file, sep="\t", header=False, index=False)

	print()

	putative_summary_all = all_report_file.replace("report", "summary")
	putative_summary_generic = generic_report_file.replace("report", "summary")

	if os.path.exists(putative_summary_all) and os.path.exists(putative_summary_generic) :
		status = subprocess.call("cat {} > {}".format(putative_summary_all, write_file.replace("report", "summary")), shell=True)
		write_summary_file(putative_summary_generic, list(report_generic.System_Id), write_file.replace("report", "summary"))

	return

##########################################################################################
##########################################################################################

def write_merge_file(generic_report_file, all_report_file, write_file) :

	"""
	Function that write a new file that is the association of the two file without line
	that are present in both

	:param generic_report_file: name of the macsyfinder report for generique.xml
	:type: str
	:param all_report_file: name of the macsyfinder report for T2SS.xml, T4P.xml, Tad.xml, Archaellum.xml
	:type: str
	:param write_file: name of the file we want to write
	:type: str
	"""

	length = len(open(all_report_file, 'rt').readlines())
	list_not_generic = []

	with open(write_file, 'w') as w_file :

		names_dataframe=['Hit_Id','Replicon_name','Position','Sequence_length','Gene','Reference_system','Predicted_system','System_Id','System_status','Gene_status','i-evalue','Score','Profile_coverage','Sequence_coverage','Begin_match','End_match']
		report_generic = pd.read_table(generic_report_file, names=names_dataframe, dtype="str")
		report_other = pd.read_table(all_report_file, names=names_dataframe, dtype="str")

		list_not_generic = report_generic[report_generic.Hit_Id.isin(report_other.Hit_Id)].System_Id.unique()

		report_generic = report_generic[~report_generic.System_Id.isin(list_not_generic)]

		report_concat = pd.concat([report_other, report_generic])
		report_concat.sort_values(['System_Id', 'Hit_Id'], inplace=True)

		report_concat.to_csv(w_file, sep="\t", header=False, index=False)

	putative_summary_all = all_report_file.replace("report", "summary")
	putative_summary_generic = generic_report_file.replace("report", "summary")

	if os.path.exists(putative_summary_all) and os.path.exists(putative_summary_generic) :
		status = subprocess.call("cat {} > {}".format(putative_summary_all, write_file.replace("report", "summary")), shell=True)
		write_summary_file(putative_summary_generic, list(report_generic.System_Id), write_file.replace("report", "summary"))

	print("\nDone!")

	return

##########################################################################################
##########################################################################################


def write_summary_file(generic_summary_file, choose_generic, write_file) :

	"""
	Function that write a new file that is the association of the two file without line
	that are present in both

	:param generic_report_file: name of the macsyfinder summary for generique.xml
	:type: str
	:param all_report_file: list of all the name of generic kept
	:type: list of str
	:param write_file: name of the file we want to write
	:type: str
	"""

	with open(write_file, "a") as w_file :
		all_summary_generic = pd.read_table(generic_summary_file, header=None)
		all_summary_generic[all_summary_generic.iloc[:,1].isin(choose_generic)].to_csv(w_file, sep="\t", header=False, index=False)

	return
