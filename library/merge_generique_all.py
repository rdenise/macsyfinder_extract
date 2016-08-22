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

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def write_merge_file(generique_report_file, all_report_file, write_file) :

	"""
	Function that write a new file that is the association of the two file without line
	that are present in both
	
	:param generique_report_file: name of the macsyfinder report for generique.xml
	:type: str
	:param all_report_file: name of the macsyfinder report for T2SS.xml, T4P.xml, Tad.xml, Archaellum.xml
	:type: str
	:param write_file: name of the file we want to write
	:type: str
	"""

	with open(all_report_file, 'r') as read_report :
		with open(write_file, 'w') as w_file :
			report_generique = np.genfromtxt(generique_report_file, dtype='string', delimiter="\t")

			for line in read_report :

				split_line = line.rstrip("\r\n").split("\t")

				if split_line[0] in report_generique[:,0] :
					index = report_generique[:,0].tolist().index(split_line[0])

					if split_line[2] == report_generique[index][2] :
						w_file.write(line)

					else :
						print "**********"
						print "ERROR::Not normal"
						print line
						print "\t".join(report_generique[index,:].tolist())

					report_generique = np.delete(report_generique, index, axis=0)

				else:
					w_file.write(line)

			for line in report_generique :
				w_file.write("\t".join(line.tolist())+"\n")
