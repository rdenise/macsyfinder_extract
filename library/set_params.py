# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Libraries used
##
##########################################################################################
##########################################################################################

import numpy as np
from Bio import SeqIO
import glob
import os, sys

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def create_folder(mypath):

	"""
	Created the folder that I need to store my result if it doesn't exist
	:param mypath: path where I want the folder (write at the end of the path)
	:type: string
	:return: Nothing
	"""

	try:
		os.makedirs(mypath)
	except OSError:
		pass

##########################################################################################
##########################################################################################

def read_protein_function(file_function):

    """
    Function that read a tabulate file with the function on the first column
    and the name of the system_protein on the other (ex : ATPase	T2SS_gspE)

    :param file_function: the full path of the file to read
    :type: string
    :return: a dictionnary with the name of the protein in key and the name of
    the function as value
    :rtype: dict
    """

    dict_function={}
    with open(file_function,"r") as r_file:
        functions=[line.split() for line in r_file]
    for function in functions:
        if function[0][0] != "#" :
            dict_function.update({protein:function[0] for protein in function[1:]})
    return dict_function

##########################################################################################

def set_cutoff(fasta_file):

    """
    Function that calculate the cutoff for the file gave as argument.

    :param fasta_file: name of the concatenated fasta file
    :type: string
    :return: a list with the lower and upper cutoff
    :rtype: list
    """

    list_len=[]
    seqiter = SeqIO.parse(open(fasta_file), 'fasta')
    for seq in seqiter :
        list_len.append(len(seq))
    numpy_len=np.array(list_len)

    upper=int(np.mean(numpy_len)+2*np.std(numpy_len))
    lower=int(np.mean(numpy_len)-2*np.std(numpy_len))

    if lower < 0 :
        lower=0

    return [lower, upper]

##########################################################################################

def set_dict_cutoff_init(listOfFasta, INFO_folder):

    """
    Function used create the cutoff dictionnary for each concatenated fasta and
    write the cutoff file in the folder of the analysis.

    :param listOfFasta: List of all the file fasta where we need to remove sequences
    :type: list of string
    :param INFO_folder: the absolute path to the info folder
    :type: string
    :return: the cutoff dictionnary set with the name of the file in key and the
    lower and upper cutoff as value
    :rtype: dict
    """

    fid=os.path.join(INFO_folder,"proteins.cutoff")
    cutoff_dict={}
    array_for_file=[]
    for fastafile in listOfFasta :
    	current_file = os.path.basename(fastafile)
        cutoff_dict[current_file]=set_cutoff(fastafile)
        array_for_file.append([current_file]+cutoff_dict[current_file])

    header="\t".join(["file", "lower_cutoff", "upper_cutoff"])
    np.savetxt(fid, np.array(array_for_file) ,delimiter="\t", fmt="%s", header=header)
    return cutoff_dict

##########################################################################################

def set_dict_cutoff(cutoff_file):

    """
    Function used create the cutoff dictionnary if the file exist (file give in
    argument).

    :param cutoff_file: File with the name of the file in first column and the
    lower and upper cutoff on the other
    :type: string
    :return: the cutoff dictionnary set with the name of the file in key and the
    lower and upper cutoff as value
    :rtype: dict
    """

    tab_numpy=np.loadtxt(cutoff_file, delimiter="\t", skiprows=1, dtype=np.object)
    cutoff_dict={line[0]:line[1:] for line in tab_numpy}

    return cutoff_dict
