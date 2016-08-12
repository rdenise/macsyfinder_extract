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
    functions=np.loadtxt(file_function, , dtype='string', delimiter="\t")
    for function in functions:
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

    upper=np.mean(numpy_len)+2*np.std(numpy_len)
    lower=np.mean(numpy_len)+2*np.std(numpy_len)

    if lower < 0 :
        lower=0

    return [lower, upper]

##########################################################################################

def set_dict_cutoff(listOfFasta):
    cutoff_dict={}
    array_for_file=[]
	for Fastafile in listOfFasta :
		current_file = Fastafile.split("/")[-1]
        cutoff_dict[current_file]=set_cutoff(fasta_file)
        array_for_file.append([current_file]+cutoff_dict[current_file])

    header=["file", "lower_cutoff", "upper_cutoff"]
    np.savetxt(fid, self.classes ,delimiter="\t", fmt="%.2f", header=header)  
