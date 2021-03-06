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
import pandas as pd

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

	return

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
##########################################################################################

def read_protein_function_reverse(file_function):

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
		for line in r_file :
			if not line.startswith("#") :
				line_split = line.split()
				dict_function[line_split[0]] = line_split[1:]
	return dict_function

##########################################################################################
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
##########################################################################################

def set_dict_cutoff(cutoff_file):

	"""
	Function used to create the cutoff dictionnary if the file exist (file give in
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

##########################################################################################
##########################################################################################

def create_dict_system(PROTEIN_FUNCTION):

	"""
	Function used to create the dictionnary that contains the name of all the
	systems in key and the list of all the protein of this systems in keys

	:param cutoff_file: dictionnary that contains the information of the function
	name
	:type: dict
	:return: dictionnary that contains the name of all the
	systems in key and the list of all the protein of this systems in keys
	:rtype: dict
	"""

	dict_system = {}

	for keys in PROTEIN_FUNCTION :
		system, *protein = keys.split("_")
		if system in dict_system :
			dict_system[system].append(keys)
		else :
			dict_system[system] = [keys]

	return dict_system

##########################################################################################
##########################################################################################

def create_dict_wanted(file_wanted):

	"""
	Function that create the dictionary of the phylum wanted.

	:param list_wanted: name of the file of phylum wanted
	:type: str
	:return: dictionary with the kingdom as key and the list of phylum as value.
	:rtype: dict
	"""

	info = np.genfromtxt(file_wanted, dtype="str", delimiter="\t")
	dict_wanted = {kingdom:[] for kingdom in np.unique(info[:,0])}
	for kingdom, phylum in info:
		dict_wanted[kingdom].append(phylum)

	return dict_wanted, info[:,1].tolist()

##########################################################################################
##########################################################################################

def create_dict_distance(file_distance):

	"""
	Function that create the dictionary of the maximal distance within systems.

	:param file_distance: name of the file of distance in systems
	:type: str
	:return: dictionary with the systems as key and the distance as value.
	:rtype: dict
	"""

	distance_df = pd.read_table(file_distance, comment="#", index_col=0, names=["distance"])

	return distance_df.distance.to_dict()

##########################################################################################
##########################################################################################

def read_systems_found(systems_found_file):

    """
    Read .names of macsyfinder_extract and created a dataframe with the rigth type on each columns

    :param systems_found_file: Name of the macsyfinder_extract systems_found.names file
    :type: str
    :return: the dataframe completed
    :rtype: pandas.Dataframe
    """

    names_dataframe=["Species_Id","Replicon_Id","System_name","System_status","System_number","Proteins","Kingdom","Phylum","Lineage"]
    summary_df = pd.read_table(systems_found_file, names=names_dataframe, comment="#")
    summary_df.ix[:,5]=summary_df.ix[:,5].apply(eval)

    return summary_df

##########################################################################################
##########################################################################################

def set_dict_protein(file_protein_wanted):

    """
    Set the dictionary of protein wanted in the analysis

    :param file_protein_wanted: the name of the file protein_function.def that is a tsv file with the name of the protein function in the first column and the name of the protein in the other ones
    :type: str
    :return: a dictionary with the name of all the protein of the analysis in value and the function of the protein as key
    :rtype: dict
    """

    with open(file_protein_wanted, "r") as r_file :
        #Je fais des lignes sans \n
        split_file = r_file.read().splitlines()
        #Je split les lignes totalement ex : [ATPase, pilA, ...]
        split_file = [line.split() for line in split_file if not line.startswith("#")]
        #Je crée mon dictionaire
        dict_protein = {line[0]:sorted(line[1:]) for line in split_file}

    return dict_protein
