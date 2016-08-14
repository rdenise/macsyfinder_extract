# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Libraries used
##
##########################################################################################
##########################################################################################

from Bio import SeqIO
import sys
import numpy as np
import rpy2.robjects as robjects
import re
import os
from set_params import *

##########################################################################################
##########################################################################################
##
##								Global variable
##
##########################################################################################
##########################################################################################

PROTEIN_FUNCTION = read_protein_function(os.path.join(os.path.dirname(os.path.realpath(__file__)), "definition_file/protein_function.def"))

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def extract_protein(fileReport, PATH_FASTA_DETECTED):

	"""
	This function is used to select the sequence identified by MacSyFinder and
	create a fasta with these sequences.

	:param fileReport: the file .report of the MacSyFinder analysis
	:type: string
	:param PATH_FASTA_DETECTED: absolute path of the folder where the rename will
	be write
	:type: string
	:return: the list of the sequence ids of the hit find by MacSyFinder, the
	list of all the new name for each sequences and the reference systems for
	each sequence.
	:rtype: list of string, list of string, list of string
	"""

	print "\n#################"
	print "# Protein extraction"
	print "#################\n"

	report_table = np.loadtxt(fileReport, dtype='string')
	number_prot = report_table.shape[0]
	index_remove = []

	for index in xrange(number_prot):
		if report_table[index][4] in PROTEIN_FUNCTION :
			system_number = report_table[index][7].split('_')[-1]

			if int(system_number) > 1 :
				report_table[index][1] = "_".join(report_table[index][1].split('_')[:-1])+"_"+system_number
			else:
				report_table[index][1] = "_".join(report_table[index][1].split('_')[:-1])
		else :
			index_remove.append(index)

	np.savetxt(os.path.join(PATH_FASTA_DETECTED, "remove_seq.seq"), report_table[np.array(index_remove),:], delimiter="\t", fmt="%s")
	report_table = np.delete(report_table, index_remove, axis=0)
	number_remove_protein = len(index_remove)

	print "There are %i proteins remove during this operation because they are not in the dictionnary" %number_remove_protein

	#NC_XXXXXX[_numero de systeme si deux systemes trouvés]_nomSysteme_D_nomProteine
	new_name = [report_table[i][1]+'_'+report_table[i][6]+'_D_'+"_".join(report_table[i][4].split('_')[1:]) for i in xrange(report_table.shape[0])]

	return report_table[:,0].tolist(), new_name, report_table[:,4].tolist()

##########################################################################################
##########################################################################################


def find_in_fasta(fileFasta, fileReport, listOfFile, PATH_FASTA_DETECTED):

	"""
	This function is used to create the fasta with MacSyFinder hits found

	:param fileFasta: name of the fasta database used in the MacSyfinder analysis
	:type: string
	:param fileReport: name of the file .report of the MacSyFinder analysis
	:type: string
	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of string
	:param PATH_FASTA_DETECTED: absolute path of the folder where the rename will
	be write
	:type: string
	:return: Nothing
	"""

	list_handle=[open(my_file,"w") for my_file in listOfFile]

	wanted, name_genes, keys_genes = extract_protein(fileReport, PATH_FASTA_DETECTED)
	seqiter = SeqIO.parse(open(fileFasta), 'fasta')

	print "\n#################"
	print "# Writing ..."
	print "#################\n"

	progression=1

	for seq in seqiter :
	  if seq.id in wanted:
			sys.stdout.write("%.2f% : %i/%i sequences wanted found\r" %(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()
			progression += 1

			index = wanted.index(seq.id)
			seq.description = ''
			seq.name = name_genes[index]
			seq.id = seq.name

			if keys_genes[index] in PROTEIN_FUNCTION :
				writing_file = re.search('[a-zA-Z0-9/_]+'+PROTEIN_FUNCTION[keys_genes[index]]+'\.fasta', "\t".join(listOfFile)).group(0)

				SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
			else :
				sys.exit("ERROR:: Function not known : "+keys_genes[index])

	print "Done!"

	#Close all file
	for open_file in list_handle:
		open_file.close()

	print "\n#################"
	print "# File wrote"
	print "#################\n"


##########################################################################################
##########################################################################################

"""
Function in R used to return the table of all the header similar in a same file

:param vect_file: list of all the file name
:type: robject.StrVector
:return: list of seq similiar in the same file
:rtype: robject.StrVector
"""

find_rename_fasta = robjects.r( ''' find_rename_fasta <- function(vect_file) {
final_vec=c()
for ( file in vect_file) {
	lines = readLines(file)
	newline=c()
	for (i in 1:length(lines)) {
		if (">" == strsplit(lines[i], "")[[1]][1]) {
			newline = c(newline, lines[i])}
		}
		final_vec = c(final_vec, names(table(newline)[table(newline) != 1]))
	}
	return(final_vec) }
''')


##########################################################################################
##########################################################################################

def rename_name_gene(listOfFile, PATH_FASTA_RENAME) :

	"""
	Function use to rename the sequence IDs of sequence with the same name in the same file

	:param new_listOfFile: list of all the file where we need to check if there are a same sequence
	id twice (or more) in the same file with absolute paths
	:type: list
	:param PATH_FASTA_RENAME: absolute path of the folder where the rename will
	be write
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# Rename protein"
	print "#################\n"

	new_listOfFile=[]

	for my_file in listOfFile :
		if os.stat(my_file).st_size != 0 :
			new_listOfFile.append(my_file)

	seq_to_rename = find_rename_fasta(new_listOfFile)
	dict_count = dict([(sequence[1:].rstrip(" "), 0) for sequence in seq_to_rename])
	progression=1

	create_folder(PATH_FASTA_RENAME)

	for my_file in new_listOfFile :

		file_name = os.path.basename(my_file)

		sys.stdout.write("%.2f% : %i/%i files renamed\r" %(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
		sys.stdout.flush()
		progression += 1

		handle = open(os.path.join(PATH_FASTA_RENAME, file_name), 'w')
		fasta_reading = SeqIO.parse(my_file, "fasta")

		for seq in fasta_reading :
			if seq.id in dict_count :
				if dict_count[seq.id] == 0 :
					dict_count[seq.id] += 1
				else :
					dict_count[seq.id] += 1
					if "NC_" in seq.id :
						seq.id = "_".join(seq.id.split("_")[:2])+"_Num"+str(dict_count[seq.id])+"_"+"_".join(seq.id.split("_")[2:])

					else :
						seq.id = seq.id.split("_")[0]+"_Num"+str(dict_count[seq.id])+"_"+"_".join(seq.id.split("_")[1:])
					seq.name = seq.id
					seq.description = ""

			SeqIO.write(seq, handle, "fasta")

		handle.close()

	print "Done!"
	return

##########################################################################################
##########################################################################################

def create_verified_fasta(listOfFile):

	"""
	Function used to extract the verified sequences for ATPase, prepilin peptidase, pilin (major and minor), IM platform

	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of string
	:return: Nothing
	"""

	PATH_TO_EXTRACT_SYSTEMS="/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/"

	print "\n#################"
	print "# Verified Fasta"
	print "#################\n"

	list_handle = [open(my_file, 'w') for my_file in listOfFile]

	info_extract = np.loadtxt(os.path.join(PATH_TO_EXTRACT_SYSTEMS,"all_seq_to_extract_annot.dat"), dtype="string")

	progression=1

	seqiter = SeqIO.parse(os.path.join(PATH_TO_EXTRACT_SYSTEMS, "all_genes_all_systems.fasta"), "fasta")

	for seq in seqiter :
		if seq.id in info_extract[:,0] :

			sys.stdout.write("%.2f% : %i/%i sequences wanted found\r" %(progression/float(info_tab.shape[0])*100, progression,info_tab.shape[0]))
			sys.stdout.flush()
			progression += 1

			position = info_extract[:,0].tolist().index(seq.id)

			if info_extract[position][1].split("_")[0] in ['T2SS','T4P', 'Tad']:

				if info_extract[position][1] in PROTEIN_FUNCTION :
					writing_file = re.search('[a-zA-Z0-9/_]+'+PROTEIN_FUNCTION[info_extract[position][1]]+'\.fasta', "\t".join(listOfFile)).group(0)

					seq.name = info_extract[position][3]+"_V_"+"_".join(info_extract[position][1].split("_")[1:])
					seq.id = seq.name
					seq.description = ''

					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")

			else :
				new_name = info_extract[position][2]+"_"+info_extract[position][1]

				if new_name in PROTEIN_FUNCTION :
					writing_file = re.search('[/a-zA-Z0-9_]*'+PROTEIN_FUNCTION[new_name]+'\.fasta', "\t".join(listOfFile)).group(0)

					seq.name = info_extract[position][3]+"_V_"+info_extract[position][1]
					seq.id = seq.name
					seq.description = ''

					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")

	print "Done!"
	return

##########################################################################################
##########################################################################################

def cut_seq_fasta_file(listOfFasta, file_cutoff=None, user_folder=None) :

	"""
	Function used to remove some sequence in the concatenated fasta, that after futher
	analysis, are considered as not good.

	:param listOfFasta: List of all the file fasta where we need to remove sequences
	:type: list of string
	:return: Nothing
	"""

	if file_cutoff is not None :
		DICT_CUTOFF=set_dict_cutoff(cutoff_file)
	else :
		DICT_CUTOFF=set_dict_cutoff_init(listOfFasta, user_folder)

	print "\n#################"
	print "# Cut concatetaned file"
	print "#################\n"

	dirname = os.path.dirname(listOfFasta[0])
	create_folder(new_path)

	for file in listOfFasta :
		current_file = os.path.basename(file)
		if current_file in DICT_CUTOFF:
			with open(os.path.join(dirname, "cut_off", current_file), "w") as writing_file :
				seqiter = SeqIO.parse(open(file), 'fasta')
				for seq in seqiter :
					if len(seq) < DICT_CUTOFF[current_file][1] and  len(seq) > DICT_CUTOFF[current_file][0] :
						SeqIO.write(seq, writing_file,"fasta")
	return
