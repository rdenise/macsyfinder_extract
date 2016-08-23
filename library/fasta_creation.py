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
##								Functions
##
##########################################################################################
##########################################################################################

def extract_protein(fileReport, INFO, PROTEIN_FUNCTION):

	"""
	This function is used to select the sequence identified by MacSyFinder and
	create a fasta with these sequences.

	:param fileReport: the file .report of the MacSyFinder analysis
	:type: string
	:param INFO: absolute path of the info_folder
	be write
	:type: string
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
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

	np.savetxt(os.path.join(INFO, "remove_report.seq"), report_table[np.array(index_remove),:], delimiter="\t", fmt="%s")
	report_table = np.delete(report_table, index_remove, axis=0)
	number_remove_protein = len(index_remove)

	print "There are %i proteins remove during this operation because they are not in the dictionnary" %number_remove_protein

	# NOTE NC_XXXXXX[_numero de systeme si deux systemes trouvés]_nomSysteme_D_nomProteine
	new_name = [report_table[i][1]+'_'+report_table[i][6]+'_D_'+"_".join(report_table[i][4].split('_')[1:]) for i in xrange(report_table.shape[0])]

	return report_table[:,0].tolist(), new_name, report_table[:,4].tolist()

##########################################################################################
##########################################################################################


def find_in_fasta(fileFasta, fileReport, listOfFile, INFO, PROTEIN_FUNCTION):

	"""
	This function is used to create the fasta with MacSyFinder hits found

	:param fileFasta: name of the fasta database used in the MacSyfinder analysis
	:type: string
	:param fileReport: name of the file .report of the MacSyFinder analysis
	:type: string
	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of string
	:param INFO: absolute path of the info_folder
	:type: string
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:return: Nothing
	"""

	list_handle=[open(my_file,"w") for my_file in listOfFile]

	wanted, name_genes, keys_genes = extract_protein(fileReport, INFO, PROTEIN_FUNCTION)
	seqiter = SeqIO.parse(open(fileFasta), 'fasta')

	print "\n#################"
	print "# Writing ..."
	print "#################\n"

	progression=1
	seq_wanted = len(wanted)

	for seq in seqiter :
	  if seq.id in wanted:
			sys.stdout.write("%.2f% : %i/%i sequences wanted found\r" %(progression/float(seq_wanted)*100, progression,seq_wanted))
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

	create_folder(PATH_FASTA_RENAME)

	new_listOfFile=[]

	for my_file in listOfFile :
		if os.stat(my_file).st_size != 0 :
			new_listOfFile.append(my_file)

	seq_to_rename = find_rename_fasta(new_listOfFile)
	dict_count = dict([(sequence[1:].rstrip(" "), 0) for sequence in seq_to_rename])
	progression=1
	number_of_file = len(new_listOfFile)

	for my_file in new_listOfFile :

		file_name = os.path.basename(my_file)

		sys.stdout.write("%.2f% : %i/%i files renamed\r" %(progression/float(number_of_file)*100, progression,number_of_file))
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
						# NOTE New name : NC_XXXXXX[_numero de systeme si deux systemes trouvés][_Num(et le nombre de fois nom trouvé)]_nomSysteme_D_nomProteine
						seq.id = "_".join(seq.id.split("_")[:2])+"_Num"+str(dict_count[seq.id])+"_"+"_".join(seq.id.split("_")[2:])

					else :
						# NOTE New name : NNNN[_numero de systeme si deux systemes trouvés][_Num(et le nombre de fois nom trouvé)]_nomSysteme_V_nomProteine
						seq.id = seq.id.split("_")[0]+"_Num"+str(dict_count[seq.id])+"_"+"_".join(seq.id.split("_")[1:])
					seq.name = seq.id
					seq.description = ""

			SeqIO.write(seq, handle, "fasta")

		handle.close()

	print "Done!"
	return

##########################################################################################
##########################################################################################

def create_verified_fasta(listOfFile, PROTEIN_FUNCTION):

	"""
	Function used to extract the verified sequences for ATPase, prepilin peptidase, pilin (major and minor), IM platform

	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of string
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
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

			sys.stdout.write("%.2f% : %i/%i sequences wanted found\r" %(progression/float(info_extract.shape[0])*100, progression,info_extract.shape[0]))
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

def write_remove_cutoff(dict_remove, INFO_folder):

	"""
	This function write a file in the information folder in markdown format to
	know which systems are removed and which sequence is in fault.

	:param dict_remove: dictionnary create by the function cut_seq_fasta_file()
	:type: dict
	:param INFO_folder: the absolute path to the info folder
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# Info concatetaned file"
	print "#################\n"

	info_cutoff_file = os.path.join(INFO_folder, "remove_cutoff.md")
	tmp = os.path.join(INFO_folder, "tmp.md")

	with open(tmp, "w") as w_file :
		for remove_system in dict_remove :
			w_file.write("# Systems remove by cutoff\n")
			if "_V_" in dict_remove[remove_system][0] :
				w_file.write("## _System %s :::: remove because %s is to long_\n" %(remove_system, dict_remove[remove_system][0]))
			else :
				w_file.write("## System %s :::: remove because %s is to long\n" %(remove_system, dict_remove[remove_system][0]))
			SeqIO.write(dict_remove[remove_system][1])
			w_file.write("# –––––––––––––––––––––––––––––––––––––––––––––––\n")
			w_file.write("# –––––––––––––––––––––––––––––––––––––––––––––––\n")
	with open(tmp, "r") as r_file:
		with open(info_cutoff_file, "w") as w_file:
			for line in r_file:
				w_file.write(line.rstrip()+"  \n")
	os.remove(tmp)
	return


##########################################################################################
##########################################################################################

def cut_seq_fasta_file(listOfFasta, PATH_FASTA_CUTOFF, INFO_folder, file_cutoff=None) :

	"""
	Function used to remove some sequence in the concatenated fasta, that after futher
	analysis, are considered as not good.

	:param listOfFasta: List of all the file fasta where we need to remove sequences
	:type: list of string
	:param PATH_FASTA_CUTOFF: path to the cutoff folder
	:type: string
    :param INFO_folder: the absolute path to the info folder
    :type: string
	:param file_cutoff: Name of the tabular file with the information for the cutoff if exist
	:type: string
	:return: Nothing
	"""

	if file_cutoff is not None :
		DICT_CUTOFF=set_dict_cutoff(cutoff_file)
	else :
		DICT_CUTOFF=set_dict_cutoff_init(listOfFasta, INFO_folder)

	print "\n#################"
	print "# Cutoff file"
	print "#################\n"

	create_folder(PATH_FASTA_CUTOFF)

	dict_remove = {}

	print "\n------------------------------------------
	print "| First read : Creation of the dictionnary"
	print "------------------------------------------\n"

	for my_file in listOfFasta :
		current_file = os.path.basename(my_file)
		if current_file in DICT_CUTOFF:

			seqiter = SeqIO.parse(my_file, 'fasta')
			number_seq = len(seqiter)
			progression = 1

			seqiter = SeqIO.parse(my_file, 'fasta')

			for seq in seqiter :
				sys.stdout.write("File : %s -> %.2f% : %i/%i sequences read\r" %(fasta_file, progression/float(number_seq)*100, progression, number_seq))
				sys.stdout.flush()
				progression += 1

				id_seq=seq.id.split("_")

				if "_D_" in seq.id :
					id_seq="_".join(id_seq[id_seq.index("D"):])
				else :
					id_seq="_".join(id_seq[id_seq_verif.index("V"):])

				if id_seq in dict_remove :
					continue
				elif len(seq) > DICT_CUTOFF[current_file][1] and  len(seq) < DICT_CUTOFF[current_file][0] :
					dict_remove[id_seq]=[seq.id,[]]
	print "Done!"

	print "\n-----------------------------"
	print "| Second read : Writing files"
	print "-----------------------------\n"

	for my_file in listOfFasta :
		current_file = os.path.basename(my_file)
		with open(os.path.join(PATH_FASTA_CUTOFF, current_file), "w") as writing_file :

			seqiter = SeqIO.parse(my_file, 'fasta')
			number_seq = len(seqiter)
			progression = 1

			seqiter = SeqIO.parse(my_file, 'fasta')
			for seq in seqiter :
				sys.stdout.write("File : %s -> %.2f% : %i/%i sequences read\r" %(fasta_file, progression/float(number_seq)*100, progression, number_seq))
				sys.stdout.flush()
				progression += 1

				id_seq=seq.id.split("_")

				if id_seq in dict_remove :
					dict_remove[id_seq][1].append(seq)
				else :
					SeqIO.write(seq, writing_file,"fasta")

	print "Done!"

	write_remove_cutoff(dict_remove, INFO_folder)

	return

##########################################################################################
##########################################################################################

def write_remove_concatenate(dict_remove, INFO_folder):

	"""
	This function write a file in the information folder in markdown format to
	know which sequences are extract and which verified sequences are the same
	of this one.

	:param dict_remove: dictionnary create by the function concatenate_detected_verified()
	:type: dict
	:param INFO_folder: the absolute path to the info folder
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# Info concatetaned file"
	print "#################\n"

	info_concatenate_file = os.path.join(INFO_folder, "remove_concatenate.md")
	tmp = os.path.join(INFO_folder, "tmp.md")

	with open(tmp, "w") as w_file :
		for remove_system in dict_remove :
			w_file.write("# Systems remove in concatenation\n")
			w_file.write("## System %s :::: %s identical\n" %(remove_system, dict_remove[remove_system][0]))
			SeqIO.write(dict_remove[remove_system][1])
			w_file.write("# –––––––––––––––––––––––––––––––––––––––––––––––\n")
			w_file.write("# –––––––––––––––––––––––––––––––––––––––––––––––\n")
	with open(tmp, "r") as r_file:
		with open(info_concatenate_file, "w") as w_file:
			for line in r_file:
				w_file.write(line.rstrip()+"  \n")
	os.remove(tmp)

	return

##########################################################################################
##########################################################################################


def concatenate_detected_verified(fasta_name, PATH_FASTA_DETECTED, PATH_FASTA_VERIFIED, INFO_folder, PATH_FASTA_CONCATENATED):

	"""
	Function that concatenate the verified and detected file and remove detected sequences
	that are already in the verified file. It write a file in the information folder in
	markdown format to know which sequences are extract and which verified sequences are
	the same of this one.

	:param fasta_name: the name of all the fasta file create ([protein_function].fasta)
	:type: list of string
	:param PATH_FASTA_DETECTED: absolute path to detected fasta folder
	:type: string
	:param PATH_FASTA_VERIFIED: absolute path to verified fasta folder
	:type: string
	:param INFO_folder: the absolute path to the info folder
	:type: string
	:param PATH_FASTA_CONCATENATED: absolute path to concatenated fasta folder
	:type: string
	:return: Nothing
	"""

	print "\n#################"
	print "# Concatetaned file"
	print "#################\n"

	for fasta_file in fasta_name :
		verified_fasta=os.path.join(fasta_file, PATH_FASTA_VERIFIED)
		detected_fasta=os.path.join(fasta_file, PATH_FASTA_DETECTED)
		concatenated_fasta=os.path.join(fasta_file, PATH_FASTA_CONCATENATED)

		os.system('cat "%s" > "%s"') %(verified_fasta, concatenated_fasta)

		list_seq_verified = list(SeqIO.parse(verified_fasta, "fasta"))
		list_id_verified = [seq.id for seq in list_seq_verified]
		list_seq_verified = [seq.seq for seq in list_seq_verified]

		seq_parser = SeqIO.parse(detected_fasta, "fasta")
		number_seq = len(seq_parser)
		progression = 1

		# NOTE Dictionaire avec en clef l'id espèce/système et en value une liste
		# NOTE ["l'id espèce/système du verifié qui correspond", [liste des sequences ATPase, IM ...]]
		dict_remove = {}

		seq_parser = SeqIO.parse(detected_fasta, "fasta")

		# IDEA Il faut tester au moins une fois pour voir si lors de la concatenation, je ne me retrouve pas avec des systems ou je n'ai pas tous enlevé. Exemple l'ATPase de X n'est pas la même que celle de Y mais l'IMplatform l'ai si c'est le cas X est a enlevé aussi pour son ATPase
		# IDEA Si idea précédente vrai alors il faut faire des fichiers temporaires des sequences que l'on garde et concatener par "cat" à la fin le fichier temporaire et son homonyme en verifié.

		with open(concatenated_fasta, "w") as w_file :
			for seq in seq_parser :

				sys.stdout.write("File : %s -> %.2f% : %i/%i sequences detected read\r" %(fasta_file, progression/float(number_seq)*100, progression,number_seq))
				sys.stdout.flush()
				progression += 1

				id_seq=seq.id.split("_")
				id_seq="_".join(id_seq[id_seq.index("D"):])

				if id_seq in dict_remove :
					dict_remove[id_seq][1].append(seq)

				elif seq.seq in list_seq_verified :
					index=list_seq_verified.index(seq.seq)

					id_seq_verif = list_id_verified[index].split("_")
					id_seq_verif = "_".join(id_seq[id_seq_verif.index("V"):])

					dict_remove[id_seq]=[id_seq_verif,[seq]]

				else :
					SeqIO.write(seq, "fasta")
		print("File : %s -> Done!" % fasta_file)

	# NOTE Dict remove complete and all concatenate write
	write_remove_concatenate(dict_remove, INFO_folder)

	return
