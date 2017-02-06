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
import shutil
import pandas as pd
from multiprocessing import Process
import glob
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
	:type: str
	:param INFO: absolute path of the info_folder
	be write
	:type: str
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.
	:type: dict
	:return: the list of the sequence ids of the hit find by MacSyFinder, the
	list of all the new name for each sequences and the reference systems for
	each sequence.
	:rtype: list of str, list of str, list of str
	"""

	print("\n-----------------")
	print("# Protein extraction")
	print("-----------------\n")

	report_table = np.genfromtxt(fileReport, dtype=str)
	number_prot = report_table.shape[0]
	index_remove = []

	for index in range(number_prot):
		if report_table[index][4] in PROTEIN_FUNCTION :
			system_number = report_table[index][7].split('_')[-1]

			if int(system_number) > 1 :
				# NOTE Pour gembases < 2015
				#report_table[index][1] = "_".join(report_table[index][1].split('_')[:-1])+"_"+system_number
				# NOTE Pour gembases > 2016
				report_table[index][1] = ".".join(report_table[index][1].split('.')[:-1])+"_"+system_number
			else:
				# NOTE Pour gembases < 2015
				#report_table[index][1] = "_".join(report_table[index][1].split('_')[:-1])
				# NOTE Pour gembases > 2016
				report_table[index][1] = ".".join(report_table[index][1].split('.')[:-1])
		else :
			index_remove.append(index)

	np.savetxt(os.path.join(INFO, "remove_report.seq"), report_table[np.array(index_remove),:], delimiter="\t", fmt="%s")
	report_table = np.delete(report_table, index_remove, axis=0)
	number_remove_protein = len(index_remove)

	print("There are {} proteins remove during this operation because they are not in the dictionnary".format(number_remove_protein))

	# NOTE AAAAKKK.B.LLLLL[_numero de systeme si deux systemes trouvés]_nomSysteme_D_nomProteine (for 2016 gembase format)
	# NOTE NC_XXXXXX[_numero de systeme si deux systemes trouvés]_nomSysteme_D_nomProteine (for 2015 gembase format)
	# NOTE XXXX[_numero de systeme si deux systemes trouvés]_nomSysteme_D_nomProteine (for 2013 gembase format)

	# XXX Je crée une table de traduction entre mon nom et le nom générique de la séquence
	with open(os.path.join(INFO, "translation_table_detected.tab"), 'w') as translate_file :
		new_name_list = []
		translate_file.write("#Name_database\tNew_name\n")

		for i in range(report_table.shape[0]) :
			new_name = "{}_{}_D_{}".format(report_table[i][1], report_table[i][6], "_".join(report_table[i][4].split('_')[1:]))
			new_name_list.append(new_name)
			translate_file.write("{}\t{}\n".format(report_table[i][0], new_name))

	return report_table[:,0].tolist(), new_name_list, report_table[:,4].tolist()

##########################################################################################
##########################################################################################


def find_in_fasta(fileFasta, fileReport, listOfFile, INFO, PROTEIN_FUNCTION):

	"""
	This function is used to create the fasta with MacSyFinder hits found

	:param fileFasta: name of the fasta database used in the MacSyfinder analysis
	:type: str
	:param fileReport: name of the file .report of the MacSyFinder analysis
	:type: str
	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of str
	:param INFO: absolute path of the info_folder
	:type: str
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:return: Nothing
	"""

	list_handle = [open(my_file,"w") for my_file in listOfFile]

	wanted, name_genes, keys_genes = extract_protein(fileReport, INFO, PROTEIN_FUNCTION)
	seqiter = SeqIO.parse(open(fileFasta), 'fasta')

	print("\n-----------------")
	print("# Writing ...")
	print("-----------------\n")

	progression=1
	seq_wanted = len(wanted)

	for seq in seqiter :
		if seq.id in wanted:
			sys.stdout.write("{:.2f}% : {}/{} sequences wanted found\r".format(progression/float(seq_wanted)*100, progression,seq_wanted))
			sys.stdout.flush()
			progression += 1

			index = wanted.index(seq.id)
			seq.description = ''
			seq.name = name_genes[index]
			seq.id = seq.name

			if keys_genes[index] in PROTEIN_FUNCTION :
				writing_file = re.search('[a-zA-Z0-9/_]+{}\.fasta'.format(PROTEIN_FUNCTION[keys_genes[index]]), "\t".join(listOfFile)).group(0)

				SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
			else :
				sys.exit("ERROR:: Function not known : {}".format(keys_genes[index]))

	print()
	print("Done!")

	# XXX Close all file
	for open_file in list_handle:
		open_file.close()

	print("\n#################")
	print("# File wrote")
	print("#################\n")

	return

##########################################################################################
##########################################################################################
def write_fasta_multithreads(subfolderFasta, listOfFile, wanted, name_genes, keys_genes, PROTEIN_FUNCTION):

	'''
	:param subfolderFasta: name of the fasta database used in the MacSyfinder analysis split by species and by process
	:type: str
	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of str
	:param wanted: name of the wanted sequence to write
	:type: list
	:param name_genes: new name of the genes after extraction (the new name I gave)
	:type: list
	:param key_genes: name of the function of the proteins to know the annotation
	:type: list
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:return: Nothing
	'''

	list_handle = [open(my_file,"w") for my_file in listOfFile]

	for fileFasta in subfolderFasta :
		seqiter = SeqIO.parse(open(fileFasta), 'fasta')

		for seq in seqiter :
			if seq.id in wanted:

				index = wanted.index(seq.id)
				seq.description = ''
				seq.name = name_genes[index]
				seq.id = seq.name

				if keys_genes[index] in PROTEIN_FUNCTION :
					writing_file = re.search('[a-zA-Z0-9/_]+{}\.fasta'.format(PROTEIN_FUNCTION[keys_genes[index]]), "\t".join(listOfFile)).group(0)
					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")

				else :
					sys.exit("ERROR:: Function not known : {}".format(keys_genes[index]))


	# XXX Close all file
	for open_file in list_handle:
		open_file.close()

	return

##########################################################################################
##########################################################################################


def find_in_fasta_multithreads(folderFasta, fileReport, listOfFile, INFO, PROTEIN_FUNCTION, nb_thread):

	"""
	This function is used to create the fasta with MacSyFinder hits found

	:param folderFasta: name of the fasta database used in the MacSyfinder analysis split by species
	:type: str
	:param fileReport: name of the file .report of the MacSyFinder analysis
	:type: str
	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of str
	:param INFO: absolute path of the info_folder
	:type: str
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:param nb_thread: number of thread choose by the user
	:type: int
	:return: Nothing
	"""

	nb_fasta = len(folderFasta)

	list_process = []

	# XXX Permet de ne pas avoir plus de threads que de fichier
	if nb_fasta//nb_thread == 0 :
		nb_thread = nb_fasta%nb_thread

	split_list_folder = [folderFasta[i:i+nb_fasta//nb_thread] for i in range(0, nb_fasta, nb_fasta//nb_thread)]

	wanted, name_genes, keys_genes = extract_protein(fileReport, INFO, PROTEIN_FUNCTION)

	print()
	print("Begin multi processing extraction ...")
	print()

	for i in range(nb_thread) :
		process_directory = os.path.join(os.path.dirname(listOfFile[0]), "tmp","process_{}".format(i))
		create_folder(process_directory)

		# XXX Je crée les fichiers fasta temporaire
		process_list_file = [os.path.join(process_directory, os.path.basename(myfile)) for myfile in listOfFile]

		# XXX Je crée et lance les process qui exécuteront la fonction d'écriture pour chaque lot de fasta
		list_process.append(Process(target=write_fasta_multithreads, args=(split_list_folder[i], process_list_file, wanted, name_genes, keys_genes, PROTEIN_FUNCTION)))
		list_process[-1].start()

		print("Begin process for {}".format(list_process[-1].pid))

	print()
	print("List of all process done : ")
	print()

	for process in list_process :
		# XXX J'attend après chaque process qu'il finisse
		process.join()
		print("Done for {} !".format(process.pid))

	for myfile in listOfFile :
		file_name = os.path.basename(myfile)
		list_tmp = glob.glob(os.path.join(os.path.dirname(listOfFile[0]), "tmp", "*", file_name))

		# XXX Je concatène toutes les séquences extraitent
		os.system("cat {} > {}".format(" ".join(list_tmp), myfile))

	# XXX Je supprime tout le dossier tmp
	shutil.rmtree(os.path.join(os.path.dirname(listOfFile[0]), "tmp"))

	print("\n#################")
	print("# File wrote")
	print("#################\n")

	return

##########################################################################################
##########################################################################################

def write_in_info(info_name, dict_info):

	"""
	Function use to rename the sequence IDs of sequence with the same name in the same file

	:param info_name: open file where it will write the information about the systems
	found in the fasta detected and verified (if exists)
	:type: file
	:param dict_info : The dictionnary that is create during the function rename_name_gene
	and that is like that : {name_species: {nameSystem_numero: {name_protein1: count, name_protein2: count} ...} ...}
	:type: dict
	:return: Nothing
	"""

	for species in dict_info:
		for key_system in dict_info[species] :
			system, numero = key_system.split("_")
			info = "{}\t{}\t{}\t{}\n".format(species, system, numero, " ".join(["{}:{}".format(key,value) for key, value in dict_info[species][key_system].items()]))
			info_name.write(info)

	return

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

def rename_name_gene(listOfFile, PATH_FASTA_RENAME, info_name, DICT_SYSTEMS, dict_info) :

	"""
	Function use to rename the sequence IDs of sequence with the same name in the same file

	:param new_listOfFile: list of all the file where we need to check if there are a same sequence
	id twice (or more) in the same file with absolute paths
	:type: list
	:param PATH_FASTA_RENAME: absolute path of the folder where the rename will
	be write
	:type: str
	:param info_name: open file where I will write the information about the systems
	found in the fasta detected and verified (if exists)
	:type: file
	:param DICT_SYSTEMS: The dictionnary that contains the name of all the
	systems in key and the list of all the protein of this systems in keys
	:type: dict
	:param dict_info: dictionary we want to fill with the information about the systems found
	:type: dict
	:return: A dictionnary that contains {name_species: {nameSystem_numero: {name_protein1: count, name_protein2: count} ...} ...}
	and a list of the systems found
	:rtype: dict, list of str
	"""

	print("\n#################")
	print("# Rename protein")
	print("#################\n")

	create_folder(PATH_FASTA_RENAME)

	list_system = []
	generic=False
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

		sys.stdout.write("{:.2f}% : {}/{} files renamed\r".format(progression/float(number_of_file)*100, progression,number_of_file))
		sys.stdout.flush()
		progression += 1

		handle = open(os.path.join(PATH_FASTA_RENAME, file_name), 'w')
		fasta_reading = SeqIO.parse(my_file, "fasta")

		for seq in fasta_reading :
			seq_id_split = seq.id.split("_")

			# NOTE Ici obligé de mettre le index_system_name avant et de reteste entre D et V car sinon pour les sequence unique ça n'existe pas
			if "_D_" in seq.id :
				index_system_name = seq_id_split.index("D")-1
			elif "_V_" in seq.id :
				index_system_name = seq_id_split.index("V")-1
			else :
				sys.exit("ERROR::Wrong seqID : {}".format(seq.id))

			if seq.id in dict_count :
				if dict_count[seq.id] == 0 :
					dict_count[seq.id] += 1
				else :
					dict_count[seq.id] += 1

					# NOTE New name : NC_XXXXXX[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_D_nomProteine
					# NOTE New name : NNNN[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_V_nomProteine
					# NOTE New name : AAAAKKK.B.LLLLL[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_D_nomProteine
					seq.id = "{}_Num{}_{}".format("_".join(seq_id_split[:index_system_name]), str(dict_count[seq.id]), "_".join(seq_id_split[index_system_name:]))

					seq.name = seq.id
					seq.description = ""

			SeqIO.write(seq, handle, "fasta")

			if re.search("_[0-9]_", seq.id) :
				numero = re.search("_[0-9]_", seq.id).group(0).replace("_", "")
			else :
				numero = "."

			if "NC_" in seq.id :
				species_name = "_".join(seq_id_split[:2])
			else :
				species_name = seq_id_split[0]

			system = seq_id_split[index_system_name]

			if system in ("Tcp", "R64", "Cof", "Bfp", "MSH", "Lng"):
				system = "T4bP"

			if system not in list_system:
				if system.lower() in ["generique", "generic"] :
					generic=system
				else :
					list_system.append(system)

			key_system = "{}_{}".format(system, numero)
			protein_name = "_".join(seq_id_split[index_system_name+2:])

			# NOTE On ajoute ici que si je ne suis pas dans les clés ou je suis dans les clés mais pas la liste de valeurs (ESCO peut avoir T2SS et T4P)
			if species_name not in dict_info :
				dict_info[species_name]={}

			if key_system not in dict_info[species_name] :
				if system.lower() in ["generique", "generic"] or "_V_" in seq.id :
					dict_info[species_name][key_system]={protein_name:0}
				else :
					dict_info[species_name][key_system] = {protein:0 for protein in DICT_SYSTEMS[system]}

			try :
				dict_info[species_name][key_system][protein_name] += 1
			except KeyError :
				dict_info[species_name][key_system][protein_name] = 1

		handle.close()

	write_in_info(info_name, dict_info)

	if generic :
		list_system = sorted(list_system)
		list_system.append(generic)

	print()
	print("Done!")

	return dict_info, list_system

##########################################################################################
##########################################################################################

def create_verified_fasta(listOfFile, PROTEIN_FUNCTION, data_fasta, info_dat, INFO):

	# NOTE Pas du tous adapté si on a un programme général si des verifiés autre que les miens c'est foutu

	"""
	Function used to extract the verified sequences for ATPase, prepilin peptidase, pilin (major and minor), IM platform

	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of str
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:data_fasta: Fasta file with the verified sequence of the systems
	:type: str
	:info_dat:File with the information about the verified systems with this information : #SeqID Gene System SystID
	:type: str
	:param INFO: absolute path of the info_folder
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# Verified Fasta")
	print("#################\n")

	list_handle = [open(my_file, 'w') for my_file in listOfFile]

	w_file = open(os.path.join(INFO, "translation_table_verified.tab"), 'w')
	w_file.write("#Name_fasta_file\tNew_name\n")

	info_extract = np.genfromtxt(info_dat, dtype=str, delimiter="\t")

	progression=1

	seqiter = SeqIO.parse(data_fasta, "fasta")

	for seq in seqiter :
		if seq.id in info_extract[:,0] :

			sys.stdout.write("{:.2f}% : {}/{} sequences wanted found\r".format(progression/float(info_extract.shape[0])*100, progression,info_extract.shape[0]))
			sys.stdout.flush()
			progression += 1

			position = info_extract[:,0].tolist().index(seq.id)

			# IDEA pour faire un truc général je pourrais par exemple crée la liste à partir du fichier de definition de protein_function.def

			if info_extract[position][1].split("_")[0] in ['T2SS','T4P', 'Tad']:

				if info_extract[position][1] in PROTEIN_FUNCTION :
					writing_file = re.search('[a-zA-Z0-9/_]+'+PROTEIN_FUNCTION[info_extract[position][1]]+'\.fasta', "\t".join(listOfFile)).group(0)

					seq.name = "{}_{}_V_{}".format(info_extract[position][-1], info_extract[position][3].split("_")[-1], "_".join(info_extract[position][1].split("_")[1:]))

					# XXX je fais aussi un fichier translation pour les séquences vérifiées
					w_file.write("{}\t{}\n".format(seq.id, seq.name))

					seq.id = seq.name
					seq.description = ''

					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
			else :
				# NOTE Permet d'avoir le bon nom si dans la colonne 2 j'ai que gspD par exemple et comme ça je reforme T2SS_gspD pour les sytèmes Com T4bP ...

				new_name = info_extract[position][2]+"_"+info_extract[position][1]

				if new_name in PROTEIN_FUNCTION :
					writing_file = re.search('[/a-zA-Z0-9_]*'+PROTEIN_FUNCTION[new_name]+'\.fasta', "\t".join(listOfFile)).group(0)

					seq.name = "{}_{}_V_{}".format(info_extract[position][-1], info_extract[position][2], info_extract[position][1])
					seq.id = seq.name
					seq.description = ''

					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")

	print()
	print("Done!")

	for open_file in list_handle :
		open_file.close()

	# XXX On ferme le fichier de translation_table
	w_file.close()

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
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# Write info cutoff file")
	print("#################\n")

	info_cutoff_file = os.path.join(INFO_folder, "remove_cutoff.md")
	tmp = os.path.join(INFO_folder, "tmp.md")

	with open(tmp, "w") as w_file :
		w_file.write("# Systems removed by cutoff\n")
		for remove_system in dict_remove :
			w_file.write("___\n")
			if "_V_" in dict_remove[remove_system][0] :
				w_file.write("#### *System {} :::: Removed*\n".format(remove_system))
				w_file.write("##### *Because {} is too {}*\n".format(dict_remove[remove_system][0],dict_remove[remove_system][2]))
			else :
				w_file.write("#### System {} :::: Removed\n".format(remove_system))
				w_file.write("##### Because {} is too {}\n".format(dict_remove[remove_system][0],dict_remove[remove_system][2]))
			w_file.write("___\n")
			SeqIO.write(dict_remove[remove_system][1], w_file, "fasta")

	with open(tmp, "r") as r_file:
		with open(info_cutoff_file, "w") as w_file:
			for line in r_file:
				w_file.write(line.rstrip().replace(">", "\>")+"  \n")
	os.remove(tmp)

	print("Done!")

	return


##########################################################################################
##########################################################################################

def cut_seq_fasta_file(listOfFasta, PATH_FASTA_CUTOFF, INFO_folder, file_cutoff=None) :

	"""
	Function used to remove some sequence in the concatenated fasta, that after futher
	analysis, are considered as not good.

	:param listOfFasta: List of all the file fasta where we need to remove sequences
	:type: list of str
	:param PATH_FASTA_CUTOFF: path to the cutoff folder
	:type: str
    :param INFO_folder: the absolute path to the info folder
    :type: str
	:param file_cutoff: Name of the tabular file with the information for the cutoff if exist or True if not
	:type: str
	:return: Nothing
	"""

	if file_cutoff == True :
		DICT_CUTOFF=set_dict_cutoff_init(listOfFasta, INFO_folder)
	else :
		DICT_CUTOFF=set_dict_cutoff(cutoff_file)


	print("\n#################")
	print("# Cutoff file")
	print("#################\n")

	create_folder(PATH_FASTA_CUTOFF)

	dict_remove = {}

	print("\n------------------------------------------")
	print("| First read : Creation of the dictionnary")
	print("------------------------------------------\n")

	for my_file in listOfFasta :
		current_file = os.path.basename(my_file)
		if current_file in DICT_CUTOFF:

			seqiter = SeqIO.parse(my_file, 'fasta')
			number_seq = len(list(seqiter))
			progression = 1

			seqiter = SeqIO.parse(my_file, 'fasta')

			for seq in seqiter :
				sys.stdout.write("File : {} -> {:.2f}% : {}/{} sequences read\r".format(current_file, progression/float(number_seq)*100, progression, number_seq))
				sys.stdout.flush()
				progression += 1

				id_seq=seq.id.split("_")

				if "_D_" in seq.id :
					id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("D")]))
				else :
					id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("V")]))

				if id_seq in dict_remove :
					continue
				elif len(seq) > DICT_CUTOFF[current_file][1] or  len(seq) < DICT_CUTOFF[current_file][0] :
					if len(seq) > DICT_CUTOFF[current_file][1] :
						dict_remove[id_seq]=[seq.id,[], "long"]
					else :
						dict_remove[id_seq]=[seq.id,[], "short"]
		print()
		print("File : {} -> Done!".format(current_file))

	print("\n-----------------------------")
	print("| Second read : Writing files")
	print("-----------------------------\n")

	for my_file in listOfFasta :
		current_file = os.path.basename(my_file)
		with open(os.path.join(PATH_FASTA_CUTOFF, current_file), "w") as writing_file :

			seqiter = SeqIO.parse(my_file, 'fasta')
			number_seq = len(list(seqiter))
			progression = 1

			seqiter = SeqIO.parse(my_file, 'fasta')
			for seq in seqiter :
				sys.stdout.write("File : {} -> {:.2f}% : {}/{} sequences read\r".format(current_file, progression/float(number_seq)*100, progression, number_seq))
				sys.stdout.flush()
				progression += 1

				id_seq=seq.id.split("_")

				if "_D_" in seq.id :
					id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("D")]))
				else :
					id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("V")]))

				if id_seq in dict_remove :
					dict_remove[id_seq][1].append(seq)
				else :
					SeqIO.write(seq, writing_file,"fasta")

		print()
		print("File : {} -> Done!".format(current_file))

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
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# Write info concatetaned file")
	print("#################\n")

	info_concatenate_file = os.path.join(INFO_folder, "remove_concatenate.md")
	tmp = os.path.join(INFO_folder, "tmp.md")

	with open(tmp, "w") as w_file :
		w_file.write("# Systems remove in concatenation\n")
		for remove_system in dict_remove :
			w_file.write("___\n")
			w_file.write("#### System {} :::: {} identical\n".format(remove_system, dict_remove[remove_system][0]))
			w_file.write("##### Because of {}\n".format(dict_remove[remove_system][2]))
			w_file.write("___\n")
			SeqIO.write(dict_remove[remove_system][1], w_file, "fasta")
	with open(tmp, "r") as r_file:
		with open(info_concatenate_file, "w") as w_file:
			for line in r_file:
				w_file.write(line.rstrip().replace(">", "\>")+"  \n")
	os.remove(tmp)

	print("Done!")

	return

##########################################################################################
##########################################################################################


def concatenate_detected_verified(fasta_name, PATH_FASTA_DETECTED, PATH_FASTA_VERIFIED, INFO_folder, PATH_FASTA_CONCATENATED, PATH_FASTA_DETECTED_SELECTED):

	"""
	Function that concatenate the verified and detected file and remove detected sequences
	that are already in the verified file. It write a file in the information folder in
	markdown format to know which sequences are extract and which verified sequences are
	the same of this one.

	:param fasta_name: the name of all the fasta file create ([protein_function].fasta)
	:type: list of str
	:param PATH_FASTA_DETECTED: absolute path to detected fasta folder
	:type: str
	:param PATH_FASTA_VERIFIED: absolute path to verified fasta folder
	:type: str
	:param INFO_folder: the absolute path to the info folder
	:type: str
	:param PATH_FASTA_CONCATENATED: absolute path to concatenated fasta folder
	:type: str
	:param PATH_FASTA_DETECTED_SELECTED: absolute path to selected fasta folder
	:type: str
	:return: Nothing
	"""

	print("\n#################")
	print("# Concatetaned file")
	print("#################\n")

	# NOTE Dictionaire avec en clef l'id espèce/système et en value une liste
	# NOTE ["l'id espèce/système du verifié qui correspond", [liste des sequences ATPase, IM ...]]
	dict_remove = {}

	print("\n------------------------------------------")
	print("| First read : Creation of the dictionnary")
	print("------------------------------------------\n")

	for fasta_file in fasta_name :
		verified_fasta=os.path.join(PATH_FASTA_VERIFIED, fasta_file)
		detected_fasta=os.path.join(PATH_FASTA_DETECTED, fasta_file)
		concatenated_fasta=os.path.join(PATH_FASTA_CONCATENATED, fasta_file)

		list_seq_verified = list(SeqIO.parse(verified_fasta, "fasta"))
		list_id_verified = [seq.id for seq in list_seq_verified]
		list_seq_verified = [seq.seq for seq in list_seq_verified]

		list_seq_detected_OK = []
		list_id_detected_OK = []

		seq_parser = SeqIO.parse(detected_fasta, "fasta")
		number_seq = len(list(seq_parser))
		progression = 1

		seq_parser = SeqIO.parse(detected_fasta, "fasta")

		# NOTE Il y avait un problème : le nom/id de l'epèce + système ne doit pas contenir le _NumX_ car ce Num fait référence au nombre de duplicat de la protéine (exemple deux ATPase gspE)
		# NOTE Quelques systèmes on des sequences qui sont similaire pour toutes les protéines sauf une exemple ESCO3 et NC_011993 qui sont identique pour tous sauf ATPase (98% seulement)

		for seq in seq_parser :

			sys.stdout.write("File : {} -> {:.2f}% : {}/{} sequences detected read\r".format(fasta_file, progression/float(number_seq)*100, progression,number_seq))
			sys.stdout.flush()
			progression += 1

			id_seq=seq.id.split("_")
			id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("D")]))

			if id_seq in dict_remove :
				continue

			elif seq.seq in list_seq_verified :
				index=list_seq_verified.index(seq.seq)

				id_seq_verif = list_id_verified[index].split("_")
				id_seq_verif = re.sub("Num[0-9]_", "", "_".join(id_seq_verif[:id_seq_verif.index("V")]))

				# NOTE dans le dictionnaire je met le système vérifié en premier, toutes les séquences du système identitique en deuxième et la séquence qui en est la cause en troisème
				dict_remove[id_seq]=[id_seq_verif,[], seq.id]

			elif seq.seq in list_seq_detected_OK :
				index=list_seq_detected_OK.index(seq.seq)

				id_seq_reference = list_id_detected_OK[index]

				# NOTE dans le dictionnaire je met le système référence en premier, toutes les séquences du système identitique en deuxième et la séquence qui en est la cause en troisème
				dict_remove[id_seq]=[id_seq_reference,[], seq.id]

			else :
				list_seq_detected_OK.append(seq.seq)
				list_id_detected_OK.append(id_seq)

		print()
		print("File : {} -> Done!".format(fasta_file))

	print("\n-----------------------------")
	print("| Second read : Writing files")
	print("-----------------------------\n")

	for fasta_file in fasta_name :
		verified_fasta=os.path.join(PATH_FASTA_VERIFIED, fasta_file)
		detected_fasta=os.path.join(PATH_FASTA_DETECTED, fasta_file)
		new_detected_fasta=os.path.join(PATH_FASTA_DETECTED_SELECTED, fasta_file)
		concatenated_fasta=os.path.join(PATH_FASTA_CONCATENATED, fasta_file)

		os.system('cat "{}" > "{}"'.format(verified_fasta, concatenated_fasta))

		seq_parser = SeqIO.parse(detected_fasta, "fasta")
		number_seq = len(list(seq_parser))
		progression = 1

		seq_parser = SeqIO.parse(detected_fasta, "fasta")

		with open(concatenated_fasta, "a") as w_file :
			with open(new_detected_fasta, "w") as nd_file :
				for seq in seq_parser :

					sys.stdout.write("File : {} -> {:.2f}% : {}/{} sequences detected read\r".format(fasta_file, progression/float(number_seq)*100, progression,number_seq))
					sys.stdout.flush()
					progression += 1

					id_seq=seq.id.split("_")
					id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("D")]))

					if id_seq in dict_remove :
						dict_remove[id_seq][1].append(seq)

					else :
						SeqIO.write(seq, w_file, "fasta")
						SeqIO.write(seq, nd_file, "fasta")
		print()
		print("File : {} -> Done!".format(fasta_file))

	# NOTE Dict remove complete and all concatenate write
	write_remove_concatenate(dict_remove, INFO_folder)

	return


##########################################################################################
##########################################################################################


def rename_seq_translation_table(INFO_folder):

	# NOTE Vérifié si c'est correct avec le verifier et le detected en comparant les séquences

	"""
	Function that will change the name of the rename sequence (it use the fact
	that the sequence are changed by place of occurence and it's the same
	in the fasta file and in the translation_table).

	:param INFO_folder: the absolute path to the info folder
	:type: str
	"""

	all_translation_files = glob.glob(os.path.join(INFO_folder, "*tab"))

	dict_count = {}

	for translation_table in all_translation_files :
		with open(os.path.join(INFO_folder, "tmp"), "w") as tmp_file, open(translation_table, 'r') as r_file :
			for line in r_file :
				split_line = line.split()

				if "#" in line :
					tmp_file.write(line)

				elif split_line[-1] in dict_count :
					dict_count[split_line[-1]] += 1
					id_split = split_line[-1].split("_")

					# NOTE Ici obligé de mettre le index_system_name avant et de reteste entre D et V car sinon pour les sequence unique ça n'existe pas
					if "_D_" in seq.id :
						index_system_name = id_split.index("D")-1
					elif "_V_" in seq.id :
						index_system_name = id_split.index("V")-1

					new_name = "{}_Num{}_{}".format("_".join(id_split[:index_system_name]), str(dict_count[seq.id]), "_".join(id_split[index_system_name:]))
					tmp_file.write(line.replace(split_line[-1], new_name))

				else :
					tmp_file.write(line)
					dict_count[line.split()[-1]] = 1

		os.rename(os.path.join(INFO_folder, "tmp"), translation_table)

	return
