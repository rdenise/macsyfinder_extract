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
import glob
from multiprocessing import Pool
import subprocess
import shlex
from set_params import *

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################


def set_name(row, dict_count) :

	"""
	Funtion that chage the name of the columns Replicon_name and return the new name

	:param row: The row of a dataframe pandas
	:type: Pandas.Series
	:param dict_count: dictionary with the name of the count of all name seen
	:type: dict
	:return: The new name
	:rtype: str
	"""

	if (int(row.System_Id.split('_')[-1]) > 1) :
	    NewNameTmp = "{}_{}".format(row.Replicon_name, row.System_Id.split('_')[-1])
	else :
	    NewNameTmp =  row.Replicon_name

	NewName = "{}_{}_D_{}".format(NewNameTmp, row.Predicted_system, "_".join(row.Gene.split('_')[1:]))

	if NewName in dict_count :
		dict_count[NewName] += 1
		NewName = "{}_Num{}_{}_D_{}".format(NewNameTmp, dict_count[NewName], row.Predicted_system, "_".join(row.Gene.split('_')[1:]))
	else :
		dict_count[NewName] = 1

	return NewName

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
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.read_protein_function()
	:type: dict
	:return: A dataframe corresponding to the report with only the protein I want and with a new columns with the new name.
	:rtype: Pandas.Dataframe
	"""

	print("\n-----------------")
	print("# Protein extraction")
	print("-----------------\n")

	# XXX Je lis le fichier report et je lui donne le bon nom de header car certain fichier ne l'ont pas ou plus
	names_dataframe=['Hit_Id','Replicon_name','Position','Sequence_length','Gene','Reference_system','Predicted_system','System_Id','System_status','Gene_status','i-evalue','Score','Profile_coverage','Sequence_coverage','Begin_match','End_match']
	report_table = pd.read_table(fileReport, names=names_dataframe, dtype="str", comment="#")

	# XXX Je fais un sous dataframe qui contient la liste de toutes les lignes que je ne veux pas et je l'écris dans un fichier
	report_table[~report_table.Gene.isin(PROTEIN_FUNCTION)].reset_index(drop=True).to_csv(os.path.join(INFO, "remove_seq.report"), sep="\t", index=False, header=False)
	number_remove_protein = report_table[~report_table.Gene.isin(PROTEIN_FUNCTION)].shape[0]
	print("There are {} proteins remove during this operation because they are not in the dictionnary".format(number_remove_protein))

	# XXX Avec ceux que je veux, je crée une première fois la colonne avec les nouveaux noms
	new_report_table=report_table[report_table.Gene.isin(PROTEIN_FUNCTION)].reset_index(drop=True)

	# NOTE New name : AAAAKKK.B.LLLLL[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_D_nomProteine (for 2016 gembase format)
	# NOTE New name : NC_XXXXXX[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_D_nomProteine (for 2015 gembase format)
	# NOTE New name : NNNN[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_V_nomProteine (for 2013 gembase format)
	dict_count = {}
	new_report_table["NewName"] = new_report_table.apply(set_name, args=[dict_count], axis=1)

	# XXX Je crée une table de traduction entre mon nom et le nom générique de la séquence et je change les nouveaux nom pour le nom final
	new_report_table.loc[:,["Hit_Id","NewName"]].to_csv(os.path.join(INFO, "translation_table_detected.tab"),sep="\t", index=False)

	return new_report_table

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
	:return: The report table with the new name at the end
	:rtype: pandas.DataFrame
	"""

	list_handle = [open(my_file,"w") for my_file in listOfFile]

	report_table = extract_protein(fileReport, INFO, PROTEIN_FUNCTION)
	seqiter = SeqIO.parse(open(fileFasta), 'fasta')

	print("\n-----------------")
	print("# Writing ...")
	print("-----------------\n")

	progression=1
	seq_wanted = report_table.shape[0]
	report_table.set_index("Hit_Id", inplace=True)

	for seq in seqiter :
		if seq.id in report_table.index:
			sys.stdout.write("{:.2f}% : {}/{} sequences wanted found\r".format(progression/float(seq_wanted)*100, progression,seq_wanted))
			sys.stdout.flush()
			progression += 1

			seq.description = ''
			seq.id = report_table.loc[seq.id, "NewName"]
			seq.name = ''

			if report_table.loc[seq.id, "Gene"] in PROTEIN_FUNCTION :
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

	return report_table

##########################################################################################
##########################################################################################
def write_fasta_multithreads(fileFasta, listOfFile, wanted, PROTEIN_FUNCTION, number_total):

	'''
	:param fileFasta: name of the fasta file for one replicon
	:type: str
	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of str
	:param wanted: the report table with the information about the new name add at the end
	:type: pandas.Dataframe
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:return: Nothing
	'''

	Replicon = os.path.splitext(os.path.basename(fileFasta))[0]

	# XXX Ici je crée les nouveaux directory de chaque sous process (le nom du replicon)

	wanted_modif = wanted.set_index("Replicon_name")
	process_directory = os.path.join(os.path.dirname(listOfFile[0]), "tmp","{}".format(Replicon))
	create_folder(process_directory)
	writing_list = [os.path.join(process_directory, os.path.basename(my_file)) for my_file in listOfFile]
	list_handle = [open(os.path.join(process_directory, os.path.basename(my_file)),"w") for my_file in listOfFile]

	seqiter = SeqIO.to_dict(SeqIO.parse(open(fileFasta), 'fasta'))

	for index, line in wanted_modif.loc[[Replicon]].iterrows() :
		seqiter[line.Hit_Id].description = ''
		seqiter[line.Hit_Id].name = ""
		seqiter[line.Hit_Id].id = line.NewName

		if line.Gene in PROTEIN_FUNCTION :
			writing_file = re.search('[a-zA-Z0-9/_\.]+{}\.fasta'.format(PROTEIN_FUNCTION[line.Gene]), "\t".join(writing_list)).group(0)
			SeqIO.write(seqiter[line.Hit_Id], list_handle[writing_list.index(writing_file)], "fasta")

		else :
			sys.exit("ERROR:: Function not known : {}".format(line.Gene))


	# XXX Close all file
	for open_file in list_handle:
		open_file.close()

	actual_number_folder = len(glob.glob(os.path.join(os.path.dirname(listOfFile[0]), "tmp", "*")))
	print("{:.2f}% : {}/{} process done".format(actual_number_folder/float(number_total)*100, actual_number_folder, number_total), end="\r", flush=True)

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
	:return: The report table with the new name at the end
	:rtype: pandas.DataFrame
	"""

	nb_fasta = len(folderFasta)

	new_report_table = extract_protein(fileReport, INFO, PROTEIN_FUNCTION)
	nb_replicon_wanted = len(set(new_report_table.Replicon_name))

	print()
	print("Multi processing extraction ...")
	print()

	number_process_done = 0

	with Pool(processes=nb_thread) as pool:
		big_list = [(fasta, listOfFile, new_report_table, PROTEIN_FUNCTION, nb_replicon_wanted) for fasta in folderFasta if os.path.splitext(os.path.basename(fasta))[0] in list(new_report_table.Replicon_name)]
		#func=partial(write_fasta_multithreads, listOfFile=listOfFile, wanted=new_report_table, PROTEIN_FUNCTION=PROTEIN_FUNCTION), iterable=folderFasta
		result = pool.starmap(write_fasta_multithreads, iterable=big_list)
		pool.close()
		pool.join()

	print()
	print("All processes : Done !")
	print()

	for myfile in listOfFile :
		file_name = os.path.basename(myfile)
		list_tmp = glob.glob(os.path.join(os.path.dirname(listOfFile[0]), "tmp", "*", file_name))

		# XXX Je concatène toutes les séquences extraitent
		print("Writing {} ... ".format(os.path.basename(myfile)), end="", flush=True)

		open(myfile, "w").close()

		for subfile in list_tmp :
			subprocess.call("cat {} >> {}".format(subfile, myfile),shell=True)

		print("Done !")

	# XXX Je supprime tout le dossier tmp
	shutil.rmtree(os.path.join(os.path.dirname(listOfFile[0]), "tmp"))

	print("\n#################")
	print("# File wrote")
	print("#################\n")

	new_report_table.to_csv(os.path.join(INFO, "report_modif","detected.report"), sep="\t", index=False)
	return new_report_table

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
	:return: A report table like that contain only the gene of the system + the system_name (REPLICON_SYSTEM_NUMBER)
	:rtype: Pandas.DataFrame
	"""


	print("\n#################")
	print("# Verified Fasta")
	print("#################\n")

	report_like = pd.DataFrame(columns=["System_Id","Gene"])

	list_handle = [open(my_file, 'w') for my_file in listOfFile]
	dict_count = {}

	w_file = open(os.path.join(INFO, "translation_table_verified.tab"), 'w')
	w_file.write("#Name_fasta_file\tNew_name\n")

	info_extract = pd.read_table(info_dat, index_col=0, names=["Gene","System","SystID","Family","Note","Note2","NewName"], comment="#")

	# XXX Ce fichier est juste là car j'ai un soucis de plus d'un système par systèmes validés
	problem = False
	if os.path.exists("/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/verifed_duplicata.txt") :
		df_translate = pd.read_table("/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/verifed_duplicata.txt", index_col=0)
		problem = True
	else :
		print("There is only one system by replicon")

	progression=1

	seqiter = SeqIO.parse(data_fasta, "fasta")
	for seq in seqiter :
		if seq.id in info_extract.index :

			sys.stdout.write("{:.2f}% : {}/{} sequences wanted found\r".format(progression/float(info_extract.shape[0])*100, progression,info_extract.shape[0]))
			sys.stdout.flush()
			progression += 1

			# IDEA pour faire un truc général je pourrais par exemple crée la liste à partir du fichier de definition de protein_function.def

			tmp_replicon = info_extract.loc[seq.id, "NewName"].split("_")[0]
			num_tmp = 1
			if problem and info_extract.loc[seq.id, "SystID"] in list(df_translate.index) and df_translate.loc[info_extract.loc[seq.id, "SystID"], "Num"] == 2:
				tmp_replicon = "{}_{}".format(tmp_replicon, 2)
				num_tmp = 2

			if info_extract.loc[seq.id, "Gene"].split("_")[0] in ['T2SS','T4P', 'Tad', "Com"]:
				if info_extract.loc[seq.id, "Gene"] in PROTEIN_FUNCTION :

					gene_tmp = info_extract.loc[seq.id, "Gene"]

					writing_file = re.search("[a-zA-Z0-9/_]+{}\.fasta".format(PROTEIN_FUNCTION[info_extract.loc[seq.id, "Gene"]]), "\t".join(listOfFile)).group(0)


					NewName = "{}_{}_V_{}".format(tmp_replicon, info_extract.loc[seq.id, "SystID"].split("_")[-1], "_".join(info_extract.loc[seq.id, "Gene"].split("_")[1:]))

					if NewName in dict_count :
						dict_count[NewName] += 1
						NewName = "{}_Num{}_{}_V_{}".format(tmp_replicon, dict_count[NewName], info_extract.loc[seq.id, "SystID"].split("_")[-1], "_".join(info_extract.loc[seq.id, "Gene"].split("_")[1:]))
					else :
						dict_count[NewName] = 1

					# XXX je fais aussi un fichier translation pour les séquences vérifiées
					w_file.write("{}\t{}\n".format(seq.id, NewName))

					seq.id = NewName
					seq.name = ""
					seq.description = ''

					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
				else :
					continue

			else :
				# NOTE Permet d'avoir le bon nom si dans la colonne 2 j'ai que gspD par exemple et comme ça je reforme T2SS_gspD pour les certains systemes ...
				gene_tmp = info_extract.loc[seq.id, "System"]+"_"+info_extract.loc[seq.id, "Gene"]

				if gene_tmp in PROTEIN_FUNCTION :
					writing_file = re.search('[/a-zA-Z0-9_]+{}\.fasta'.format(PROTEIN_FUNCTION[gene_tmp]), "\t".join(listOfFile)).group(0)

					NewName = "{}_{}_V_{}".format(tmp_replicon, info_extract.loc[seq.id, "System"].split("_")[-1], info_extract.loc[seq.id, "Gene"])

					if NewName in dict_count :
						dict_count[NewName] += 1
						NewName = "{}_Num{}_{}_V_{}".format(tmp_replicon, dict_count[NewName], info_extract.loc[seq.id, "System"].split("_")[-1], info_extract.loc[seq.id, "Gene"])
					else :
						dict_count[NewName] = 1

					# XXX je fais aussi un fichier translation pour les séquences vérifiées
					w_file.write("{}\t{}\n".format(seq.id, NewName))

					seq.id = NewName
					seq.name = ''
					seq.description = ''
					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
				else :
					continue


			NewName_split = NewName.split("_")
			index_V = NewName_split.index("V")
			report_like.loc[-1, "NewName"] = NewName
			report_like.loc[-1, "System_Id"] = "{}_{}_{}".format(NewName_split[0], NewName_split[index_V-1], num_tmp)
			report_like.loc[-1, "Gene"] = gene_tmp
			report_like.index = report_like.index+1

	report_like.reset_index(drop=True, inplace=True)

	print()
	print("Done!")

	for open_file in list_handle :
		open_file.close()

	# XXX On ferme le fichier de translation_table
	w_file.close()
	report_like.to_csv(os.path.join(INFO, "report_modif", "verified.report"),sep="\t", index=False)
	return report_like

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
			id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("D")+1]))

			if id_seq in dict_remove :
				continue

			elif seq.seq in list_seq_verified :
				index=list_seq_verified.index(seq.seq)

				id_seq_verif = list_id_verified[index].split("_")
				id_seq_verif = re.sub("Num[0-9]_", "", "_".join(id_seq_verif[:id_seq_verif.index("V")+1]))

				# NOTE dans le dictionnaire je met le système vérifié en premier, toutes les séquences du système identitique en deuxième et la séquence qui en est la cause en troisième
				dict_remove[id_seq]=[id_seq_verif,[], seq.id]

			elif seq.seq in list_seq_detected_OK :
				index=list_seq_detected_OK.index(seq.seq)

				id_seq_reference = list_id_detected_OK[index]

				# NOTE dans le dictionnaire je met le système référence en premier, toutes les séquences du système identitique en deuxième et la séquence qui en est la cause en troisième
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
					id_seq=re.sub("Num[0-9]_", "", "_".join(id_seq[:id_seq.index("D")+1]))

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
