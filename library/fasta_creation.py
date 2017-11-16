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
import re
import os
import shutil
import pandas as pd
import glob
from multiprocessing import Pool
import subprocess
import shlex
import time
from set_params import *

from macsypy.config import Config
from macsypy.system import system_bank
from macsypy.gene import gene_bank
from macsypy.system_parser import SystemParser

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def get_right_name(gene_name, gene_system):
	"""
	Function that translate the name of all the genes homologues found in th analysis.
	Could by more generic by asking the gene for the user or with a xml with the name of all the gene in the description not in homologs only.

	:param gene_name: name of the gene to test
	:type: str
	:param gene_system: the instance of the system predicted of the gene
	:type: macsypy.system.System
	:return: the good name for the gene
	:rtype: str

	"""

	#if gene_name.split("_")[0] !=  gene_system.name:
	ref_gene = gene_system.get_gene_ref(gene_system.get_gene(gene_name))

	if ref_gene and ref_gene.exchangeable :
		return ref_gene.name

		'''
		genes = gene_system.mandatory_genes + gene_system.accessory_genes
		for gene in genes :
			if gene.exchangeable :
				homologs = gene.get_homologs()
				for homolog in homologs :
					if gene_name == homolog.name :
						return gene.name'''
	return gene_name


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

	if "generic" != row.Predicted_system :
		gene_tmp = "_".join(get_right_name(row.Gene, system_bank[row.Predicted_system]).split("_")[1:])
	else :
		gene_tmp = "_".join(row.Gene.split("_")[1:])

	#On met ici le numero du system après le replicon name
	NewNameTmp = "{}_{}".format(row.Replicon_name, row.System_Id.split('_')[-1])

	NewName = "{}_{}_D_{}".format(NewNameTmp, row.Predicted_system, gene_tmp)

	if NewName in dict_count :
		dict_count[NewName] += 1
		NewName = "{}_Num{}_{}_D_{}".format(NewNameTmp, dict_count[NewName], row.Predicted_system, gene_tmp)
	else :
		dict_count[NewName] = 1

	return NewName


##########################################################################################
##########################################################################################


def extract_protein(fileReport, INFO, PROTEIN_FUNCTION, config_file):

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
	:param config_file: the name of the config file of macsyfinder
	:type: str
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

	# NOTE New name : AAAAKKK.B.LLLLL_numero de systeme[_Num(nombre de fois nom trouvé)]_nomSysteme_D_nomProteine (for 2016 gembase format)
	# NOTE New name : NC_XXXXXX[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_D_nomProteine (for 2015 gembase format)
	# NOTE New name : NNNN[_numero de systeme si deux systemes trouvés][_Num(nombre de fois nom trouvé)]_nomSysteme_V_nomProteine (for 2013 gembase format)

	list_systems = new_report_table.Predicted_system.unique().tolist()
	list_systems.remove("generic")
	out_tmp = "tmp_{}_fasta".format(time.strftime("%Y%m%d"))
	SystemParser(Config(cfg_file=config_file, out_dir=out_tmp), system_bank, gene_bank).parse(list_systems)

	dict_count = {}
	new_report_table["NewName"] = new_report_table.apply(set_name, args=[dict_count], axis=1)
	new_report_table["NewName"] = new_report_table.apply(lambda x: x.NewName if ("_Num" in x.NewName or dict_count[x.NewName] == 1) else "{}_Num1_{}".format("_".join(x.NewName.split("_")[:2]), "_".join(x.NewName.split("_")[2:])), axis=1)

	# XXX Je crée une table de traduction entre mon nom et le nom générique de la séquence et je change les nouveaux nom pour le nom final plus le gene pour avoir l'information du changement de nom
	new_report_table.loc[:,["Hit_Id","NewName", "Gene"]].to_csv(os.path.join(INFO, "translation_table_detected.tab"),sep="\t", index=False)

	return new_report_table

##########################################################################################
##########################################################################################

def find_in_fasta(fileFasta, fileReport, listOfFile, INFO, PROTEIN_FUNCTION, config_file):

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
	:param config_file: the name of the config file of macsyfinder
	:type: str
	:return: The report table with the new name at the end
	:rtype: pandas.DataFrame
	"""

	list_handle = [open(my_file,"w") for my_file in listOfFile]

	report_table = extract_protein(fileReport, INFO, PROTEIN_FUNCTION, config_file)
	shutil.rmtree("tmp_{}".format(time.strftime("%Y%m%d")))
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
	:param number_total: total number of Replicon
	:type: int
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

def find_in_fasta_multithreads(folderFasta, fileReport, listOfFile, INFO, PROTEIN_FUNCTION, nb_thread, config_file):

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
	:param config_file: the name of the config file of macsyfinder
	:type: str
	:return: The report table with the new name at the end
	:rtype: pandas.DataFrame
	"""

	nb_fasta = len(folderFasta)

	new_report_table = extract_protein(fileReport, INFO, PROTEIN_FUNCTION, config_file)
	nb_replicon_wanted = len(set(new_report_table.Replicon_name))

	print()
	print("Multi processing extraction ...")
	print()

	number_process_done = 0

	with Pool(processes=nb_thread) as pool:
		# Je fais une liste avec le tous les arguments pour chaques fasta dans le dossier de genome
		big_list = [(fasta, listOfFile, new_report_table, PROTEIN_FUNCTION, nb_replicon_wanted) for fasta in folderFasta if os.path.splitext(os.path.basename(fasta))[0] in list(new_report_table.Replicon_name)]
		#func=partial(write_fasta_multithreads, listOfFile=listOfFile, wanted=new_report_table, PROTEIN_FUNCTION=PROTEIN_FUNCTION), iterable=folderFasta

		# Je lance mon multiprocesseur sur chaque association argument/fasta de ma grande liste
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

def set_validated_newname(row, dict_count):

	"""
	Function that set the name of validated sequences

	:param row: The row of the dataframe that contain all the information of the information about the validated sequences
	:type: pandas.Series
	:param dict_count: a dictionnary that will contain all the newname set to know if the newname is set yet.
	:type: dict
	:return: the newname set
	:rtype: str
	"""

	REPLICON_SYSTEM_NUMBER = "{}_{}".format(row.Replicon_name, row.System_Id.split("_")[-2])
	NewName = "{}_{}_V_{}".format(REPLICON_SYSTEM_NUMBER, row.System_Id.split("_")[1], "_".join(row.Gene.split("_")[1:]))

	if NewName in dict_count :
		dict_count[NewName] += 1
		NewName = "{}_Num{}_{}_V_{}".format(REPLICON_SYSTEM_NUMBER, dict_count[NewName], row.System_Id.split("_")[1], "_".join(row.Gene.split("_")[1:]))
	else :
		dict_count[NewName] = 1

	return NewName

##########################################################################################
##########################################################################################

def create_validated_fasta(listOfFile, PROTEIN_FUNCTION, data_fasta, info_dat, INFO):

	# NOTE Pas du tous adapté si on a un programme général si des verifiés autre que les miens c'est foutu

	"""
	Function used to extract the validated sequences for ATPase, prepilin peptidase, pilin (major and minor), IM platform

	:param listOfFile: list of all the file where the sequences will be write (one for each kind of protein)
	:type: list of str
	:param PROTEIN_FUNCTION: dictionnary return by the function set_params.set_dict_cutoff
	:type: dict
	:data_fasta: Fasta file with the validated sequence of the systems
	:type: str
	:info_dat:File with the information about the validated systems with this information : #SeqID Gene System SystID
	:type: str
	:param INFO: absolute path of the info_folder
	:type: str
	:return: A report table like that contain only the gene of the system + the system_name (REPLICON_SYSTEM_NUMBER)
	:rtype: Pandas.DataFrame
	"""


	print("\n#################")
	print("# Validated Fasta")
	print("#################\n")

	report_like = pd.DataFrame(columns=["NewName","Hit_Id","Replicon_name","Sequence_length","Gene", "Reference_system", "Predicted_system","System_Id"])

	list_handle = [open(my_file, 'w') for my_file in listOfFile]
	dict_count = {}

	w_file = open(os.path.join(INFO, "translation_table_validated.tab"), 'w')
	w_file.write("#Name_fasta_file\tNew_name\n")

	# NOTE Ici à l'importation des données, je ne donne pas de nom à la première colonne du .dat donc il prend la première colonne en tant qu'index
	info_extract = pd.read_table(info_dat, index_col=0, names=["Replicon_name","Gene","System_name","System_Id","Family","In_gembases", "Species_name", "Kingdom", "Phylum", "Notes"], comment="#")

	"""
	# XXX Ce fichier est juste là car j'ai un soucis de plus d'un système par systèmes validés
	problem = False
	if os.path.exists("/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/verifed_duplicata.txt") :
		df_translate = pd.read_table("/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/verifed_duplicata.txt", index_col=0)
		problem = True
	else :
		print("There is only one system by replicon")
	"""

	progression=1

	info_extract = info_extract[info_extract.Gene.isin(PROTEIN_FUNCTION)]
	df_newname = info_extract.apply(set_validated_newname, args=[dict_count], axis=1)
	df_newname = df_newname.apply(lambda x: x if ("_Num" in x or dict_count[x] == 1) else "{}_Num1_{}".format("_".join(x.split("_")[:2]), "_".join(x.split("_")[2:])))
	df_newname.reset_index().to_csv(w_file, sep="\t", index=False, header=False)

	#Pour facilité le loc après
	df_newname = df_newname.to_frame()
	df_newname.columns = ["NewName"]

	#df_newname.to_excel("pb_count_2.xlsx")

	seqiter = SeqIO.parse(data_fasta, "fasta")
	for seq in seqiter :
		if seq.id in df_newname.index :

			sys.stdout.write("{:.2f}% : {}/{} sequences wanted found\r".format(progression/float(df_newname.shape[0])*100, progression,df_newname.shape[0]))
			sys.stdout.flush()
			progression += 1

			# IDEA pour faire un truc général je pourrais par exemple crée la liste à partir du fichier de definition de protein_function.def
			writing_file = re.search("[a-zA-Z0-9/_]+{}\.fasta".format(PROTEIN_FUNCTION[info_extract.loc[seq.id, "Gene"]]), "\t".join(listOfFile)).group(0)

			old_name = seq.id

			seq.id = df_newname.loc[seq.id, 'NewName']
			seq.name = ""
			seq.description = ''

			SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")

			NewName_split = seq.id.split("_")
			index_V = NewName_split.index("V")
			report_like.loc[-1, "NewName"] = seq.id
			report_like.loc[-1, "Hit_Id"] = old_name
			report_like.loc[-1, "Sequence_length"] = len(seq.seq)
			report_like.loc[-1, "System_Id"] = info_extract.loc[old_name, "System_Id"]
			report_like.loc[-1, "Replicon_name"] = info_extract.loc[old_name, "Replicon_name"]
			report_like.loc[-1, 'Predicted_system'] = info_extract.loc[old_name, "System_name"]
			report_like.loc[-1, 'Reference_system'] = info_extract.loc[old_name, "System_name"]
			report_like.loc[-1, "Gene"] = info_extract.loc[old_name, "Gene"]
			report_like.index = report_like.index+1


	report_like.reset_index(drop=True, inplace=True)

	print()
	print("Done!")

	for open_file in list_handle :
		open_file.close()

	# XXX On ferme le fichier de translation_table
	w_file.close()
	report_like.to_csv(os.path.join(INFO, "report_modif", "validated_tmp.report"),sep="\t", index=False)

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
	know which sequences are extract and which validated sequences are the same
	of this one.

	:param dict_remove: dictionnary create by the function concatenate_detected_validated()
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


def concatenate_detected_validated(fasta_name, PATH_FASTA_DETECTED, PATH_FASTA_validated, INFO_folder, PATH_FASTA_CONCATENATED, PATH_FASTA_DETECTED_SELECTED):

	"""
	Function that concatenate the validated and detected file and remove detected sequences
	that are already in the validated file. It write a file in the information folder in
	markdown format to know which sequences are extract and which validated sequences are
	the same of this one.

	:param fasta_name: the name of all the fasta file create ([protein_function].fasta)
	:type: list of str
	:param PATH_FASTA_DETECTED: absolute path to detected fasta folder
	:type: str
	:param PATH_FASTA_validated: absolute path to validated fasta folder
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

	# NOTE Dictionnaire avec en clef l'id espèce/système et en value une liste : pour le critère 100% similaire
	# NOTE ["l'id espèce/système du verifié qui correspond", [liste des sequences ATPase, IM ...]]
	dict_remove = {}

	# NOTE Dictionnaire avec en clef l'id espèce/système et en value une liste : pour le critère au moins 1 ATPase et 1 IMplatform
	# NOTE ["l'id espèce/système du verifié qui correspond", [liste des sequences ATPase, IM ...]]
	dict_remove_ATP_IM = {}

	print("\n------------------------------------------")
	print("| First read : Creation of the dictionnary")
	print("------------------------------------------\n")

	for fasta_file in fasta_name :
		validated_fasta=os.path.join(PATH_FASTA_validated, fasta_file)
		detected_fasta=os.path.join(PATH_FASTA_DETECTED, fasta_file)
		concatenated_fasta=os.path.join(PATH_FASTA_CONCATENATED, fasta_file)

		list_seq_validated = list(SeqIO.parse(validated_fasta, "fasta"))
		list_id_validated = [seq.id for seq in list_seq_validated]
		list_seq_validated = [seq.seq for seq in list_seq_validated]

		# NOTE liste des sequence et identifiant qui sont à garder
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
			id_seq=re.sub("Num[0-9]+_", "", "_".join(id_seq[:id_seq.index("D")+1]))

			if id_seq in dict_remove :
				continue

			elif seq.seq in list_seq_validated :
				index=list_seq_validated.index(seq.seq)

				id_seq_verif = list_id_validated[index].split("_")
				id_seq_verif = re.sub("Num[0-9]+_", "", "_".join(id_seq_verif[:id_seq_verif.index("V")+1]))

				# NOTE dans le dictionnaire je met le système vérifié en premier, toutes les séquences du système identitique en deuxième et la séquence qui en est la cause en troisième
				dict_remove[id_seq]=[id_seq_verif,[], seq.id]

			elif seq.seq in list_seq_detected_OK :
				index=list_seq_detected_OK.index(seq.seq)

				id_seq_reference = list_id_detected_OK[index]

				# NOTE dans le dictionnaire je met le système référence en premier, toutes les séquences du système identitique en deuxième et la séquence qui en est la cause en troisième
				dict_remove[id_seq] = [id_seq_reference,[], seq.id]

			else :
				list_seq_detected_OK.append(seq.seq)
				list_id_detected_OK.append(id_seq)

		print()
		print("File : {} -> Done!".format(fasta_file))

	print("\n-----------------------------")
	print("| Second read : Writing files")
	print("-----------------------------\n")

	for fasta_file in fasta_name :
		validated_fasta = os.path.join(PATH_FASTA_validated, fasta_file)
		detected_fasta = os.path.join(PATH_FASTA_DETECTED, fasta_file)
		new_detected_fasta = os.path.join(PATH_FASTA_DETECTED_SELECTED, fasta_file)
		concatenated_fasta = os.path.join(PATH_FASTA_CONCATENATED, fasta_file)

		os.system('cat "{}" > "{}"'.format(validated_fasta, concatenated_fasta))

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
					id_seq=re.sub("Num[0-9]+_", "", "_".join(id_seq[:id_seq.index("D")+1]))

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

def concatenate_reduce(list_file_concatenated, PATH_FASTA_CONCATENATED_REMOVE_POOR, report_df_full, DICT_IMPORTANT_PROTEINS, info_folder, liste_detected_file = False, path_detected = False):

	"""
	Function that concatenate the validated and detected file and remove detected sequences
	that are already in the validated file. It write a file in the information folder in
	markdown format to know which sequences are extract and which validated sequences are
	the same of this one.

	:param list_file_concatenated: list of the concatenated fasta files
	:type: list of str
	:param PATH_FASTA_CONCATENATED_REMOVE_POOR: absolute path to concatenated fasta folder whete to write the new files
	:type: str
	:param report_df_full: the dataframe that contain the information about each proteins of the analysis
	:type: pandas.DataFrame
	:param DICT_IMPORTANT_PROTEINS: dictionnary of the functions in key and proteins wanted in value
	:type: dict
	:param info_folder: The path to the information folder
	:type: str
	:param liste_detected_file: list of the detected fasta files
	:type: str
	:param path_detected: absolute path to detected fasta folder whete to write the new files
	:type: str
	:return: Nothing
	"""

	print("\n#################################################")
	print("# Remove systems without ATPase and/or IMplatform")
	print("#################################################\n")

	new_report = report_df_full.copy()

	with open(os.path.join(info_folder, "remove_not_complete_systems.report"), "w") as w_file :
		for function in DICT_IMPORTANT_PROTEINS :
			w_file.write("# These systems are remove because they have not at least protein of with the function {}\n".format(function))
			list_wanted_system = new_report.System_Id[new_report.Gene.isin(DICT_IMPORTANT_PROTEINS[function])].tolist()
			new_report[~(new_report.System_Id.isin(list_wanted_system))].to_csv(w_file, sep="\t", index=False, header=False)
			new_report = new_report[new_report.System_Id.isin(list_wanted_system)]


	if liste_detected_file :
		type_file = "concatenated"
	else :
		type_file = "detected"

	print("\n################################")
	print("# Write the new {} files".format(type_file))
	print("################################\n")


	for myfile in list_file_concatenated :

		seq_parser = SeqIO.parse(myfile, "fasta")
		number_seq = len(list(seq_parser))
		progression = 1

		seq_parser = SeqIO.parse(myfile, format="fasta")
		with open(os.path.join(PATH_FASTA_CONCATENATED_REMOVE_POOR, os.path.basename(myfile)), 'w') as w_file :
			for seq in seq_parser :

				print("File : {} -> {:.2f}% : {}/{} sequences read".format(os.path.basename(myfile), progression/float(number_seq)*100, progression,number_seq), end='\r', flush=True)
				progression += 1

				if seq.id in new_report.NewName.tolist() :
					SeqIO.write(seq, w_file, format="fasta")

		print()
		print("File : {} -> Done!".format(os.path.basename(myfile)))

	if liste_detected_file :

		print("\n#############################")
		print("# Write the new detected files")
		print("#############################\n")


		for myfile in liste_detected_file :

			seq_parser = SeqIO.parse(myfile, "fasta")
			number_seq = len(list(seq_parser))
			progression = 1

			seq_parser = SeqIO.parse(myfile, format="fasta")
			with open(os.path.join(path_detected, os.path.basename(myfile)), 'w') as w_file :
				for seq in seq_parser :

					print("File : {} -> {:.2f}% : {}/{} sequences read".format(os.path.basename(myfile), progression/float(number_seq)*100, progression,number_seq), end='\r', flush=True)
					progression += 1

					if seq.id in new_report.NewName.tolist() :
						SeqIO.write(seq, w_file, format="fasta")

			print()
			print("File : {} -> Done!".format(os.path.basename(myfile)))

	print()
	print("Done!")
	print()

	return
