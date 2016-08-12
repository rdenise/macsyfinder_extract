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
import progressbar
import rpy2.robjects as robjects
import re
import os

##########################################################################################
##########################################################################################
##
##								Global variable
##
##########################################################################################
##########################################################################################

PROTEIN_FUNCTION = {"T2SS_gspE":"ATPase",
					"T4P_pilB":"ATPase",
					"T4P_pilT_pilU":"ATPase",
					"Tad_tadA":"ATPase",
					"T2SS_gspF":"IMplatform",
					"T4P_pilC":"IMplatform",
					"Tad_tadB":"IMplatform",
					"Tad_tadC":"IMplatform",
					"T2SS_gspH":"minorPilin",
					"T2SS_gspI":"minorPilin",
					"T2SS_gspJ":"minorPilin",
					"T2SS_gspK":"minorPilin",
					"T4P_pilI_pilV":"minorPilin",
					"Tad_tadE":"minorPilin",
					"Tad_tadF":"minorPilin",
					"T2SS_gspG":"majorPilin",
					"T4P_pilAE":"majorPilin",
					"Tad_flp":"majorPilin",
					"T2SS_gspO":"prepilinPeptidase",
					"T4P_pilD":"prepilinPeptidase",
					"Tad_tadV":"prepilinPeptidase",
					"T2SS_gspD":"secretin",
					"T4P_pilQ":"secretin",
					"Tad_rcpA":"secretin",
					"Archaellum_FlaI":"ATPase",
					"Archaellum_FlaJ":"IMplatform",
					"Archaellum_FlaB":"majorPilin",
					"Archaellum_FlaK":"prepilinPeptidase"}


PATH_FASTA_VERIFIED = "/Users/rdenise/Documents/Analysis_tree/fasta_verify/rename/"

PATH_FASTA_DETECTED = "/Users/rdenise/Documents/Analysis_tree/fasta_detected/21_04_16/"

PATH_FASTA_DETECTED_CUTOFF = PATH_FASTA_DETECTED+"cut_off/"

PATH_FASTA_DETECTED_CLUSTERED = PATH_FASTA_DETECTED_CUTOFF+"clustered/"

PATH_FASTA_CONCATENATED = "/Users/rdenise/Documents/Analysis_tree/fasta_concatenated/21_04_16/"

PATH_FASTA_RENAME = "/Users/rdenise/Documents/Analysis_tree/fasta_detected/21_04_16/rename_fasta/"

PATH_SCRIPT_PYTHON = "/Users/rdenise/Documents/script_python"

DICT_CUTOFF = {"ATPase.fasta" : 690,
				"prepilinPeptidase.fasta" : 400,
				"IMplatform.fasta" : 450,
				"majorPilin.fasta" : 0,
				"minorPilin.fasta" : 500,
				"secretin.fasta" : 900,
				"Archaellum.fasta" : 0}

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



def extract_protein(fileReport):
	
	"""
	This function is used to select the sequence identified by MacSyFinder and create a fasta with these sequences. 
	
	:param fileReport: the file .report of the MacSyFinder analysis  
	:type: string  
	:return: the list of the sequence ids of the hit find by MacSyFinder, the list of all the new name for each sequences and the reference systems for each sequence.  
	:rtype: list of string, list of string, list of string  
	"""
	
	print "\n#################"
	print "# Protein extraction"
	print "#################\n"
	
	report_table = np.loadtxt(fileReport, dtype='string')
	number_prot = report_table.shape[0]
	index_remove = []
	
	for j in xrange(number_prot):
		if report_table[j][4] in PROTEIN_FUNCTION :
			system_number = report_table[j][7].split('_')[-1]
			if int(system_number) > 1 :
				report_table[j][1] = "_".join(report_table[j][1].split('_')[:-1])+"_"+system_number
			else:
				report_table[j][1] = "_".join(report_table[j][1].split('_')[:-1])
		else :
			index_remove.append(j)
	
	np.savetxt(os.path.join(PATH_FASTA_DETECTED, "remove_seq.seq"), report_table[np.array(index_remove),:], delimiter="\t", fmt="%s")
	report_table = np.delete(report_table, index_remove, axis=0)
	number_remove_protein = len(index_remove)
	print "There are %i proteins remove during this operation because they are not in the dictionnary" %number_remove_protein
		
	new_name = [report_table[i][1]+'_'+report_table[i][6]+'_D_'+"_".join(report_table[i][4].split('_')[1:]) for i in xrange(report_table.shape[0])]
	
	"""
	print "\n#################"
	print "# Selection protein"
	print "#################\n"	
	
	remove_from_new_name = []
	
	current_index = 0
	
	
	with open(PATH_FASTA_DETECTED+"view_delete_protein.log", "w") as log_file :
		for current_name in new_name:
			
			if current_index in remove_from_new_name :
				pass
			
			else :
				all_same = np.where(np.array(new_name) == current_name)[0].tolist()
				
				try :
					choice = np.argmin(report_table[all_same,10].astype(np.float))
				except IndexError :
					#print report_table[all_same,10].astype(np.float), all_same, current_name
					print all_same, current_name, len(new_name), len(report_table)
					sys.exit(0)
			
				seq_choice = all_same[choice]
				
				del all_same[choice]
			
				log_file.write("###\nChoose :\n")
				log_file.write("\t".join(report_table[seq_choice,:].tolist())+"\n")
				log_file.write("---\nDelete :\n")
				if all_same == [] :
					log_file.write("None\n")
				else :
					for seq_delete in all_same :
						log_file.write("\t".join(report_table[seq_delete,:].tolist())+"\n")
						remove_from_new_name.append(seq_delete)
			
			current_index += 1
	
	new_name = np.delete(np.array(new_name), np.array(remove_from_new_name)).tolist()		
	report_table = np.delete(report_table, np.array(remove_from_new_name), axis=0)
	"""		
	return report_table[:,0].tolist(), new_name, report_table[:,4].tolist()

##########################################################################################
##########################################################################################


def find_in_fasta(fileFasta, fileReport) :
	
	"""
	This function is used to create the fasta with MacSyFinder hits found
	
	:param fileFasta: name of the fasta database used in the MacSyfinder analysis  
	:type: string  
	:param fileReport: name of the file .report of the MacSyFinder analysis  
	:return: Nothing
	"""
	
	ATPase_file = open(PATH_FASTA_DETECTED+"ATPase.fasta", "w")
	Prep_pep_file = open(PATH_FASTA_DETECTED+'prepilinPeptidase.fasta',"w")
	IM_file = open(PATH_FASTA_DETECTED+'IMplatform.fasta','w')
	pilin_maj_file = open(PATH_FASTA_DETECTED+'majorPilin.fasta','w')
	pilin_min_file = open(PATH_FASTA_DETECTED+'minorPilin.fasta','w')
	secretin_file = open(PATH_FASTA_DETECTED+'secretin.fasta','w')

	
	wanted, name_genes, keys_genes = extract_protein(fileReport)
	seqiter = SeqIO.parse(open(fileFasta), 'fasta')                                    
	
	print "\n#################"
	print "# Writing ..."
	print "#################\n"
	
	bar = progressbar.ProgressBar(maxval=len(wanted), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	progression=1
	bar.start()
	
	for seq in seqiter :
	  if seq.id in wanted:
			bar.update(progression)                                                       
			progression+=1 
	    
			index = wanted.index(seq.id)
			seq.description = ''
			seq.name = name_genes[index]
			seq.id = seq.name
			
			if keys_genes[index] in PROTEIN_FUNCTION :
				if PROTEIN_FUNCTION[keys_genes[index]] == "ATPase" :
					SeqIO.write(seq , ATPase_file, "fasta")
				elif PROTEIN_FUNCTION[keys_genes[index]] == "prepilinPeptidase" :
					SeqIO.write(seq, Prep_pep_file, "fasta")
				elif PROTEIN_FUNCTION[keys_genes[index]] == 'IMplatform' :
					SeqIO.write(seq, IM_file, "fasta")
				elif PROTEIN_FUNCTION[keys_genes[index]] == 'majorPilin' :
					SeqIO.write(seq, pilin_maj_file, "fasta")
				elif PROTEIN_FUNCTION[keys_genes[index]] == 'minorPilin' :
					SeqIO.write(seq, pilin_min_file, "fasta")
				elif PROTEIN_FUNCTION[keys_genes[index]] == 'secretin' :
					SeqIO.write(seq, secretin_file, "fasta")
			else :
				sys.exit("ERROR:: Function Not Know : "+keys_genes[index])
	
	#Close all file
	ATPase_file.close()
	Prep_pep_file.close()
	IM_file.close()
	pilin_maj_file.close()
	pilin_min_file.close()
	secretin_file.close()
	
	bar.finish()
	
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

def rename_name_gene(listOfFile) :
	
	"""
	Function use to rename the sequence ids of sequence with the same name in the same file
	
	:param new_listOfFile: list of all the file where we need to check if there are a same sequence 
	id twice (or more) in the same file with absolute paths 
	:type: list  
	:return: Nothing
	"""
	
	print "\n#################"
	print "# Rename protein"
	print "#################\n"	
	
	new_listOfFile=[]
	
	for file in listOfFile :
		if os.stat(file).st_size != 0 :
			new_listOfFile.append(file)
	
	seq_to_rename = find_rename_fasta(new_listOfFile)
	dict_count = dict([(sequence[1:].rstrip(" "), 0) for sequence in seq_to_rename])
	bar = progressbar.ProgressBar(maxval=len(new_listOfFile), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	progression=1
	bar.start()
	
	create_folder(PATH_FASTA_RENAME)
	
	for file in new_listOfFile :
		
		file_name = os.path.basename(file)
		
		bar.update(progression)
		progression += 1

		handle = open(PATH_FASTA_RENAME+file_name, 'w')
		fasta_reading = SeqIO.parse(file, "fasta")
		for seq in fasta_reading :
			if seq.id in dict_count :
				if dict_count[seq.id] == 0 :
					dict_count[seq.id] += 1
				else :
					dict_count[seq.id] += 1
					if "NC_" in seq.id :
						seq.id = "_".join(seq.id.split("_")[:2])+"_"+str(dict_count[seq.id])+"_"+"_".join(seq.id.split("_")[2:])
					else :
						seq.id = seq.id.split("_")[0]+"_"+str(dict_count[seq.id])+"_"+"_".join(seq.id.split("_")[1:])
					seq.name = seq.id
					seq.description = ""
			SeqIO.write(seq, handle, "fasta")
		handle.close()
	bar.finish()
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
	
	print "\n#################"
	print "# Verified Fasta"
	print "#################\n"
	
	list_handle = []
	
	for file in listOfFile :
		list_handle.append(open(file, 'w'))
	
	info_extract = np.loadtxt("/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/all_seq_to_extract_annot.dat", dtype="string")
	
	bar = progressbar.ProgressBar(maxval=info_extract.shape[0], widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	progression=1
	bar.start()
	
	seqiter = SeqIO.parse("/Users/rdenise/Documents/de_sophie_a_remi/pour_remi/experiment_validated_systems/all_genes_all_systems.fasta", "fasta")
	
	for seq in seqiter :
		if seq.id in info_extract[:,0] :
			
			bar.update(progression)
			progression += 1
			
			position = info_extract[:,0].tolist().index(seq.id)
			
			if info_extract[position][1].split("_")[0] in ['T2SS','T4P', 'Tad']:
					
				if info_extract[position][1] in PROTEIN_FUNCTION :
					writing_file = re.search('[a-zA-Z/_]+'+PROTEIN_FUNCTION[info_extract[position][1]]+'\.fasta', "\t".join(listOfFile)).group(0)
					
					seq.name = info_extract[position][3]+"_V_"+"_".join(info_extract[position][1].split("_")[1:])
					seq.id = seq.name
					seq.description = ''
					
					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
						
			else :
				new_name = info_extract[position][2]+"_"+info_extract[position][1]
				
				if new_name in PROTEIN_FUNCTION :
					writing_file = re.search('[/a-zA-Z_]*'+PROTEIN_FUNCTION[new_name]+'\.fasta', "\t".join(listOfFile)).group(0)
					
					seq.name = info_extract[position][3]+"_V_"+info_extract[position][1]
					seq.id = seq.name
					seq.description = ''
					
					SeqIO.write(seq, list_handle[listOfFile.index(writing_file)], "fasta")
					
	bar.finish()
	
##########################################################################################
##########################################################################################

def cut_seq_fasta_file(listOfFasta) :
	
	"""
	Function used to remove some sequence in the concatenated fasta, that after futher
	analysis, are considered as not good.
	
	:param listOfFasta: List of all the file fasta where we need to remove sequences
	:type: list of string
	:return: Nothing
	"""
	
	print "\n#################"
	print "# Cut concatetaned file"
	print "#################\n"
		
	new_path = os.path.dirname(listOfFasta[0])+"/cut_off/"
	create_folder(new_path)
	
	for file in listOfFasta :
		current_file = file.split("/")[-1]
		
		if current_file in DICT_CUTOFF:
			
			if DICT_CUTOFF[current_file] != 0 :
				
				with open(new_path+current_file, "w") as writing_file :
					seqiter = SeqIO.parse(open(file), 'fasta') 
				
					for seq in seqiter :
						if len(seq) < DICT_CUTOFF[current_file] :
							SeqIO.write(seq, writing_file,"fasta")
			else :
				os.system(" ".join(["cp",file,new_path+current_file]))			
				

##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################

#create_folder(PATH_FASTA_DETECTED)

FASTA = sys.argv[1]
REPORT = sys.argv[2]

#find_in_fasta(FASTA, REPORT)

list_file = robjects.StrVector(("ATPase.fasta",'prepilinPeptidase.fasta','IMplatform.fasta','majorPilin.fasta','minorPilin.fasta', 'secretin.fasta'))

#create_verified_fasta(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep=''))

#list_file_detected = robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep='')

rename_name_gene(robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep=''))

#cut_seq_fasta_file(list_file_detected)

#execfile(PATH_SCRIPT_PYTHON+"selection_uclust.py")

#print "\n#################"
#print "# Concatetaned file"
#print "#################\n"

#create_folder(PATH_FASTA_CONCATENATED)

#cmd_cat = 'for file in `ls %s` ; do cat %s$file %s$file > %s$file ; done' % (PATH_FASTA_VERIFIED, PATH_FASTA_DETECTED,PATH_FASTA_VERIFIED,PATH_FASTA_CONCATENATED)
#os.system(cmd_cat)

#cut_seq_fasta_file(robjects.r['paste'](PATH_FASTA_CONCATENATED, list_file, sep=''))

#rename_name_gene(robjects.r['paste'](PATH_FASTA_CONCATENATED, list_file, sep=''))

## VIEUX TRUC INUTILE 
##list_rename = find_rename_fasta(robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep=''))


## SERT AUSSI A RIEN CAR PAS DE DOUBLON ET DANS VERIFIED DEJA FAIT
##rename_name_gene(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep=''))
##rename_name_gene(list_file_detected)
