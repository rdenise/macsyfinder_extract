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
import time

##########################################################################################
##########################################################################################
##
##								Global variable
##
##########################################################################################
##########################################################################################

PATH_FASTA_DETECTED = "/Users/rdenise/Documents/Analysis_tree/fasta_detected/%s" % time.strftime("%d_%m_%y")

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


def extract_protein(fileReport, Protein):
	
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
		if report_table[j][4] == Protein :
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
	print "There are %i proteins remove during this operation because they are not the protein you want %s" %(number_remove_protein, Protein)
		
	new_name = [report_table[i][1]+'_'+report_table[i][6]+'_D_'+"_".join(report_table[i][4].split('_')[1:]) for i in xrange(report_table.shape[0])]
	
	return report_table[:,0].tolist(), new_name, report_table[:,4].tolist()


##########################################################################################


def find_in_fasta(fileFasta, fileReport, Protein) :
	
	"""
	This function is used to create the fasta with MacSyFinder hits found
	
	:param fileFasta: name of the fasta database used in the MacSyfinder analysis  
	:type: string  
	:param fileReport: name of the file .report of the MacSyFinder analysis  
	:return: Nothing
	"""
	
	final_file = open(os.path.join(PATH_FASTA_DETECTED,Protein+".fasta"), "w")
	
	wanted, name_genes, keys_genes = extract_protein(fileReport, Protein)
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
			
			if keys_genes[index] == Protein :
				SeqIO.write(seq , final_file, "fasta")
			else :
				sys.exit("ERROR:: Function Not Know : "+keys_genes[index])
	
	#Close all file
	final_file.close()
	
	bar.finish()
	
	print "\n#################"
	print "# File wrote"
	print "#################\n"

##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################

create_folder(PATH_FASTA_DETECTED)

fasta_file = "/Users/rdenise/Documents/Analysis_macsyfinder/database/Prokaryote_0615/alltogether.prot"
report_file = "/Users/rdenise/Documents/Analysis_macsyfinder/18_05_16/analysis_all_newdatabase_macsyfinder_18_05_16/macsyfinder.report"
protein = "T4SS_virb4"
find_in_fasta(fasta_file, report_file, protein)