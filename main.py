# -*- coding: utf-8 -*-

from fasta_creation import *
import argparse

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

#cmd_cat = 'for file in `ls %s` ; do cat "%s$file" "%s$file" > "%s$file" ; done' % (PATH_FASTA_VERIFIED, PATH_FASTA_DETECTED,PATH_FASTA_VERIFIED,PATH_FASTA_CONCATENATED)
#os.system(cmd_cat)

#cut_seq_fasta_file(robjects.r['paste'](PATH_FASTA_CONCATENATED, list_file, sep=''))

#rename_name_gene(robjects.r['paste'](PATH_FASTA_CONCATENATED, list_file, sep=''))

## VIEUX TRUC INUTILE
##list_rename = find_rename_fasta(robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep=''))


## SERT AUSSI A RIEN CAR PAS DE DOUBLON ET DANS VERIFIED DEJA FAIT
##rename_name_gene(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep=''))
##rename_name_gene(list_file_detected)
