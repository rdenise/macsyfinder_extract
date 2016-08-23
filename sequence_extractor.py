# -*- coding: utf-8 -*-
import argparse
from textwrap import dedent
import sys, os

sys.path.insert(0, os.path.abspath('library'))

from set_params import *
from fasta_creation import *
from merge_generique_all import *

##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""
     *            *               *                   *
*           *               *   *   *  *    **                *   *
  **     *    *   *  *     *                    *               *
    __  __  *              ____ *        *  *  *    **     *
|| |  \/  | __ _  ___  || / ___| _   _  ||   ___ _         _        *
|| | |\/| |/ _` |/ __| || \___ \| | | | ||  | __(_)_ _  __| |___ _ _
|| | |  | | (_| | (__  ||  ___) | |_| | ||  | _|| | ' \/ _` / -_) '_|
|| |_|  |_|\__,_|\___| || |____/ \__, | ||  |_| |_|_||_\__,_\___|_|
           *             *       |___/         *                   *
 *      *   * *     *   **         *   *  *           *
  *      *         *        *    *              *
             *                           *  *           *     *
||  ____                                         _____      _                  _
|| / ___|  ___  __ _ _   _  ___ _ __   ___ ___  | ____|_  _| |_ _ __ __ _  ___| |_ ___  _ __
|| \___ \ / _ \/ _` | | | |/ _ \ '_ \ / __/ _ \ |  _| \ \/ / __| '__/ _` |/ __| __/ _ \| '__|
||  ___) |  __/ (_| | |_| |  __/ | | | (_|  __/ | |___ >  <| |_| | | (_| | (__| || (_) | |
|| |____/ \___|\__, |\__,_|\___|_| |_|\___\___| |_____/_/\_\\\__|_|  \__,_|\___|\__\___/|_|
||                |_|


""") )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-s",'--seqdata',
 							required=True,
							metavar="<file>",
							dest="seqData",
							help="File with all the fasta sequences used for macsyfinder analysis")
general_option.add_argument("-r",'--reportfile',
 							required=True,
							metavar="<file>",
							dest="reportFile",
							help="Report file from the macsyfinder analysis")
general_option.add_argument("-d",'--defFile',
 							required=True,
							metavar="<file>",
							dest="defFile",
							help="Tabular file with the the function name in first column and name of the protein in the other")
general_option.add_argument("-pre",'--prefix',
 							default=None,
							dest="prefix",
							metavar='<PREFIX>',
							help="Using <PREFIX> for output files (default: reportFile directory)")

extraction_option = parser.add_argument_group(title = "Extraction options")
extraction_option.add_argument("-cut",'--cutoff',
							metavar="<CUTOFF_FILE>",
							nargs='?',
							dest="cutoff",
							default=None,
							help="Option to remove sequences that are way to much longer or shorter beside of the other. (If a file isn't gave it will calculate the cutoff and generate a file)")
extraction_option.add_argument("-veriFile",'--verifiedFile',
							metavar="<VERIFIED_FASTA_FILE>",
							dest="veriFile",
							default=None,
							help="File with the sequence in fasta of the verified file")
extraction_option.add_argument("-veriData",'--verifiedData',
							metavar="<VERIFIED_DATA>",
							dest="veriData",
							default=None,
							help="File with the information about the verified systems with this information : #SeqID	Gene	System	SystID	Family	Note")
extraction_option.add_argument("-conc",'--concatenate',
							dest="concat",
							default=None,
							help="Allow to concatenate detected sequences and verified sequences")

merge_option = parser.add_argument_group(title = "Merge report options")
merge_option.add_argument("-m",'--merge',
							metavar=("<OTHER_REPORT>", "<GENERIQUE_REPORT>"),
							nargs=2,
							dest="merge",
							default=None,
							help="Merge the generique .report and the other systems together without add a systems generique that is in the OTHER_REPORT. Write the new .report in GENERIQUE_REPORT directory")

args = parser.parse_args()

if args.merge :
	fileAll = os.path.abspath(args.merge[0])
	fileGenerique = os.path.abspath(args.merge[1])

	fileWrite = os.path.join(os.path.dirname(fileGenerique),"merge_macsyfinder.report")

	write_merge_file(fileGenerique, fileAll, fileWrite)

	print "\n#################"
	print "# File merged"
	print "#################\n"

	sys.exit(0)



if not args.prefix :
	PREFIX = os.path.join(os.path.abspath(args.reportFile),"extraction_%s" %(time.strftime("%d_%m_%y")))
else :
	PREFIX = args.prefix

create_folder(PREFIX)

#Creation of an information folder for each sequence remove or file generate (as cutoff, ...)
INFO=os.path.join(PREFIX,"info_folder")
create_folder(INFO)

FASTA = args.seqData
REPORT = args.reportFile

PROTEIN_FUNCTION = read_protein_function(args.defFile)

#List of all the function name with .fasta add at the end
all_function = [function+".fasta" for function in PROTEIN_FUNCTION.keys()]
list_file = robjects.StrVector(all_function)

#Si j'ai l'option verifié mise en place je demande les deux options
if args.veriFile or args.veriData :
	if not (args.veriFile and args.veriData):
		parser.error("you MUST provided a verified fasta file and a annotation data file. If you want verified fasta")
	else :
		PATH_FASTA_VERIFIED = os.path.join(PREFIX, "fasta_verified")
		create_folder(PATH_FASTA_VERIFIED)
		create_verified_fasta(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep=''), PROTEIN_FUNCTION)
		rename_name_gene(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep=''))


# Première liste de fichiers détectés
PATH_FASTA_DETECTED = os.path.join(PREFIX, "fasta_detected", "raw")
create_folder(PATH_FASTA_DETECTED)
list_file_detected = robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep='')

find_in_fasta(FASTA, REPORT, list_file_detected, INFO, PROTEIN_FUNCTION)

# Deuxième liste de fichiers détectés après que tous les nom soit renomé
PATH_FASTA_RENAME = os.path.join(PREFIX, "fasta_detected", "rename")
rename_name_gene(list_file_detected, PATH_FASTA_RENAME)
PATH_FASTA_DETECTED = PATH_FASTA_RENAME
list_file_detected = robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep='')

if args.concat :
	if args.cutoff :
		PATH_FASTA_CONCATENATED = os.path.join(PREFIX, "fasta_concatenated", "raw")
	else :
		PATH_FASTA_CONCATENATED = os.path.join(PREFIX, "fasta_concatenated")

	create_folder(PATH_FASTA_CONCATENATED)

	concatenate_detected_verified(list_file, PATH_FASTA_DETECTED, PATH_FASTA_VERIFIED, INFO, PATH_FASTA_CONCATENATED)

# Deuxieme liste de fichiers concaténés ou detectés après cutoff
if args.cutoff :
	if args.concat :
		PATH_FASTA_CONCATENATED_CUTOFF = os.path.join(PREFIX, "fasta_concatenated", "cut_off")
		cut_seq_fasta_file(list_file_detected, PATH_FASTA_CONCATENATED_CUTOFF, INFO, file_cutoff=args.cutoff)
	else :
		PATH_FASTA_DETECTED_CUTOFF = os.path.join(PREFIX, "fasta_detected", "cut_off")
		cut_seq_fasta_file(list_file_detected, PATH_FASTA_DETECTED_CUTOFF, INFO, file_cutoff=args.cutoff)

print "\n#################"
print "# End"
print "#################\n"
