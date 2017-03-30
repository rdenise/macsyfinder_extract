# -*- coding: utf-8 -*-
import argparse
from textwrap import dedent
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'library'))

from set_params import *
from fasta_creation import *
from merge_generic_all import *
from statistique_count import *

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
							metavar="<file(s)>",
							dest="seqData",
							help="File with all the fasta sequences used for macsyfinder analysis. If you want to multiprocessing put the name of the folder split by species")
general_option.add_argument("-r",'--reportfile',
							metavar="<file>",
							dest="reportFile",
							help="Report file from the macsyfinder analysis")
general_option.add_argument("-d",'--defFile',
 							required=True,
							metavar="<file>",
							dest="defFile",
							help="Tabular file with the the function name in first column and name of the protein in the other")
general_option.add_argument("-o",'--output',
 							default=None,
							dest="output",
							metavar='<OUTPUT>',
							help="Using <OUTPUT> for output files (default: reportFile directory)")

extraction_option = parser.add_argument_group(title = "Extraction options")
extraction_option.add_argument("-cut",'--cutoff',
							metavar="<CUTOFF_FILE>",
							nargs='?',
							const=True,
							dest="cutoff",
							default=None,
							help="Option to remove sequences that are way to much longer or shorter beside of the other. (If there is no file, it will calculate the cutoff and generate a file)")
extraction_option.add_argument("-veriFile",'--verifiedFile',
							metavar="<VERIFIED_FASTA_FILE>, <VERIFIED_DATA>",
							nargs=2,
							dest="veriFile",
							default=None,
							help="File with the sequence in fasta of the verified file and the file with the information about the verified systems with this information : #SeqID	Gene	System	SystID	Family	Note")
extraction_option.add_argument("-conc",'--concatenate',
							action='store_true',
							dest="concat",
							help="Allow to concatenate detected sequences and verified sequences")
extraction_option.add_argument("-w", "--worker_proc",
                            metavar="<NUMBER_OF_THREADS>",
                            dest="number_proc",
                            default=1,
                            help="Number of processor you want to used in multi_threading")

merge_option = parser.add_argument_group(title = "Merge report options")
merge_option.add_argument("-m",'--merge',
							metavar=("<OTHER_REPORT>", "<GENERIC_REPORT>"),
							nargs=2,
							dest="merge",
							default=None,
							help="Merge the generic .report and the other systems together without add a systems generic that is in the OTHER_REPORT. Write the new .report in GENERIC_REPORT directory")

stat_option = parser.add_argument_group(title = "Count and statistics options")
stat_option.add_argument("-stats",'--statistics',
							metavar=("<ANNOTATION_TABLE>"),
							dest="stats",
							default=None,
							help="Do count and statistics of the systems found. Need a annotation table to know the kingdom and phylum of the species where the systems is found. Need a annotation table")
stat_option.add_argument("-wanted",'--phylum_wanted',
							metavar=("<PHYLUM_LIST>"),
							dest="wanted",
							default=None,
							help="List of all the phylum that we want information about the number of systems found with one phylum by line. (default: All the phylum found)")
stat_option.add_argument("-so",'--stats_only',
							dest="stats_only",
							action='store_true',
							help="Remove all the steps of searching through database. [ONLY available if the previous analysis is done in the same folder]")


# TODO Peut être entre tous les systèmes prendre celui qui a le plus de gènes ? et non pas le premier. Mais du coup fait faire plus de calculs.

args = parser.parse_args()

if args.merge :
	fileAll = os.path.abspath(args.merge[0])
	fileGeneric = os.path.abspath(args.merge[1])

	fileWrite = os.path.join(os.path.dirname(os.path.dirname(fileGeneric)),"merge_macsyfinder.report")

	write_merge_file(fileGeneric, fileAll, fileWrite)

	print("\n#################")
	print("# File merged")
	print("#################\n")

	sys.exit(0)



if not args.output :
	OUTPUT = os.path.join(os.path.abspath(args.reportFile),"extraction_{}".format(time.strftime("%d_%m_%y")))
else :
	OUTPUT = args.output

create_folder(OUTPUT)

# XXX Creation of an information folder for each sequence remove or file generate (as cutoff, ...)
INFO=os.path.join(OUTPUT,"info_folder")
create_folder(INFO)
create_folder(os.path.join(INFO, "report_modif"))

if args.stats_only :
	if os.path.isfile(os.path.join(INFO, "report_modif" "detected.report")) :
		sys.exit("You need to specify the folder of a previous analysis with the option '-o' [PATH/TO/PREVIOUS/ANALYSIS]")
else :
	if args.seqData and args.reportFile :
		FASTA = args.seqData
		REPORT = args.reportFile
	else :
		sys.exit("sequence_extractor.py: error: the following arguments are required: -s/--seqdata, -r/--reportfile, -d/--defFile")

PROTEIN_FUNCTION = read_protein_function(args.defFile)

# XXX List of all the function name with .fasta add at the end
all_function = np.unique([function+".fasta" for function in PROTEIN_FUNCTION.values()]).tolist()

# XXX Si j'ai l'option verifié mise en place je demande les deux options.
if args.veriFile :
	veriFile, veriData = args.veriFile
	PATH_FASTA_VERIFIED = os.path.join(OUTPUT, "fasta_verified")
	list_file_verified = [os.path.join(PATH_FASTA_VERIFIED, function) for function in all_function]
	create_folder(PATH_FASTA_VERIFIED)
	report_df_verified = create_verified_fasta(list_file_verified, PROTEIN_FUNCTION, veriFile, veriData, INFO)

if not args.stats_only :
	# XXX Première liste de fichiers détectés
	print("\n#################")
	print("# Detected Fasta")
	print("#################\n")

	PATH_FASTA_DETECTED = os.path.join(OUTPUT, "fasta_detected", "raw")
	create_folder(PATH_FASTA_DETECTED)
	list_file_detected = [os.path.join(PATH_FASTA_DETECTED, function) for function in all_function]

	# XXX Je teste si je suis en multi_thread ou pas
	if args.number_proc == 1  :
		report_df_detected = find_in_fasta(FASTA, REPORT, list_file_detected, INFO, PROTEIN_FUNCTION)
	else :
		report_df_detected = find_in_fasta_multithreads(glob.glob(os.path.join(FASTA, "*")), REPORT, list_file_detected, INFO, PROTEIN_FUNCTION, int(args.number_proc))


if args.concat :
	# NOTE Je crée un dossier qui va contenir les fichiers détectés avec juste les séquences non identique au verifiées.
	PATH_FASTA_DETECTED_SELECTED = os.path.join(OUTPUT, "fasta_detected", "selected_concatenation")
	create_folder(PATH_FASTA_DETECTED_SELECTED)

	if args.cutoff :
		PATH_FASTA_CONCATENATED = os.path.join(OUTPUT, "fasta_concatenated", "raw")
	else :
		PATH_FASTA_CONCATENATED = os.path.join(OUTPUT, "fasta_concatenated")

	create_folder(PATH_FASTA_CONCATENATED)
	concatenate_detected_verified(all_function, PATH_FASTA_DETECTED, PATH_FASTA_VERIFIED, INFO, PATH_FASTA_CONCATENATED, PATH_FASTA_DETECTED_SELECTED)
	list_file_concatenated = [os.path.join(PATH_FASTA_CONCATENATED, function) for function in all_function]

# XXX Deuxieme liste de fichiers concaténés ou detectés après cutoff
if args.cutoff :
	if args.concat :
		PATH_FASTA_CONCATENATED_CUTOFF = os.path.join(OUTPUT, "fasta_concatenated", "cut_off")
		cut_seq_fasta_file(list_file_concatenated, PATH_FASTA_CONCATENATED_CUTOFF, INFO, file_cutoff=args.cutoff)
	else :
		PATH_FASTA_DETECTED_CUTOFF = os.path.join(OUTPUT, "fasta_detected", "cut_off")
		cut_seq_fasta_file(list_file_detected, PATH_FASTA_DETECTED_CUTOFF, INFO, file_cutoff=args.cutoff)

# XXX Dernière étape les stats

if args.stats_only :
	# XXX Debug si j'ai enregister les df_report dans des fichiers séparés.
	if os.path.isfile(os.path.join(INFO, "report_modif", "verified.report")) :
		report_df_verified = pd.read_table(os.path.join(INFO, "report_modif", "verified.report"), comment="#")
	report_df_detected = pd.read_table(os.path.join(INFO, "report_modif", "detected.report"), comment="#")

if args.stats :

	print("\n#################")
	print("# Count and Stats")
	print("#################\n")

	DICT_SYSTEMS = create_dict_system(PROTEIN_FUNCTION)

	# XXX Je crée le fichier ou je met mes systemes de mon analyse
	info_file = open(os.path.join(INFO,"systems_found.names"), "w")
	info_file.write("# {}\n".format("\t".join(["Species_Id","Replicon_Id","System_name","System_number","Proteins", "Kingdom", "Phylum", "Lineage"])))

	if os.path.isfile(os.path.join(INFO, "report_modif", "verified.report")) :
		# XXX Pour les verifiés
		info_file.write("## Verified systems\n")
		df_info_verified = set_df_info_system(report_df_verified, info_file, args.stats, DICT_SYSTEMS, "V")

	# XXX Pour les detectés
	info_file.write("## Detected systems\n")
	df_info_detected = set_df_info_system(report_df_detected, info_file, args.stats, DICT_SYSTEMS, "D")

	info_file.close()

	DICT_SPECIES = {kingdom:np.unique(df_info_detected.Phylum[df_info_detected.Kingdom == kingdom]) for kingdom in set(df_info_detected.Kingdom)}

	# BUG Je vais avoir besoin d'utiliser la fonction set_df_info_system() pour avoir des info dont j'ai besoin dans les fonctions suivantes

	if args.wanted :
		DICT_SPECIES_WANTED, LIST_WANTED = create_dict_wanted(args.wanted)
	else :
		LIST_WANTED = np.unique(df_info_detected.Phylum).tolist()
		DICT_SPECIES_WANTED = DICT_SPECIES

	PATH_TO_DATAFRAME = os.path.join(OUTPUT,"statistics")
	create_folder(os.path.join(PATH_TO_DATAFRAME, "xlsx"))
	create_folder(os.path.join(PATH_TO_DATAFRAME, "tsv"))
	create_folder(os.path.join(PATH_TO_DATAFRAME, "data_color"))
	create_folder(os.path.join(PATH_TO_DATAFRAME, "figure"))

	# XXX fichier ou je met des infos utiles
	out_file=open(os.path.join(PATH_TO_DATAFRAME, "stats.info"), 'w')

	print('Creating DataFrame ...')


	# XXX dataframe avec les counts pour chaques types de protéines je pense pas que je l'utilise après donc pas la peine de récupéré le dataframe
	# NOTE Ici pour les deux dataframes j'inclus les verifiés car ils font partis du DICT_INFO. Je pense que je vais enlever les verifiés car ils sont inutiles ici.
	df_count = count_all(df_info_detected, PROTEIN_FUNCTION, PATH_TO_DATAFRAME, DICT_SPECIES)

	# XXX Dataframe avec le compte pour tous les systemes donc utile pour la suite des "stats"
	df_systems = systems_count(df_info_detected, PATH_TO_DATAFRAME, LIST_WANTED, DICT_SPECIES_WANTED)

	print()
	print('Figure in process ...')

	# XXX premiere figure
	proportion_phylum(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, df_info_detected)

	# XXX Deuxieme figure
	proportion_systems(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, out_file)

	# XXX Troisieme figure
	proportion_proteobacteria(os.path.join(PATH_TO_DATAFRAME, "figure"), df_info_detected, df_systems.columns)

	# XXX Les dataframes en couleurs
	dataframe_color(os.path.join(PATH_TO_DATAFRAME, "data_color"), df_systems, DICT_SPECIES_WANTED, LIST_WANTED, args.stats, out_file)


	out_file.close()

	print()
	print("Done!")


print("\n#################")
print("# End")
print("#################\n")
