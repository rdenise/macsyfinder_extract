# -*- coding: utf-8 -*-
import argparse
from textwrap import dedent
import sys, os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'library'))

from set_params import *
from fasta_creation import *
from merge_generic_all import *
from statistique_count import *
from modification_report import *

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
							metavar="<file>",
							dest="defFile",
							help="Tabular file with the the function name in first column and name of the protein in the other")
general_option.add_argument("-o",'--output',
 							default=None,
							dest="output",
							metavar='<OUTPUT>',
							help="Using <OUTPUT> for output files (default: reportFile directory)")
general_option.add_argument("-annot",'--annotation',
							metavar=("<ANNOTATION_TABLE>", "<SYSTEMES_DISTANCE>"),
							nargs=2,
							dest="annotation",
							default=None,
							help="Use a annotation table to know the kingdom and phylum of the species where the systems is found and a tabulate file that contain the name of the system in fisrt columns and the maximum distance in the second. Write a report table with more information on it")
general_option.add_argument("-c",'--config_file',
							metavar="<file>",
							dest="config_file",
							help="The file .conf generate by macsyfinder")

extraction_option = parser.add_argument_group(title = "Extraction options")
extraction_option.add_argument("-cut",'--cutoff',
							metavar="<CUTOFF_FILE>",
							nargs='?',
							const=True,
							dest="cutoff",
							default=None,
							help="Option to remove sequences that are way to much longer or shorter beside of the other. (If there is no file, it will calculate the cutoff and generate a file)")
extraction_option.add_argument("-veriFile",'--validatedFile',
							metavar=("<VALIDATED_FASTA_FILE>", "<VALIDATED_DATA>"),
							nargs=2,
							dest="veriFile",
							default=None,
							help="File with the sequence in fasta of the validated file and the file with the information about the validated systems with this information : #SeqID	Gene	System	SystID	Family	Note")
extraction_option.add_argument("-conc",'--concatenate',
							action='store_true',
							dest="concat",
							help="Allow to concatenate detected sequences and validated sequences")
extraction_option.add_argument("-w", "--worker_proc",
                            metavar="<NUMBER_OF_THREADS>",
                            dest="number_proc",
                            default=1,
                            help="Number of processor you want to used in multi_threading")
extraction_option.add_argument("-rm_poor",'--remove_poor',
							metavar="<FUNCTION>",
							nargs='+',
							dest="remove_poor",
							default=None,
							help="Option to remove systems that do not have at least the proteins with the function gave in arguments")


stat_option = parser.add_argument_group(title = "Count and statistics options")
stat_option.add_argument("-stats",'--statistics',
							dest="stats",
							action='store_true',
							help="Do count and statistics of the systems found.")
stat_option.add_argument("-wanted",'--phylum_wanted',
							metavar=("<PHYLUM_LIST>"),
							dest="wanted",
							default=None,
							help="List of all the phylum that we want information about the number of systems found with one phylum by line. (default: All the phylum found)")
stat_option.add_argument("-so",'--stats_only',
							dest="stats_only",
							action='store_true',
							help="Remove all the steps of searching through database. [ONLY available if the previous analysis is done in the same folder]")


subparsers = parser.add_subparsers()
merge_parser = subparsers.add_parser('merge', help='Merge report command : Merge the generic .report and the other systems together without add a systems generic that is in the OTHER_REPORT. Write the new .report in GENERIC_REPORT directory')

merge_parser.add_argument("-or",'--other_report',
							metavar=("<OTHER_REPORT>"),
							dest="other_report",
							default=None,
							required=True,
							help="File with the 'classic' systems")
merge_parser.add_argument("-ar",'--archaea_report',
							metavar=("<ARCHAEA_REPORT>"),
							dest="archaea_report",
							default=None,
							help="File with the 'generic' systems")
merge_parser.add_argument("-gr",'--generic_report',
							metavar=("<GENERIC_REPORT>"),
							dest="generic_report",
							default=None,
							help="File with the 'generic' systems")

# TODO Peut être entre tous les systèmes prendre celui qui a le plus de gènes ? et non pas le premier. Mais du coup fait faire plus de calculs.

args = parser.parse_args()

if "other_report" in args :

	fileAll = os.path.abspath(args.other_report)

	if "archaea_report" not in args and "generic_report" not in args :
		sys.exit(dedent("""usage: sequence_extractor.py merge -or <OTHER_REPORT>,
		[-ar <ARCHAEA_REPORT>],
		[-gr >GENERIC_REPORT>]
		sequence_extractor.py: error: one of this arguments are required: -ar/--archaea_report, -gr/--generic_report"""))


	if "archaea_report" in args and args.archaea_report:

		print("----------------------")
		print("| Merging with archeae")
		print("----------------------")

		fileArchaea = os.path.abspath(args.archaea_report)
		fileWrite = os.path.join(os.path.dirname(os.path.dirname(fileArchaea)),"merge_macsyfinder_archaea.report")

		write_merge_file(fileArchaea, fileAll, fileWrite)
		fileAll = fileWrite

	if "generic_report" in args and args.generic_report :

		print("----------------------")
		print("| Merging with generic")
		print("----------------------")

		fileGeneric = os.path.abspath(args.generic_report)
		fileWrite = os.path.join(os.path.dirname(os.path.dirname(fileGeneric)),"merge_macsyfinder.report")

		write_merge_file(fileGeneric, fileAll, fileWrite)


	print("\n#################")
	print("# File merged")
	print("#################\n")

	sys.exit(0)

if not args.defFile or not args.annotation or not args.config_file :
	sys.exit(dedent("""usage: sequence_extractor.py [-h] [-s <files>] [-r <file>] -d <file>
                             [-o <OUTPUT>] -annot <ANNOTATION_TABLE>
                             <SYSTEMES_DISTANCE> -c <file>
                             [-cut [<CUTOFF_FILE>]]
                             [-veriFile <VALIDATED_FASTA_FILE> <VALIDATED_DATA>]
                             [-conc] [-w <NUMBER_OF_THREADS>]
                             [-rm_poor <FUNCTION> [<FUNCTION> ...]]
                             [-m <OTHER_REPORT> <GENERIC_REPORT>] [-stats]
                             [-wanted <PHYLUM_LIST>] [-so]
                             {merge} ...
sequence_extractor.py: error: the following arguments are required: -d/--defFile, -annot/--annotation, -c/--config_file"""))


if not args.output :
	OUTPUT = os.path.join(os.path.abspath(args.reportFile),"extraction_{}".format(time.strftime("%Y%m%d")))
else :
	OUTPUT = args.output

create_folder(OUTPUT)


if args.annotation :
	ANNOTATION = args.annotation[0]
	FILE_DISTANCE = args.annotation[1]

# XXX Creation of an information folder for each sequence remove or file generate (as cutoff, ...)
INFO=os.path.join(OUTPUT,"info_folder")
create_folder(INFO)
create_folder(os.path.join(INFO, "report_modif"))

if args.stats_only and args.reportFile:
	if os.path.isfile(os.path.join(INFO, "report_modif" "detected.report")) :
		sys.exit("You need to specify the folder of a previous analysis with the option -o/--output [PATH/TO/PREVIOUS/ANALYSIS], the report file -r/--reportfile [REPORTFILE] and -annot/--annotation [ANNOTATION_TABLE]")
elif args.stats and not args.annotation :
	sys.exit("You need to -annot/--annotation [ANNOTATION_TABLE] option to do stats")
else :
	if args.seqData and args.reportFile:
		FASTA = args.seqData
	else :
		sys.exit("sequence_extractor.py: error: the following arguments are required: -s/--seqdata, -r/--reportfile, -d/--defFile, -annot/--annotation")

REPORT = args.reportFile
PROTEIN_FUNCTION = read_protein_function(args.defFile)

# XXX List of all the function name with .fasta add at the end
all_function = np.unique([function+".fasta" for function in PROTEIN_FUNCTION.values()]).tolist()

# XXX Si j'ai l'option verifié mise en place je demande les deux options.
if args.veriFile :
	veriFile, veriData = args.veriFile
	PATH_FASTA_VALIDATED = os.path.join(OUTPUT, "fasta_validated")
	list_file_validated = [os.path.join(PATH_FASTA_VALIDATED, function) for function in all_function]
	create_folder(PATH_FASTA_VALIDATED)
	report_df_validated = create_validated_fasta(list_file_validated, PROTEIN_FUNCTION, veriFile, veriData, INFO)

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
		report_df_detected = find_in_fasta(FASTA, REPORT, list_file_detected, INFO, PROTEIN_FUNCTION, args.config_file)
	else :
		report_df_detected = find_in_fasta_multithreads(glob.glob(os.path.join(FASTA, "*")), REPORT, list_file_detected, INFO, PROTEIN_FUNCTION, int(args.number_proc), args.config_file)

	shutil.rmtree("tmp_{}_fasta".format(time.strftime("%Y%m%d")), ignore_errors=True)
	shutil.rmtree("tmp_{}_fasta".format(time.strftime("%Y%m%d")), ignore_errors=True)

# XXX Permet d'avoir un report plus complet avec plus d'information sur mes systemes seulement
if not args.stats_only and args.annotation:

	print("\n########################")
	print("# Create annotation files")
	print("########################\n")


	DICT_SYSTEMS = create_dict_system(PROTEIN_FUNCTION)


	print("\n------------------------")
	print("| Create found.names file")
	print("------------------------\n")

	# XXX Je crée le fichier ou je met mes systemes de mon analyse
	info_file = open(os.path.join(INFO,"systems_found.names"), "w")
	info_file.write("# {}\n".format("\t".join(["Species_Id","Replicon_Id","System_name","System_status","System_number","Proteins", "Kingdom", "Phylum", "Lineage"])))

	if os.path.isfile(os.path.join(INFO, "report_modif", "validated_tmp.report")) :
		# XXX Pour les verifiés
		info_file.write("## Validated systems\n")
		df_info_validated = set_df_info_system(report_df_validated, info_file, ANNOTATION, DICT_SYSTEMS, "V")

	# XXX Pour les detectés
	info_file.write("## Detected systems\n")
	df_info_detected = set_df_info_system(report_df_detected, info_file, ANNOTATION, DICT_SYSTEMS, "D")

	info_file.close()

	print()
	print("Done!")
	print()

	distance_max = create_dict_distance(FILE_DISTANCE)

    # XXX Creation of the modified report with information about the if the gene is in tandem, loner, multi copy, the kingdom of the species, the phylum...
	if os.path.isfile(os.path.join(INFO, "report_modif", "validated_tmp.report")) :
		merge_detected_validated_report(os.path.join(INFO, "report_modif", "validated_tmp.report"), REPORT, os.path.join(INFO, "report_modif", "detected.report"), os.path.join(INFO,"systems_found.names"), pd.read_table(ANNOTATION, index_col=0), PROTEIN_FUNCTION, os.path.join(INFO, "report_modif"), distance_max)
	else :
		names_dataframe = ['Hit_Id','Replicon_name','Position','Sequence_length','Gene','Reference_system','Predicted_system','System_Id','System_status','Gene_status','i-evalue','Score','Profile_coverage','Sequence_coverage','Begin_match','End_match']
		REPORT_df = pd.read_table(REPORT, names=names_dataframe)
		REPORT_df_modif = pd.read_table(os.path.join(INFO, "report_modif", "detected.report"))

		fill_reportdf_detected(REPORT_df, PROTEIN_FUNCTION, REPORT_df_modif, pd.read_table(ANNOTATION, index_col=0), os.path.join(INFO, "report_modif"), distance_max, names_dataframe)


if args.concat :
	# NOTE Je crée un dossier qui va contenir les fichiers détectés avec juste les séquences non identique au verifiées.
	PATH_FASTA_DETECTED_SELECTED = os.path.join(OUTPUT, "fasta_detected", "selected_concatenation")
	create_folder(PATH_FASTA_DETECTED_SELECTED)

	if args.cutoff or args.remove_poor:
		PATH_FASTA_CONCATENATED = os.path.join(OUTPUT, "fasta_concatenated", "raw")
	else :
		PATH_FASTA_CONCATENATED = os.path.join(OUTPUT, "fasta_concatenated")

	create_folder(PATH_FASTA_CONCATENATED)
	concatenate_detected_validated(all_function, PATH_FASTA_DETECTED, PATH_FASTA_VALIDATED, INFO, PATH_FASTA_CONCATENATED, PATH_FASTA_DETECTED_SELECTED)
	list_file_concatenated = [os.path.join(PATH_FASTA_CONCATENATED, function) for function in all_function]

	# NOTE je met un nouveau path pour les fichier detectés si jamais je dois revenir sur ces fichiers
	PATH_FASTA_DETECTED = PATH_FASTA_DETECTED_SELECTED
	list_file_detected = [os.path.join(PATH_FASTA_DETECTED, function) for function in all_function]


# XXX Permet d'enlevé les systemes sans ATPase et/ou IMplatform
if args.remove_poor :
	# NOTE Je crée un dossier qui va contenir les fichiers détectés avec juste les séquences selectionés
	PATH_FASTA_DETECTED_SELECTED = os.path.join(OUTPUT, "fasta_detected", "selected_remove_poor")
	create_folder(PATH_FASTA_DETECTED_SELECTED)

	dict_all_function = read_protein_function_reverse(args.defFile)
	DICT_IMPORTANT_PROTEINS = {function:dict_all_function[function] for function in dict_all_function if function in args.remove_poor}

	report_df_full = pd.read_table(os.path.join(INFO, "report_modif","report_modif_annotation_full.report"))

	if args.concat :
		PATH_FASTA_CONCATENATED_REMOVE_POOR = os.path.join(OUTPUT, "fasta_concatenated", "selected_remove_poor")
		create_folder(PATH_FASTA_CONCATENATED_REMOVE_POOR)
		concatenate_reduce(list_file_concatenated, PATH_FASTA_CONCATENATED_REMOVE_POOR, report_df_full, DICT_IMPORTANT_PROTEINS, INFO, liste_detected_file = list_file_detected, path_detected = PATH_FASTA_DETECTED_SELECTED)

		# NOTE nouveau path concatenated
		PATH_FASTA_CONCATENATED = PATH_FASTA_CONCATENATED_REMOVE_POOR
		list_file_concatenated = [os.path.join(PATH_FASTA_CONCATENATED, function) for function in all_function]
	else :
		concatenate_reduce(list_file_detected, PATH_FASTA_DETECTED_SELECTED, report_df_full, DICT_IMPORTANT_PROTEINS, INFO)

	# NOTE nouveau path detected
	PATH_FASTA_DETECTED = PATH_FASTA_DETECTED_SELECTED
	list_file_detected = [os.path.join(PATH_FASTA_DETECTED, function) for function in all_function]


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
	if os.path.isfile(os.path.join(INFO, "report_modif", "validated.report")) :
		report_df_validated = pd.read_table(os.path.join(INFO, "report_modif", "validated.report"), comment="#")
	report_df_detected = pd.read_table(os.path.join(INFO, "report_modif", "detected.report"), comment="#")

if args.stats or args.stats_only:

	print("\n#################")
	print("# Count and Stats")
	print("#################\n")


	DICT_SPECIES = {kingdom:np.unique(df_info_detected.Phylum[df_info_detected.Kingdom == kingdom]).tolist() for kingdom in set(df_info_detected.Kingdom)}
	DICT_SPECIES['Archaea'].append('Other')
	DICT_SPECIES['Bacteria'].append('Other')

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
	print("\n- Creating count DataFrame for the proteines...", end='\t')
	df_count = count_all(df_info_detected, PROTEIN_FUNCTION, PATH_TO_DATAFRAME, DICT_SPECIES)
	print("Done !")

	# XXX Dataframe avec le compte pour tous les systemes donc utile pour la suite des "stats"
	print("\n- Creating count DataFrame for the systems...", end='\t')
	df_systems = systems_count(df_info_detected, PATH_TO_DATAFRAME, LIST_WANTED, DICT_SPECIES_WANTED)
	print("Done !")

	print()
	print('Figure in process ...')

	# XXX premiere figure
	print("\n- Creating the figure with the proportion of phylum...", end='\t')
	proportion_phylum(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, df_info_detected)
	print("Done !")

	# XXX Deuxieme figure
	print("\n- Creating the figure with the proportion of systems...", end='\t')
	proportion_systems(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, out_file)
	print("Done !")

	# XXX Troisieme figure
	print("\n- Creating the figure with the proportion of Proteobacteria...", end='\t')
	proportion_proteobacteria(os.path.join(PATH_TO_DATAFRAME, "figure"), df_info_detected, df_systems.columns)
	print("Done !")

	# XXX Les dataframes en couleurs
	print("\n- Creating the color in the dataframe of count of systems...", end='\n')
	dataframe_color(os.path.join(PATH_TO_DATAFRAME, "data_color"), df_systems, DICT_SPECIES_WANTED, LIST_WANTED, ANNOTATION, out_file)
	print("Done !")

	# XXX Information sur les detection de systems validés
	if args.veriFile :

		print("\n------------------")
		print("| Validated System Stats")
		print("------------------\n")

		print("\n- Creating the figure with the identification of validated systems in the detected systems...", end='\t')

		validated_stats(veriData, REPORT, args.config_file, os.path.join(PATH_TO_DATAFRAME, "figure"), INFO)
		print("Done !")

	out_file.close()

	# XXX Figure sur la distribution de toutes les protéines des systèmes
	print("\n- Creating the figure of the distribution of the kind of proteins in the systems...", end='\t')
	df_found = read_systems_found(os.path.join(INFO,"systems_found.names"))
	dict_protein = set_dict_protein(args.defFile)
	plot_protein_distribution(df_found, os.path.join(PATH_TO_DATAFRAME, "figure"), dict_protein)
	print("Done !")

	print()
	print("Done!")


print("\n#################")
print("# End")
print("#################\n")
