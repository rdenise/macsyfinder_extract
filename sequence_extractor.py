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
 							required=True,
							metavar="<file(s)>",
							dest="seqData",
							help="File with all the fasta sequences used for macsyfinder analysis. If you want to multiprocessing put the name of the folder split by species")
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

FASTA = args.seqData
REPORT = args.reportFile

PROTEIN_FUNCTION = read_protein_function(args.defFile)
DICT_SYSTEMS = create_dict_system(PROTEIN_FUNCTION)

# XXX Dictionnaire qui va recuperer les associations nom d'espèce : nom de systemes pour éviter les doublons
DICT_INFO_VERIFIED = {}
DICT_INFO_DETECTED = {}

# XXX List of all the function name with .fasta add at the end
all_function = np.unique([function+".fasta" for function in PROTEIN_FUNCTION.values()]).tolist()
# BUG Plus besoin non plus car j'ai plus besoin de robjects
list_file = robjects.StrVector(all_function)

# XXX Je crée le fichier ou je met mes systemes de mon analyse
info_file = open(os.path.join(INFO,"systems_found.names"), "w")
# BUG Enlevéer si nouvelle fonction en stats ?
info_file.write("# Species_name\tSystem_name\tSystem_number\tProteins\n")

# XXX Si j'ai l'option verifié mise en place je demande les deux options.
if args.veriFile :
	veriFile, veriData = args.veriFile
    # BUG Pareil ici ?
	info_file.write("## Verified systems\n")
	PATH_FASTA_VERIFIED = os.path.join(OUTPUT, "fasta_verified", "raw")
	create_folder(PATH_FASTA_VERIFIED)
	# BUG L'appel devient pour l'instant :
	create_verified_fasta(os.path.join(PATH_FASTA_VERIFIED, all_function), PROTEIN_FUNCTION, veriFile, veriData, INFO)
	#create_verified_fasta(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep='/'), PROTEIN_FUNCTION, veriFile, veriData, INFO)
	# BUG Plus besoin de ça il faut que je renome dans la fonction précédente
	PATH_FASTA_RENAME = os.path.join(OUTPUT, "fasta_verified", "rename")

	# NOTE Pas besoin de récupérer le dictionnaire et la liste des systèmes verifiés car je ne veux pas faire de stats dessus. Et qu'ici je ne peux pas séparer les detectés des verifiés qui sont identiques.
	# BUG Fonction a disparue donc plus besoin
	rename_name_gene(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep='/'), PATH_FASTA_RENAME, info_file, DICT_SYSTEMS, DICT_INFO_VERIFIED)
	PATH_FASTA_VERIFIED = PATH_FASTA_RENAME


# XXX Première liste de fichiers détectés
print("\n#################")
print("# Detected Fasta")
print("#################\n")

info_file.write("## Detected systems\n")
PATH_FASTA_DETECTED = os.path.join(OUTPUT, "fasta_detected", "raw")
create_folder(PATH_FASTA_DETECTED)
list_file_detected = os.path.join(PATH_FASTA_DETECTED, all_function)

# XXX Je teste si je suis en multi_thread ou pas
if args.number_proc == 1  :
	find_in_fasta(FASTA, REPORT, list_file_detected, INFO, PROTEIN_FUNCTION)
else :
	find_in_fasta_multithreads(glob.glob(os.path.join(FASTA, "*")), REPORT, list_file_detected, INFO, PROTEIN_FUNCTION, int(args.number_proc))

# XXX Deuxième liste de fichiers détectés après que tous les nom soit renomé
PATH_FASTA_RENAME = os.path.join(OUTPUT, "fasta_detected", "rename")
# BUG car je n'ai plus rename donc je ne renvoie rien les informations sont maintenant retrouner par la fonction set_df_info_system()
DICT_INFO_DETECTED, list_systems_detected = rename_name_gene(list_file_detected, PATH_FASTA_RENAME, info_file, DICT_SYSTEMS, DICT_INFO_DETECTED)
PATH_FASTA_DETECTED = PATH_FASTA_RENAME
list_file_detected = robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep='/')

info_file.close()

# XXX Creation de la tble de translation renomée.
# BUG Plus besoin les sequence sont bien nommé dès le début
rename_seq_translation_table(INFO)

if args.concat :
	# NOTE Je crée un dossier qui va contenir les fichiers détectés avec juste les séquences non identique au verifiées.
	PATH_FASTA_DETECTED_SELECTED = os.path.join(OUTPUT, "fasta_detected", "selected_concatenation")
	create_folder(PATH_FASTA_DETECTED_SELECTED)

	if args.cutoff :
		PATH_FASTA_CONCATENATED = os.path.join(OUTPUT, "fasta_concatenated", "raw")
	else :
		PATH_FASTA_CONCATENATED = os.path.join(OUTPUT, "fasta_concatenated")

	create_folder(PATH_FASTA_CONCATENATED)
	concatenate_detected_verified(list_file, PATH_FASTA_DETECTED, PATH_FASTA_VERIFIED, INFO, PATH_FASTA_CONCATENATED, PATH_FASTA_DETECTED_SELECTED)
	list_file_concatenated = robjects.r['paste'](PATH_FASTA_CONCATENATED, list_file, sep='/')

# XXX Deuxieme liste de fichiers concaténés ou detectés après cutoff
if args.cutoff :
	if args.concat :
		PATH_FASTA_CONCATENATED_CUTOFF = os.path.join(OUTPUT, "fasta_concatenated", "cut_off")
		cut_seq_fasta_file(list_file_concatenated, PATH_FASTA_CONCATENATED_CUTOFF, INFO, file_cutoff=args.cutoff)
	else :
		PATH_FASTA_DETECTED_CUTOFF = os.path.join(OUTPUT, "fasta_detected", "cut_off")
		cut_seq_fasta_file(list_file_detected, PATH_FASTA_DETECTED_CUTOFF, INFO, file_cutoff=args.cutoff)

# XXX Dernière étape les stats
if args.stats :

	print("\n#################")
	print("# Count and Stats")
	print("#################\n")

	PROTEIN_FUNCTION_ONLY = {"_".join(key.split("_")[1:]):value for key, value in PROTEIN_FUNCTION.items()}
	INFO_STATS = np.genfromtxt(args.stats, dtype=str, delimiter="\t", comments="##")
	DICT_SPECIES = {kingdom:np.unique(INFO_STATS[INFO_STATS[:,3] == kingdom,4]) for kingdom in np.unique(INFO_STATS[:,3]).tolist()}

	# BUG Je vais avoir besoin d'utiliser la fonction set_df_info_system() pour avoir des info dont j'ai besoin dans les fonctions suivantes

	if args.wanted :
		DICT_SPECIES_WANTED, LIST_WANTED = create_dict_wanted(args.wanted)
	else :
		LIST_WANTED = np.unique(INFO_STATS[:,4]).tolist()
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
	df_count = count_all(INFO_STATS, DICT_INFO_DETECTED, PROTEIN_FUNCTION_ONLY, list_systems_detected, PATH_TO_DATAFRAME, DICT_SPECIES)

	# XXX Dataframe avec le compte pour tous les systemes donc utile pour la suite des "stats"
	df_systems = systems_count(INFO_STATS, DICT_INFO_DETECTED, PROTEIN_FUNCTION_ONLY, list_systems_detected, PATH_TO_DATAFRAME, LIST_WANTED, DICT_SPECIES_WANTED)

	print()
	print('Figure in process ...')

	# XXX premiere figure
	proportion_phylum(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, list_systems_detected)

	# XXX Deuxieme figure
	proportion_systems(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, out_file)

	# XXX Troisieme figure
	proportion_proteobacteria(os.path.join(PATH_TO_DATAFRAME, "figure"), df_systems)

	# XXX Les dataframes en couleurs
	dataframe_color(os.path.join(PATH_TO_DATAFRAME, "data_color"), df_systems, DICT_SPECIES_WANTED, LIST_WANTED, INFO_STATS, out_file)


	out_file.close()

	print()
	print("Done!")


print("\n#################")
print("# End")
print("#################\n")
