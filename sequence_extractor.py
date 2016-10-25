# -*- coding: utf-8 -*-
import argparse
from textwrap import dedent
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'library'))

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
							const=True,
							dest="cutoff",
							default=None,
							help="Option to remove sequences that are way to much longer or shorter beside of the other. (If there is no file, it will calculate the cutoff and generate a file)")
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
							action='store_true',
							dest="concat",
							help="Allow to concatenate detected sequences and verified sequences")

merge_option = parser.add_argument_group(title = "Merge report options")
merge_option.add_argument("-m",'--merge',
							metavar=("<OTHER_REPORT>", "<GENERIQUE_REPORT>"),
							nargs=2,
							dest="merge",
							default=None,
							help="Merge the generique .report and the other systems together without add a systems generique that is in the OTHER_REPORT. Write the new .report in GENERIQUE_REPORT directory")

stat_option = parser.add_argument_group(title = "Count and statistics options")
stat_option.add_argument("-stats",'--statistics',
							metavar=("<ANNOTATION_TABLE>"),
							nargs=1,
							dest="stats",
							default=None,
							help="Do count and statistics of the systems found. Need a annotation table to know the kingdom and phylum of the species where the systems is found")
stat_option.add_argument("-wanted",'--phylum_wanted',
							metavar=("<PHYLUM_LIST>"),
							nargs=1,
							dest="wanted",
							default=None,
							help="List of all the phylum that we want information about the number of systems found with one phylum by line. (default: All the phylum found)")


args = parser.parse_args()

if args.merge :
	fileAll = os.path.abspath(args.merge[0])
	fileGenerique = os.path.abspath(args.merge[1])

	fileWrite = os.path.join(os.path.dirname(fileGenerique),"merge_macsyfinder.report")

	write_merge_file(fileGenerique, fileAll, fileWrite)

	print("\n#################")
	print("# File merged")
	print("#################\n")

	sys.exit(0)



if not args.prefix :
	PREFIX = os.path.join(os.path.abspath(args.reportFile),"extraction_{}".format(time.strftime("%d_%m_%y")))
else :
	PREFIX = args.prefix

create_folder(PREFIX)

# XXX Creation of an information folder for each sequence remove or file generate (as cutoff, ...)
INFO=os.path.join(PREFIX,"info_folder")
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
list_file = robjects.StrVector(all_function)

# XXX Je crée le fichier ou je met mes systemes de mon analyse
info_file = open(os.path.join(INFO,"systems_found.names"), "w")
info_file.write("# Species_name\tSystem_name\tSystem_number\tProteins\n")

# XXX Si j'ai l'option verifié mise en place je demande les deux options.
if args.veriFile or args.veriData :
	if not (args.veriFile and args.veriData):
		parser.error("you MUST provided a verified fasta file and a annotation data file. If you want verified fasta")
	else :
		info_file.write("## Verified systems")
		PATH_FASTA_VERIFIED = os.path.join(PREFIX, "fasta_verified", "raw")
		create_folder(PATH_FASTA_VERIFIED)
		create_verified_fasta(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep='/'), PROTEIN_FUNCTION, args.veriFile, args.veriData)
		PATH_FASTA_RENAME = os.path.join(PREFIX, "fasta_verified", "rename")

		# NOTE Pas besoin de récupérer le dictionnaire et la liste des systèmes verifiés car je ne veux pas faire de stats dessus. Et qu'ici je ne peux pas séparer les detectés des verifiés qui sont identiques.
		rename_name_gene(robjects.r['paste'](PATH_FASTA_VERIFIED, list_file, sep='/'), PATH_FASTA_RENAME, info_file, DICT_SYSTEMS, DICT_INFO_VERIFIED)
		PATH_FASTA_VERIFIED = PATH_FASTA_RENAME


# XXX Première liste de fichiers détectés
info_file.write("## Detected systems")
PATH_FASTA_DETECTED = os.path.join(PREFIX, "fasta_detected", "raw")
create_folder(PATH_FASTA_DETECTED)
list_file_detected = robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep='/')

find_in_fasta(FASTA, REPORT, list_file_detected, INFO, PROTEIN_FUNCTION)

# XXX Deuxième liste de fichiers détectés après que tous les nom soit renomé
PATH_FASTA_RENAME = os.path.join(PREFIX, "fasta_detected", "rename")
DICT_INFO_DETECTED, list_systems_detected = rename_name_gene(list_file_detected, PATH_FASTA_RENAME, info_file, DICT_SYSTEMS, DICT_INFO_DETECTED)
PATH_FASTA_DETECTED = PATH_FASTA_RENAME
list_file_detected = robjects.r['paste'](PATH_FASTA_DETECTED, list_file, sep='/')

info_file.close()

if args.concat :
    # NOTE Je crée un dossier qui va contenir les fichiers détectés avec juste les séquences non identique au verifiées.
    create_folder(os.path.join(PATH_FASTA_DETECTED, "selected_concatenation"))

	if args.cutoff :
		PATH_FASTA_CONCATENATED = os.path.join(PREFIX, "fasta_concatenated", "raw")
	else :
		PATH_FASTA_CONCATENATED = os.path.join(PREFIX, "fasta_concatenated")

	create_folder(PATH_FASTA_CONCATENATED)
	concatenate_detected_verified(list_file, PATH_FASTA_DETECTED, PATH_FASTA_VERIFIED, INFO, PATH_FASTA_CONCATENATED)
	list_file_concatenated = robjects.r['paste'](PATH_FASTA_CONCATENATED, list_file, sep='/')

# XXX Deuxieme liste de fichiers concaténés ou detectés après cutoff
if args.cutoff :
	if args.concat :
		PATH_FASTA_CONCATENATED_CUTOFF = os.path.join(PREFIX, "fasta_concatenated", "cut_off")
		cut_seq_fasta_file(list_file_concatenated, PATH_FASTA_CONCATENATED_CUTOFF, INFO, file_cutoff=args.cutoff)
	else :
		PATH_FASTA_DETECTED_CUTOFF = os.path.join(PREFIX, "fasta_detected", "cut_off")
		cut_seq_fasta_file(list_file_detected, PATH_FASTA_DETECTED_CUTOFF, INFO, file_cutoff=args.cutoff)

# XXX Dernière étape les stats
if args.stats :

	print("\n#################")
	print("# Count and Stats")
	print("#################\n")

	INFO_STATS = np.genfromtxt(args.stats, dtype=str, delimiter="\t", comments="##")
	DICT_SPECIES = {kingdom:np.unique(info_tab[info_tab[:,2] == kingdom,3]) for kingdom in np.unique(info_tab[:,2])}

	if args.wanted :
		LIST_WANTED = read_list_wanted(args.wanted)
		DICT_SPECIES_WANTED = create_dict_wanted()
	else :
		LIST_WANTED = np.unique(INFO_STATS[:,3]).tolist()
		DICT_SPECIES_WANTED = DICT_SPECIES

	PATH_TO_DATAFRAME = os.path.join(PREFIX,"statistics")
	create_folder(os.path.join(PATH_TO_DATAFRAME, "xlsx"))
	create_folder(os.path.join(PATH_TO_DATAFRAME, "tsv"))
	create_folder(os.path.join(PATH_TO_DATAFRAME, "data_color"))
	create_folder(os.path.join(PATH_TO_DATAFRAME, "figure"))

	# XXX fichier ou je met des infos utiles
	out_file=open(PATH_TO_DATAFRAME, "stats.info")

    print('Creating DataFrame ...')

	# XXX dataframe avec les counts pour chaques types de protéines je pense pas que je l'utilise après donc pas la peine de récupéré le dataframe
	# NOTE Ici pour les deux dataframes j'inclus les verifiés car ils font partis du DICT_INFO. Je pense que je vais enlever les verifiés car ils sont inutiles ici.
	df_count = count_all(INFO_STATS, DICT_INFO, PROTEIN_FUNCTION, list_file_detected, PATH_TO_DATAFRAME, DICT_SPECIES)

	# XXX Dataframe avec le compte pour tous les systemes donc utile pour la suite des "stats"
	df_systems = systems_count(INFO_STATS, DICT_INFO, PROTEIN_FUNCTION, list_file_detected, PATH_TO_DATAFRAME, LIST_WANTED, DICT_SPECIES_WANTED)

    print('Figure in process ...')
	
	# XXX premiere figure
	proportion_phylum(os.path.join(PATH_TO_DATAFRAME, "figure"), df_count, LIST_SYSTEMS)

	# XXX Deuxieme figure
	proportion_systems(os.path.join(PATH_TO_DATAFRAME, "figure"), df_systems, out_file)

	# XXX Troisieme figure
	proportion_proteobacteria(os.path.join(PATH_TO_DATAFRAME, "figure"), df_systems)

	# XXX Les dataframes en couleurs
	dataframe_color(os.path.join(PATH_TO_DATAFRAME, "data_color"), df_systems, DICT_SPECIES_WANTED, LIST_WANTED)


	out_file.close()

    print("Done!")


print("\n#################")
print("# End")
print("#################\n")
