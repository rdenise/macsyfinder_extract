# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Library
##
##########################################################################################
##########################################################################################

import numpy as np
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import json
from weasyprint import HTML
from Bio import SeqIO

from set_params import *

##########################################################################################
##########################################################################################
##
##								Figure Style
##
##########################################################################################
##########################################################################################

sns.set_style("ticks")
sns.set_context("talk")

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def set_df_info_system(report_df, w_file, INFO_TAB, DICT_SYSTEMS, status) :

	"""
	Function that create a dataframe and set all the information about the systeme and write it in a file

	:param report_file: report dataframe (or report like dataframe) that contain at least the column "System_Id" and "Gene"
	:type: pandas.DataFrame
	:param w_file: handle of the file we want to write (that will contain the information created by this function)
	:type: open file
	:param INFO_TAB: name of the INFO_TAB file that contain all the information about the phylum and lineage of each prokaryote of gembases
	:type: str
	:param DICT_SYSTEMS: The dictionnary that contains the name of all the systems in key and the list of all the protein of this systems in values
	:type: dict
	:return: a dataframe with the information about ["Species_Id","Replicon_Id","System_name","System_number","Proteins", "Kingdom", "Phylum", "Lineage"]
	:rtype: pandas.Dataframe
	"""

	# NOTE FINIR DE LE METTRE EN  FONCTION done je pense

	value_counts_series = report_df.groupby("System_Id").Gene.value_counts()
	info_tab = pd.read_table(INFO_TAB, index_col=0, names=["Taxon_id", "Name", "Kingdom", "Phylum", "Lineage", "NC_ids"])

	df_info_system = pd.DataFrame(index=value_counts_series.index.levels[0], columns=["Species_Id","Replicon_Id","System_name", "System_status", "System_number","Proteins", "Kingdom", "Phylum", "Lineage"])
	for my_index in df_info_system.index :
		df_info_system.set_value(my_index, "Proteins", value_counts_series.loc[my_index].to_dict())

		Replicon, System, Number = my_index.split("_")
		if int(Number) == 1 :
		    df_info_system.loc[my_index, "System_number"] = "."
		else :
		    df_info_system.loc[my_index, "System_number"] = Number

		if System in ("Tcp", "R64", "Cof", "Bfp", "Lng"):
		    System = "T4bP"

		if System.lower() not in ["generic", "generique", "t4bp"] :
			protein_in_system = list(df_info_system.loc[my_index, "Proteins"].keys())
			key_system = protein_in_system[0].split("_")[0]
			theroric_protein_system = DICT_SYSTEMS[key_system]

			proteins_to_add = set(theroric_protein_system) - set(protein_in_system)
			df_info_system.loc[my_index, "Proteins"].update({key:0 for key in proteins_to_add})


		df_info_system.loc[my_index, "System_name"] = System
		df_info_system.loc[my_index, "Replicon_Id"] = Replicon
		df_info_system.loc[my_index, "Species_Id"] = ".".join(Replicon.split(".")[:-1])

		df_info_system.loc[my_index, "Kingdom"] = info_tab.loc[df_info_system.loc[my_index, "Species_Id"], "Kingdom"]
		df_info_system.loc[my_index, "Phylum"] = info_tab.loc[df_info_system.loc[my_index, "Species_Id"], "Phylum"]
		df_info_system.loc[my_index, "Lineage"] = info_tab.loc[df_info_system.loc[my_index, "Species_Id"], "Lineage"]

	df_info_system["System_status"] = status
	df_info_system.to_csv(w_file, sep="\t", index=False, header=False)
	return df_info_system

##########################################################################################
##########################################################################################


def tuple_like_all(species_dict, protein_list):

	"""
	Function that creating the multi index for pandas dataframe

	:param species_dict: dictionary that contain in key the name of the kingdom and in value all the phylum of the analysis
	:type: dict
	:param protein_list: name of all the proteins function of the analysis
	:type: list of str
	:return: the multi index for pandas dataframe
	:rtype: pandas.indexes.multi.MultiIndex
	"""

	list_finale = []
	for kingdom in ['Bacteria', 'Archaea'] :
		for phylum in species_dict[kingdom]:
			for protein in protein_list:
				list_finale.append((kingdom,phylum,protein,"Unique"))
				list_finale.append((kingdom,phylum,protein,"Total"))

			list_finale.append((kingdom,phylum,"Total_system", ""))

		list_finale.append((kingdom, "Total_system", "", ""))

	list_finale.append(("Summary_total","","",""))

	return pd.MultiIndex.from_tuples(list_finale)


##########################################################################################
##########################################################################################

def count_all(df_info_system, protein_function, PATH_TO_DATAFRAME, speciesDict) :

	"""
	Function that will read the df_info_system and count all information about system, protein, species...
	It will write a file with this information and return a dataframe.
	It's a count by protein so I'll have the number of ATPase in T2SS for Proteobacteria

	:param df_info_system: the dataframe with all the information about ["Species_Id","Replicon_Id","System_name","System_number","Proteins", "Kingdom", "Phylum", "Lineage"]
	:type: pandas.DataFrame
	:param protein_function: dictionnary with the name of all the protein studied
	:type: dict
	:param PATH_TO_DATAFRAME: absolute path to the results folder
	:type: str
	:param speciesDict: Dictionary that contain in key the kingdoms name and in value the list
	of all the phylum
	:type: dict
	:return: a dataframe with all this count of the systems
	:rtype: pandas.Dataframe
	"""

	LIST_SYSTEMS = np.unique(df_info_system.System_name).tolist()
	proteinlist = set(protein_function.values())

	miindex = tuple_like_all(speciesDict, proteinlist)
	df_count_system = pd.DataFrame(0, index=miindex, columns=LIST_SYSTEMS+['Total'])

	list_seen=[]

	for index, species_line in df_info_system.iterrows() :

		if species_line.Phylum in ['Archaea','Bacteria'] :
			species_line.loc["Phylum"] = 'Other'

		df_count_system.loc[(species_line.Kingdom,species_line.Phylum,"Total_system",""),species_line.System_name] += 1
		df_count_system.loc[(species_line.Kingdom,"Total_system","",""),species_line.System_name] += 1
		df_count_system.loc[("Summary_total","","",""),species_line.System_name] += 1

		for protein_name in species_line.Proteins :
			if species_line.Proteins[protein_name] != 0 :
				df_count_system.loc[(species_line.Kingdom,species_line.Phylum,protein_function[protein_name],'Unique'),species_line.System_name] += 1
			df_count_system.loc[(species_line.Kingdom,species_line.Phylum,protein_function[protein_name],'Total'),species_line.System_name] += species_line.Proteins[protein_name]

	df_count_system.Total = df_count_system.sum(axis=1)

	df_count_system = df_count_system[(df_count_system.Total > 0)]
	df_count_system.index.names = ["Kingdom", "Phylum", "Function", "Count"]
	df_count_system.to_csv(os.path.join(PATH_TO_DATAFRAME,"tsv","data_frame_count_all.tsv"), sep="\t", float_format="%.1f", index_label=False, encoding="utf-8" )
	df_count_system.to_excel(os.path.join(PATH_TO_DATAFRAME,"xlsx","data_frame_count_all.xlsx"), index_label=False)
	return df_count_system




##########################################################################################
##########################################################################################

def tuple_like_systems(species_dict):

	"""
	Function that creating the multi index for pandas dataframe

	:param species_dict: dictionary that contain in key the name of the kingdom and in value all the phylum of the analysis
	:type: dict
	:return: the multi index for pandas dataframe
	:rtype: pandas.indexes.multi.MultiIndex
	"""

	list_finale = []
	for kingdom in ['Bacteria', 'Archaea'] :
		for phylum in species_dict[kingdom]:
			list_finale.append((kingdom,phylum))

		list_finale.append((kingdom, "Other"))
		list_finale.append((kingdom, "Total_system"))

	list_finale.append(("Summary_total",""))

	return pd.MultiIndex.from_tuples(list_finale)


##########################################################################################
##########################################################################################

def systems_count(df_info_system, PATH_TO_DATAFRAME, list_wanted, speciesDict) :

	"""
	Function that will read the DICT_INFO and count all information about system, protein, species...
	It will write a file with this information and return a dataframe.
	It's a count by systems so I'll have the number of T2SS for Proteobacteria

	:param df_info_system: the dataframe with all the information about ["Species_Id","Replicon_Id","System_name","System_number","Proteins", "Kingdom", "Phylum", "Lineage"]
	:type: pandas.DataFrame
	:param PATH_TO_DATAFRAME: absolute path to the results folder
	:type: str
	:param list_wanted: list of all the phylum wanted
	:type: list of str
	:param speciesDict: Dictionary that contain in key the kingdoms name and in value the list
	of all the phylum for all the phylum wanted
	:type: dict
	:return: a dataframe with all this count of the systems
	:rtype: pandas.Dataframe
	"""

	LIST_SYSTEMS = np.unique(df_info_system.System_name).tolist()

	miindex = tuple_like_systems(speciesDict)
	df_count_system = pd.DataFrame(0, index=miindex, columns=LIST_SYSTEMS)

	# XXX on prend ici que le nombre d'espèce avec au moins 1 systemes (les doublons détectés pas comptés)

	for kingdom in speciesDict :
		for phylum in speciesDict[kingdom] :
			df_count_system.loc[(kingdom, phylum)] = df_info_system[((df_info_system.System_number == ".") | (df_info_system.System_name == "generic")) & (df_info_system.Lineage.str.contains(phylum))].System_name.value_counts()

		df_count_system.loc[(kingdom,"Other")] = df_info_system[((df_info_system.System_number == ".") | (df_info_system.System_name == "generic")) & ~(df_info_system.Lineage.str.contains("|".join(list_wanted)))  & (df_info_system.Kingdom == kingdom)].System_name.value_counts()
		df_count_system.loc[(kingdom,"Total_system")] = df_info_system[((df_info_system.System_number == ".") | (df_info_system.System_name == "generic")) & (df_info_system.Kingdom == kingdom)].System_name.value_counts()

	df_count_system.loc[("Summary_total","")] = df_info_system[((df_info_system.System_number == ".") | (df_info_system.System_name == "generic"))].System_name.value_counts()
	df_count_system.fillna(0, inplace=True)
	df_count_system = df_count_system.astype(int)

	df_count_system.index.names = ["Kingdom", "Phylum"]
	df_count_system.to_csv(os.path.join(PATH_TO_DATAFRAME,"tsv","data_frame_count_systems.tsv"), sep="\t", index_label=False, encoding="utf-8")
	df_count_system.to_excel(os.path.join(PATH_TO_DATAFRAME,"xlsx","data_frame_count_systems.xlsx"), index_label=False)
	return df_count_system

##########################################################################################
##########################################################################################

def proportion_phylum(PATH_TO_FIGURE, df, report_df):

	"""
	Function that plot the proportion of each phylum in all the systems found

	:param PATH_TO_FIGURE: The absolute path to the figure folder
	:type: str
	:param df: The dataframe with the counts with the detail of the phylogeny
	:type: pandas.Dataframe
	:param LIST_SYSTEMS: list with the name of all the species of the analysis (ex : T2SS, T4P ...)
 	:type: list of str
	:return: Nothing
	"""

	sub_df_for_plot=df.xs("Total_system", level="Function")
	name_legend = np.array([i[1] for i in sub_df_for_plot.index])

	# XXX Plot en lui même
	plot = (sub_df_for_plot.loc[:,np.unique(report_df.System_name)]/df.xs("Summary_total", level="Kingdom").loc[:,np.unique(report_df.System_name)].ix[0]).T.plot(kind="bar", stacked=True, cmap='Paired', rot=0)

	# XXX Mise en forme
	plt.legend(name_legend,loc='center left', bbox_to_anchor=(1.0, 0.5))
	plt.subplots_adjust(right=0.75)
	plt.ylabel("Proportion")
	plt.xlabel("System")
	plt.title("Proportion of Phylum in the studied system")
	ax_phylum = plt.gca()

	# XXX Sauvegarde du dictionnaire de couleurs pour une future utilisation qui peut être utilisé plus tard pour les couleurs
	# NOTE Je pense pas en avoir besoin
	#info_legend = ax_phylum.get_legend_handles_labels()
	#taille_infolegend =len(info_legend[1])
	#colors_dict = {info_legend[1][i]:info_legend[0][i][0].get_facecolor() for i in range(taille_infolegend)}
	#f = open(os.path.join(PATH_TO_FIGURE, ".color_dict.json"), 'w')
	#json.dump(colors_dict, f)
	#f.close()

	# XXX Sauvegarde de la figure
	plt.savefig(os.path.join(PATH_TO_FIGURE,"proportion_phylum.pdf"))
	plt.close('all')

	return

##########################################################################################
##########################################################################################

def proportion_systems(PATH_TO_FIGURE, df, w_file):

	"""
	Function that plot the proportion of each systems found

	:param PATH_TO_FIGURE: The absolute path to the figure folder
	:type: str
	:param df: The dataframe with the counts with the detail of the phylogeny
	:type: pandas.Dataframe
	:param w_file: Open file where the textual information will be write
 	:type: file
	:return: Nothing
	"""

	# XXX On met un tableau recapitulatif en texte
	w_file.write("Percentage within kingdoms\n")
	w_file.write("--------------------------\n\n")
	mini_tab = df.xs("Total_system", level="Phylum").iloc[:,:-1].div(df.xs("Total_system", level="Phylum").Total, axis='index')*100
	mini_tab.to_string(w_file)

	# XXX Et le même en figure
	ax = mini_tab.plot(kind='bar', rot=0)

	# XXX Mise en forme
	plt.xlabel("System")
	plt.ylabel("Percentage")
	ax.set_xticklabels(['Bacteria', 'Archaea'])
	plt.title("Percentage of studied system found in each kingdom clades studied")
	plt.savefig(os.path.join(PATH_TO_FIGURE,"proportion_found_each_kingdom.pdf"))
	plt.close('all')

	# XXX Proportion générale
	ax = (df.xs("Summary_total", level="Kingdom").iloc[:,:-1]/df.xs("Summary_total", level="Kingdom").Total.ix[0]).plot(kind='bar')

	# XXX Mise en forme
	plt.xlabel("System")
	plt.ylabel("Proportion")
	ax.set_xticklabels("")
	plt.title("Proportion of studied system found")
	plt.savefig(os.path.join(PATH_TO_FIGURE,"proportion_found.pdf"))
	plt.close('all')

	return

##########################################################################################
##########################################################################################

def proportion_proteobacteria(PATH_TO_FIGURE, df_found, LIST_SYSTEMS):

	"""
	Function that plot the proportion of each systems found

	:param PATH_TO_FIGURE: The absolute path to the figure folder
	:type: str
	:param df_found: The dataframe with the information about the  detected systems
	:type: pandas.Dataframe
	:param LIST_SYSTEMS: list of the name of all the systems of the analysis (ex : T2SS, T4P ...)
	:type: list of str
	:return: Nothing
	"""

	# XXX Creation du dataframe intermediare
	df_figure_2 = pd.DataFrame(0, index=['Proteobacteria (𝛼, 𝛽, 𝛾)','Rest'] ,columns=LIST_SYSTEMS)
	df_figure_2.loc['Proteobacteria (𝛼, 𝛽, 𝛾)'] = df_found[(df_found.Lineage.str.contains("Gammaproteobacteria")) | (df_found.Lineage.str.contains("Alphaproteobacteria")) | (df_found.Lineage.str.contains("Betaproteobacteria"))].System_name.value_counts()
	df_figure_2.loc['Rest'] = df_found[~((df_found.Lineage.str.contains("Gammaproteobacteria")) | (df_found.Lineage.str.contains("Alphaproteobacteria")) | (df_found.Lineage.str.contains("Betaproteobacteria")))].System_name.value_counts()
	df_figure_2.fillna(0, inplace=True)

	# XXX Le plot avec les couleurs choisies
	colors = [(0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0), (0.3997693305214246, 0.6478123867044262, 0.80273742044673246, 1.0)]
	plot = df_figure_2.T.plot(kind="bar", stacked=True, color=colors, rot=0, legend=False)
	sns.despine(offset=10, trim=True)

	# XXX Ajout des valeurs en haut du plot et on ne prend que les valeurs finale
	index_true_patches = int(len(plot.patches)/2)
	true_patches = plot.patches[index_true_patches:]
	for p in true_patches:
		x=p.get_bbox().get_points()[:,0]
		y=p.get_bbox().get_points()[1,1]
		plot.annotate('{}'.format(int(y)), (x.mean(), y),
				ha='center', va='bottom')

	plt.savefig(os.path.join(PATH_TO_FIGURE,"numbers_proteobacteria_rest.pdf"))
	plt.close('all')

	return

##########################################################################################
##########################################################################################

def do_gradient(data, count_df, color='red'):

    """
    Create a gradient in the background in a Series or DataFrame
	:param data: a Series from an apply of a DataFrame
	:type: pd.Series
	:param color: The spectrum of color we want a gradient
	:type: str
	:return: a Styler with the color in each cells
	:rtype: pd.DataFrame.Styler
    """
    cm = sns.light_palette(color, as_cmap=True)
    div_df = data.div(count_df.Count, axis=0)
    return ['background-color: #{:02X}{:02X}{:02X}'.format(*cm(value, bytes=1)[:-1]) for value in div_df]

##########################################################################################
##########################################################################################

def count_series(list_wanted, INFO_TAB):

    """
    Create a list with the count of each phylum wanted

    :param list_wanted: list of all the phylum wanted
    :type: list of str
    :param INFO_TAB: the name of a tabulate file with the information about the replicons
    :type: str
    :return: a Series that contain the count of each phylum wanted
    :rtype: pd.Series
    """

    info_tab = pd.read_table(INFO_TAB, index_col=0, names=["Taxon_id", "Name", "Kingdom", "Phylum", "Lineage", "NC_ids"],skiprows=1)

    series_wanted = pd.Series(0, index=list_wanted)

    for wanted in list_wanted :
        series_wanted.loc[wanted] = info_tab[info_tab.Lineage.str.contains(wanted)].shape[0]


    return series_wanted


##########################################################################################
##########################################################################################

def dataframe_color(PATH_TO_HTLM, df, dict_wanted, list_wanted, INFO_TAB, w_file):

	"""
	Function that plot the proportion of each systems found

	:param PATH_TO_HTLM: The absolute path to the figure folder
	:type: str
	:param df: The dataframe with the counts with the detail of the phylogeny
	:type: pandas.Dataframe
	:param dict_wanted: Dictionary that contain in key the kingdoms name and in value the list
	of all the phylum for all the phylum wanted
	:type: dict
	:param list_wanted: list of all the phylum wanted
	:type: list of str
	:param INFO_TAB: Name of the annotation table with the information about the database
	:type: str
	:param w_file: Open file where the textual information will be write
 	:type: file
	:return: Nothing
	"""

	series_wanted = count_series(list_wanted, INFO_TAB)
	count_df=pd.DataFrame(series_wanted.values, index=pd.MultiIndex.from_tuples([(kingdom,phylum) for kingdom in ['Bacteria', 'Archaea'] for phylum in dict_wanted[kingdom] ]), columns=['Count'])

	# XXX On met le tableau de compte dans le fichier texte ce qui peut être plus rapide pour trouver combien j'ai de génomes des espèces voulut dans ma database
	w_file.write("\n\nNumber of genomes wanted in the database used\n")
	w_file.write("---------------------------------------------\n\n")
	# XXX Je transpose pas car sinon c'est trop grand en ligne, plus lisible en colonne
	count_df.to_string(w_file)


	# XXX On met un tableau recapitulatif en texte
	w_file.write("\n\nPercentage Proteobacteria (𝛼, 𝛽, 𝛾) versus other\n")
	w_file.write("------------------------------------------------\n\n")
	mini_tab = pd.DataFrame(0, index=['Proteobacteria (𝛼, 𝛽, 𝛾)','Rest'] ,columns=["Count"])
	mini_tab.loc['Proteobacteria (𝛼, 𝛽, 𝛾)'] = count_df.loc['Bacteria', 'Gammaproteobacteria'] + count_df.loc['Bacteria', 'Betaproteobacteria'] + count_df.loc['Bacteria', 'Alphaproteobacteria']
	mini_tab.loc['Rest'] = count_df.sum() - mini_tab.loc['Proteobacteria (𝛼, 𝛽, 𝛾)']
	mini_tab.T.to_string(w_file)

	count_df.index.names = ["Kingdom", "Phylum"]

	drop_df = df.drop([('Bacteria', 'Other'), ('Archaea', 'Other'), "Summary_total", ('Bacteria', 'Total_system'), ('Archaea', 'Total_system')])

	maxi = max([len(n) for n in drop_df.columns])

	drop_df_rename = drop_df.rename(columns=lambda x: "{:_^{number}}".format(x, number=maxi))
	style_df = drop_df_rename.style.apply(do_gradient, args=(count_df,))

	# XXX Ecriture de la dataframe en HTML
	my_file = open(os.path.join(PATH_TO_HTLM,'dataframe_color.html'), 'w')
	my_file.write(style_df._repr_html_())
	my_file.close()

	# XXX La même en pdf mais ecriture mal lu par illustrator // Maintenant semble fonctionner, matplotlib est passé à sa version 2.0
	html_df=HTML(string=style_df._repr_html_())
	html_df.write_pdf(os.path.join(PATH_TO_HTLM,"dataframe_color.pdf"))

	return
