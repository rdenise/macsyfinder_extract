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
##								Functions
##
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

def count_all(info_tab, DICT_INFO, protein_function, LIST_SYSTEMS, PATH_TO_DATAFRAME, speciesDict) :

	"""
	Function that will read the DICT_INFO and count all information about system, protein, species...
	It will write a file with this information and return a dataframe.
	It's a count by protein so I'll hav the number of ATPase in T2SS for Proteobacteria

	:param info_tab: the genfromtxt of a tabulate file with the information about the replicons
	:type: numpay.ndarray
	:param DICT_INFO: dictionnary that contains all the information about all the
	systems found (create by rename_name_gene)
	:type: str
	:param protein_function: dictionnary with the name of all the protein studied
	:type: dict
	:params LIST_SYSTEMS: name of all the systems of the analysis
	:type: list of str
	:param PATH_TO_DATAFRAME: absolute path to the results folder
	:type: str
	:param speciesDict: Dictionary that contain in key the kingdoms name and in value the list
	of all the phylum
	:type: dict
	:return: a dataframe with all this count of the systems
	:rtype: pandas.Dataframe
	"""

	proteinlist = set(protein_function.values())
	number_of_line = (len(speciesDict["Bacteria"])+len(speciesDict["Archaea"]))*len(proteinlist)*2+(len(speciesDict["Bacteria"])+len(speciesDict["Archaea"]))+3

	info_tab = {line[0]:line[1:] for line in info_tab}

	miindex = tuple_like_all(speciesDict, proteinlist)
	number_systems = len(LIST_SYSTEMS)+1
	df_count_system = pd.DataFrame(np.zeros((number_of_line,number_systems), dtype=int), index=miindex, columns=LIST_SYSTEMS+['Total'])

	list_seen=[]

	for species_id in DICT_INFO :
		for key_systems in DICT_INFO[species_id] :
			system_name, number = key_systems.split("_")
			kingdom = info_tab[species_id][1]
			phylum = info_tab[species_id][2]

			if phylum in ['Archaea','Bacteria'] :
				phylum = 'Other'

			df_count_system.loc[(kingdom,phylum,"Total_system",""),system_name] += 1
			df_count_system.loc[(kingdom,"Total_system","",""),system_name] += 1
			df_count_system.loc[("Summary_total","","",""),system_name] += 1

			for protein_name in DICT_INFO[species_id][key_systems] :
				df_count_system.loc[(kingdom,phylum,protein_function[protein_name],'Unique'),system_name] += 1
				df_count_system.loc[(kingdom,phylum,protein_function[protein_name],'Total'),system_name] += DICT_INFO[species_id][key_systems][protein_name]

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

def systems_count(info_tab, DICT_INFO, protein_function, LIST_SYSTEMS, PATH_TO_DATAFRAME, list_wanted, speciesDict) :

	"""
	Function that will read the DICT_INFO and count all information about system, protein, species...
	It will write a file with this information and return a dataframe.
	It's a count by systems so I'll have the number of T2SS for Proteobacteria

	:param info_tab: the genfromtxt of a tabulate file with the information about the replicons
	:type: numpay.ndarray
	:param DICT_INFO: dictionnary that contains all the information about all the
	systems found (create by rename_name_gene)
	:type: str
	:param protein_function: dictionnary with the name of all the protein studied
	:type: dict
	:params LIST_SYSTEMS: name of all the systems of the analysis
	:type: list of str
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

	number_of_line = len(list_wanted)+5   # IDEA Ici si je veux le rendre général le 5 c'est "nombre de kingdom*2+1"
	proteinlist = protein_function.keys()

	# NOTE Pas besoin de ces lignes ici car je n'ai pas de verifié normalement
	#LIST_SYSTEMS.remove('T4bP')
	#LIST_SYSTEMS.remove('Com')

	info_tab = {line[0]:line[1:] for line in info_tab}

	miindex = tuple_like_systems(speciesDict)
	number_systems = len(LIST_SYSTEMS)
	df_count_system = pd.DataFrame(np.zeros((number_of_line, number_systems), dtype=int), index=miindex, columns=LIST_SYSTEMS)

	list_seen=[]

	for species_id in DICT_INFO :
		for key_systems in DICT_INFO[species_id] :
			system_name, number = key_systems.split("_")
			kingdom = info_tab[species_id][1]

			phylum = info_tab[species_id][2]
			if phylum not in list_wanted :
				phylum = "Other"

			df_count_system.loc[(kingdom,phylum),system_name] += 1
			df_count_system.loc[(kingdom,"Total_system"),system_name] += 1
			df_count_system.loc[("Summary_total",""),system_name] += 1


	df_count_system.index.names = ["Kingdom", "Phylum"]
	df_count_system.to_csv(os.path.join(PATH_TO_DATAFRAME,"tsv","data_frame_count_systems.tsv"), sep="\t", float_format="%.1f", index_label=False, encoding="utf-8" )
	df_count_system.to_excel(os.path.join(PATH_TO_DATAFRAME,"xlsx","data_frame_count_systems.xlsx"), index_label=False)
	return df_count_system

##########################################################################################
##########################################################################################

def proportion_phylum(PATH_TO_FIGURE, df, LIST_SYSTEMS):

	"""
	Function that plot the proportion of each phylum in all the systems found

	:param PATH_TO_FIGURE: The absolute path to the figure folder
	:type: str
	:param df: The dataframe with the counts with the detail of the phylogeny
	:type: pandas.Dataframe
	:param LIST_SYSTEMS: list with the name of all the species
 	:type: list of str
	:return: Nothing
	"""

	sub_df_for_plot=df.xs("Total_system", level="Function")
	name_legend = np.array([i[1] for i in sub_df_for_plot.index])

	# XXX Plot en lui même
	fig = plt.figure()
	plot = (sub_df_for_plot.loc[:,LIST_SYSTEMS]/df.xs("Summary_total", level="Kingdom").loc[:,LIST_SYSTEMS].ix[0]).T.plot(kind="bar", stacked=True, cmap='Paired', rot=0)

	# XXX Mise en forme
	plt.legend(name_legend,loc='center left', bbox_to_anchor=(1.0, 0.5))
	plt.subplots_adjust(right=0.75)
	plt.ylabel("Proportion")
	plt.xlabel("System")
	plt.title("Proportion of Phylum in the studied system")
	ax_phylum = plt.gca()

	# XXX Sauvegarde du dictionnaire de couleurs pour une future utilisation qui peut être utilisé plus tard pour les couleurs
	info_legend = ax_phylum.get_legend_handles_labels()
	taille_infolegend =len(info_legend[1])
	colors_dict = {info_legend[1][i]:info_legend[0][i][0].get_facecolor() for i in range(taille_infolegend)}
	f = open(os.path.join(PATH_TO_FIGURE, ".color_dict.json"), 'w')
	json.dump(colors_dict, f)

	# XXX Sauvegarde de la figure
	plt.savefig(os.path.join(PATH_TO_FIGURE,"proportion_phylum.pdf"))

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
	w_file.write("Percentage dans les kingdoms\n")
	w_file.write("----------------------------\n\n")
	mini_tab = df.xs("Total_system", level="Phylum").iloc[:,:5].div(df.xs("Total_system", level="Phylum").Total, axis='index')*100
	mini_tab.to_string(w_file)

	# XXX Et le même en figure
	fig = plt.figure()
	ax = mini_tab.plot(kind='bar', rot=0)

	# XXX Mise en forme
	plt.xlabel("System")
	plt.ylabel("Percentage")
	ax.set_xticklabels(["Archaea", "Bacteria"])
	plt.title("Percentage of studied system found in each kingdom clades studied")
	plt.savefig(os.path.join(PATH_TO_FIGURE,"proportion_found_each_kingdom.pdf"))

	# XXX Proportion générale
	fig = plt.figure()
	ax = (df.xs("Summary_total", level="Kingdom").iloc[:,:5]/df.xs("Summary_total", level="Kingdom").Total.ix[0]).plot(kind='bar')

	# XXX Mise en forme
	plt.xlabel("System")
	plt.ylabel("Proportion")
	ax.set_xticklabels("")
	plt.title("Proportion of studied system found")
	plt.savefig(os.path.join(PATH_TO_FIGURE,"proportion_found.pdf"))

	return

##########################################################################################
##########################################################################################

def proportion_proteobacteria(PATH_TO_FIGURE, df):

	"""
	Function that plot the proportion of each systems found

	:param PATH_TO_FIGURE: The absolute path to the figure folder
	:type: str
	:param df: The dataframe with the counts with the detail of the phylogeny
	:type: pandas.Dataframe
	:return: Nothing
	"""

	# XXX Creation du dataframe intermediare
	df_figure_2 = pd.DataFrame(index=pd.MultiIndex.from_tuples([('Bacteria', 'Proteobacteria'),('Rest', '')]) ,columns=df.columns)
	df_figure_2.loc['Bacteria', 'Proteobacteria'] = df.loc['Bacteria', 'Gammaproteobacteria'] + df.loc['Bacteria', 'Betaproteobacteria'] + df.loc['Bacteria', 'Alphaproteobacteria']
	df_figure_2.loc['Rest', ''] = df.loc['Summary_total', ''] - df_figure_2.loc['Bacteria', 'Proteobacteria']
	df_figure_2.sum(axis=1)

	# XXX Le plot avec les couleurs choisies
	fig = plt.figure()
	colors = [(0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0), (0.3997693305214246, 0.6478123867044262, 0.80273742044673246, 1.0)]
	plot = df_figure_2.T.plot(kind="bar", stacked=True, color=colors, rot=0, legend=False)
	plt.savefig(os.path.join(PATH_TO_FIGURE,"numbers_proteobacteria_rest.pdf"))

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

def count_series(genome_list, list_wanted, info_tab):
	"""
	Create a list with the count of each phylum wanted

	:param genome_list: name of a file that contains a list of all the genome in the initial database
	:type: str
	:param list_wanted: list of all the phylum wanted
	:type: list of str
	:param info_tab: the genfromtxt of a tabulate file with the information about the replicons
	:type: numpay.ndarray
	:return: a Series that contain the count of each phylum wanted
	:rtype: pd.Series
	"""
	max_columns=[]
	with open(genome_list, 'r') as r_file:
	     for line in r_file :
		        max_columns.append(len(line.split()))
	max_columns = max(max_columns)

	dict_position = {phylum:list_wanted.index(phylum) for phylum in list_wanted}
	info_tab = {line[0]:line[1:] for line in info_tab}

	len_wanted = len(list_wanted)
	series_wanted = pd.Series(data=np.zeros(len_wanted), index=list_wanted)

	genome_list = pd.read_csv(genome_list, sep=' ', names=[i for i in range(max_columns)], header=None, engine='python')
	for species in genome_list.iloc[:,1] :
		if info_tab[species][2] in list_wanted :
			series_wanted[info_tab[species][2]] += 1

	return series_wanted


##########################################################################################
##########################################################################################

def dataframe_color(PATH_TO_HTLM, df, dict_wanted, list_wanted, INFO_STATS, genome_list):

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
	:param info_tab: the genfromtxt of a tabulate file with the information about the replicons
	:type: numpay.ndarray
	:param genome_list: name of a file that contains a list of all the genome in the initial database
	:type: str
	:return: Nothing
	"""

	count_df=pd.DataFrame(count_series(genome_list, list_wanted, INFO_STATS).values, index=pd.MultiIndex.from_tuples([(kingdom,phylum) for kingdom in ['Bacteria', 'Archaea'] for phylum in dict_wanted[kingdom] ]), columns=['Count'])
	count_df.index.names = ["Kingdom", "Phylum"]

	drop_df = df.drop([('Bacteria', 'Other'), ('Archaea', 'Other'), "Summary_total", ('Bacteria', 'Total_system'), ('Archaea', 'Total_system')])

	maxi = max([len(n) for n in drop_df.columns])

	drop_df_rename = drop_df.rename(columns=lambda x: x+"_"*(maxi-len(x)))
	style_df = drop_df_rename.style.apply(do_gradient, args=(count_df,))

	# XXX Ecriture de la dataframe en HTML
	my_file = open(os.path.join(PATH_TO_HTLM,'dataframe_color.html'), 'w')
	my_file.write(style_df._repr_html_())
	my_file.close()

	# XXX La même en pdf mais ecriture mal lu par illustrator
	html_df=HTML(string=style_df._repr_html_())
	html_df.write_pdf(os.path.join(PATH_TO_HTLM,"dataframe_color.pdf"))

	return
