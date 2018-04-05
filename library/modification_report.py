# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Library
##
##########################################################################################
##########################################################################################

import os
import pandas as pd
import numpy as np

##########################################################################################
##########################################################################################
##
##								Functions
##
##########################################################################################
##########################################################################################

def read_systems_found(systems_found_file):

    """
    Read .names of macsyfinder_extract and created a dataframe with the rigth type on each columns

    :param systems_found_file: Name of the macsyfinder_extract systems_found.names file
    :type: str
    :return: the dataframe completed
    :rtype: pandas.Dataframe
    """

    names_dataframe=["Species_Id","Replicon_name","System_Id", "System_name", "System_status","System_number","Proteins","Kingdom","Phylum","Lineage", "Species"]
    summary_df = pd.read_table(systems_found_file, names=names_dataframe, comment="#")
    summary_df["Proteins"] = summary_df.apply(lambda x : eval(x.Proteins), axis=1)
    summary_df["System_number"] = summary_df.apply(lambda x : str(int(x.System_number)), axis=1)

    return summary_df

##########################################################################################
##########################################################################################


def dict_tandem(report_df) :

    """
    Function that create a series with the information about if the genes is in tandem in the cluster

    :param report_df: The dataframe with the information about the gene, the system_id and the position of the gene in the chromosome
    :type: pandas.DataFrame
    :return: a pandas.Series with Yes or No if the gene is in tandem
    :rtype: pandas.Series
    """

    big_dict = {}
    for group in report_df.groupby("System_Id") :
        tandem = []
        sub_df = group[1]
        for index, line in sub_df.iterrows() :
            if tandem == [] :
                tandem.append([index, (line.Gene, line.Position)])
            elif (tandem[-1][1][0] == line.Gene) and (line.Position - tandem[-1][1][1] == 1):
                tandem.append([index, (line.Gene, line.Position)])
            elif len(tandem) > 1:
                big_dict.update({myindex:"Yes" for myindex, other in tandem})
                tandem=[[index, (line.Gene, line.Position)]]
            else :
                big_dict[tandem[0][0]] = "No"
                tandem=[[index, (line.Gene, line.Position)]]
        if len(tandem) > 1:
            big_dict.update({myindex:"Yes" for myindex, other in tandem})
        else :
            big_dict[tandem[0][0]] = "No"
    return pd.Series(big_dict)

##########################################################################################
##########################################################################################

def dict_loner(report_df, distance_max):

    """
    Function that create a series with the information about if the genes is loner in the cluster => ['pilD'] = loner ; ['pilD','pilD'] != loner

    :param report_df: The dataframe with the information about the gene, the system_id and the position of the gene in the chromosome
    :type: pandas.DataFrame
    :param distance_max: A dictionnary with the maximum distance between the gene for each systems
    :type: dict
    :return: a pandas.Series with Yes or No if the gene is in loner
    :rtype: pandas.Series
    """

    big_dict = {}
    for group in report_df.groupby("System_Id") :
        cluster = []
        sub_df = group[1]
        for index, line in sub_df.iterrows() :
            if cluster == [] :
                cluster.append([index, line.Position])
            elif line.Position - cluster[-1][1] <= distance_max[line.Predicted_system]:
                cluster.append([index, line.Position])
            elif len(cluster) > 1:
                big_dict.update({myindex:"No" for myindex, posiotion in cluster})
                cluster=[[index, line.Position]]
            else :
                big_dict[cluster[0][0]] = "Yes"
                cluster=[[index, line.Position]]
        if len(cluster) > 1:
            big_dict.update({myindex:"No" for myindex, position in cluster})
        else :
            big_dict[cluster[0][0]] = "Yes"
        cluster=[index]
    return pd.Series(big_dict)

##########################################################################################
##########################################################################################

def dict_loner_uniq(report_df, distance_max) :

    """
    Function that create a series with the information about if the genes is loner in the cluster => ['pilD'] = loner ; ['pilD','pilD'] = loner

    :param report_df: The dataframe with the information about the gene, the system_id and the position of the gene in the chromosome
    :type: pandas.DataFrame
    :param distance_max: A dictionnary with the maximum distance between the gene for each systems
    :type: dict
    :return: a pandas.Series with Yes or No if the gene is in loner
    :rtype: pandas.Series
    """

    big_dict = {}
    for group in report_df.groupby("System_Id") :
        cluster = []
        sub_df = group[1]
        for index, line in sub_df.iterrows() :
            if cluster == [] :
                cluster.append([index, (line.Position, line.Gene)])
            elif line.Position - cluster[-1][1][0] <= distance_max[line.Predicted_system]:
                cluster.append([index, (line.Position, line.Gene)])
            elif len(cluster) > 1:
                big_dict.update({myindex:"No" for myindex, other in cluster})
                cluster=[[index, (line.Position, line.Gene)]]
            else :
                big_dict[cluster[0][0]] = "Yes"
                cluster=[[index, (line.Position, line.Gene)]]
        if len(cluster) > 1:
            if np.unique([p[1] for i, p in cluster]).shape[0] == 1 :
                big_dict.update({myindex:"Yes" for myindex,other in cluster})
            else :
                big_dict.update({myindex:"No" for myindex, other in cluster})
        else :
            big_dict[cluster[0][0]] = "Yes"
    return pd.Series(big_dict)

##########################################################################################
##########################################################################################

def fill_reportdf_detected(report_df, protein_wanted, report_df_modif, info_tab, folder, distance_max, names_dataframe) :

    """
    Function that set new columns in the report dataframe : NewName, MultiCopy, Status of the system, the kingdom and the phylum

    :param report_df: The report dataframe from the merge report file from the macsyfinder analysis
    :type: pandas.Dataframe
    :param protein_wanted: the list of all the gene name wanted
    :type: list of str
    :param report_df_modif: the dataframe modified by the extract_protein function
    :type: pandas.Dataframe
    :param info_tab: The dataframe with information about the kingdom and the phylum of each species
    :type: pandas.DataFrame
    :param folder: The name of the folder where the information will be write
    :type: str
    :param distance_max: A dictionnary with the maximum distance between the gene for each systems
    :type: dict
    :param names_dataframe: the list with the name of the columns of the report dataframe
    :type: list of str
    :return: a dataframe with the columns wanted added
    :rtype: pandas.DataFrame
    """

    print("-----------------------")
    print("|Create report detected")
    print("-----------------------")

    # adding information about the cluster
    report_df['Tandem'] = dict_tandem(report_df)
    report_df["Loner"] = dict_loner(report_df, distance_max)
    report_df["Loner_unique"] = dict_loner_uniq(report_df, distance_max)

    sub_report = report_df[report_df.Gene.isin(protein_wanted)].reset_index(drop=True)

    #adding new name (plus test if the two dataframes are in the same order)
    if sum(sub_report.Hit_Id == report_df_modif.Hit_Id) == sub_report.shape[0] :
        sub_report['NewName'] = report_df_modif['NewName']
    else :
        sub_report = sub_report.sort_values(['System_Id','Position']).reset_index(drop=True)
        report_df_modif = report_df_modif.sort_values(['System_Id','Position']).reset_index(drop=True)
        sub_report['NewName'] = report_df_modif['NewName']

    #addind if the gene have multicopy
    new_sub_report = sub_report.set_index(["System_Id", "Gene"])
    new_sub_report["Multi_copy"] = sub_report.groupby("System_Id").apply(lambda x : x.groupby("Gene").apply(lambda y: "Yes" if y.shape[0]>1 else "No"))

    #reorder and reset index for sub_report
    sub_report = new_sub_report.reset_index()[["NewName"]+names_dataframe+["Tandem","Loner", "Loner_unique","Multi_copy"]].sort_values(['System_Id', 'Position']).reset_index(drop=True)

    #set the sequence with the max score
    sub_report["Max_Score"] = "No"
    sub_report.loc[sub_report.groupby(["System_Id", "Gene"])["Score"].idxmax().values,"Max_Score"] = "Yes"
    sub_report["Min_Score"] = "No"
    sub_report.loc[sub_report.groupby(["System_Id", "Gene"])["Score"].idxmin().values,"Min_Score"] = "Yes"

    #set the sequence with the min ievalue
    sub_report["Min_ievalue"] = "No"
    sub_report.loc[sub_report.groupby(["System_Id", "Gene"])["i-evalue"].idxmin().values,"Min_ievalue"] = "Yes"
    sub_report["Max_ievalue"] = "No"
    sub_report.loc[sub_report.groupby(["System_Id", "Gene"])["i-evalue"].idxmax().values,"Max_ievalue"] = "Yes"

    if sub_report.loc[sub_report.Max_Score != sub_report.Min_ievalue].shape[0] == 2 :
        sub_report.loc[sub_report.Max_Score != sub_report.Min_ievalue, "Min_ievalue"] = "Yes"

    # adding information about the status of the system
    sub_report['System_status'] = "D"

    #adding phylum and kingdom
    infoloc = info_tab.loc
    sub_report['Kingdom'] = sub_report.apply(lambda x : infoloc[x.Replicon_name[:-5], "kingdom"], axis=1)
    sub_report['Phylum'] = sub_report.apply(lambda x : infoloc[x.Replicon_name[:-5], "phylum"], axis=1)
    sub_report['Lineage'] = sub_report.apply(lambda x : infoloc[x.Replicon_name[:-5], "lineage"], axis=1)
    sub_report['Species_name'] = sub_report.apply(lambda x : infoloc[x.Replicon_name[:-5], "name"], axis=1)

    #sub_report.to_excel(os.path.join(folder, "report_modif_loner_tandem_detected.xlsx"),index=False)
    sub_report.to_csv(os.path.join(folder, "detected.report"),index=False, sep="\t")

    print()
    print("Done!")
    print()

    return sub_report

##########################################################################################
##########################################################################################

def fill_reportdf_validated(df_report_validated, df_found, folder):

    """
    Function that set new columns in the report dataframe : Tandem, Loner, Loner unique,
    NewName, MultiCopy, Max score, Min evalue, Status of the system, the kingdom and the phylum

    :param df_report_validated: the dataframe modified by the create_validated_fasta function
    :type: pandas.Dataframe
    :param df_found: The dataframe create by the function set_df_info_system
    :type: pandas.DataFrame
    :param folder: The name of the folder where the information will be write
    :type: str
    :return: a dataframe with the columns wanted added with the validated sequences
    :rtype: pandas.DataFrame
    """

    print("-----------------------")
    print("|Create report validated")
    print("-----------------------")

    # On met le nom du replicon DEJA
    #df_report_validated["Replicon_name"] = df_report_validated.apply(lambda x : x.System_Id.split("_")[0], axis=1)

    # On dit quel est le nom du système de chaques sequences dans la colonne et on renome T4bP les systemes qui ne sont pas T2SS, T4P, Tad, Com, MSH et le reste T4bP
    #df_report_validated["Predicted_system"] = df_report_validated.apply(lambda x : x.System_Id.split("_")[1] if x.System_Id.split("_")[1] in ("T4P","T2SS", "Tad", "Com", "MSH") else "T4bP", axis=1)
    #df_report_validated["Reference_system"] = df_report_validated["Predicted_system"]

    # On dit que le status sont validés
    df_report_validated["System_status"] = "V"

    # On merge les informations de df_found dand le report pour avoir le nom du kingdom, phylum.
    df_report_validated = df_report_validated.merge(df_found[["Replicon_name", "Kingdom", "Phylum", "Lineage", "Species"]].drop_duplicates(), on="Replicon_name")

    #Place les multicopy
    df_report_validated['Multi_copy'] = "No"

    new_df_report_validated = df_report_validated.set_index(["System_Id", "Gene"])
    #a verifé
    #df_report_validated.to_csv(os.path.join(folder, "validated_look.report"),index=False, sep="\t")
    new_df_report_validated["Multi_copy"] = df_report_validated.groupby("System_Id").apply(lambda x : x.groupby("Gene").apply(lambda y: "Yes" if y.shape[0]>1 else "No"))

    #reorder and reset index for sub_report
    #a finir normalement c'est fini
    df_report_validated = new_df_report_validated.reset_index()[['NewName', 'Hit_Id', 'Replicon_name', "Sequence_length", 'Gene', 'Reference_system', 'Predicted_system', "System_Id", "System_status","Multi_copy", 'Kingdom', 'Phylum', "Lineage", "Species"]].sort_values(['System_Id']).reset_index(drop=True)
    df_report_validated.rename(columns={'Species': 'Species_name'}, inplace=True)

    #df_report_validated.to_excel(os.path.join(folder, "report_modif_validated.xlsx"),index=False)
    df_report_validated.to_csv(os.path.join(folder, "validated.report"),index=False, sep="\t")
    # XXX On remove ici le fichier qui est maintenant inutile
    os.remove(os.path.join(folder, "validated_tmp.report"))

    print()
    print("Done!")
    print()

    return df_report_validated

##########################################################################################
##########################################################################################


def merge_detected_validated_report(report_file_validated, report_file_detected, report_detected_modif, file_found, info_tab, protein_def, info_folder, DISTANCE_DICT) :

    """
    Function that set new columns in the report dataframe : Tandem, Loner, Loner unique,
    NewName, MultiCopy, Max score, Min evalue, Status of the system, the kingdom and the phylum

    :param report_df_validated: the name dataframe modified by the create_validated_fasta function
    :type: str
    :param report_file_detected: The name report dataframe from the merge report file from the macsyfinder analysis
    :type: str
    :param report_detected_modif: the name of the file the dataframe modified by the extract_protein function
    :type: str
    :param file_found: The name of dataframe create by the function set_df_info_system
    :type: str
    :param info_tab: The dataframe with information about the kingdom and the phylum of each species
    :type: pandas.DataFrame
    :param protein_def: the dictionary of the names of the gene wanted in index and the function in value
    :type: dict
    :param info_folder: The name of the folder where the information will be write
    :type: str
    :param DISTANCE_DICT: the dictionnary with the maximal distance inside systeme
    :type: dict
    :return: a dataframe with the two reports modified merge
    :rtype: pandas.DataFrame
    """

    # On traite d'abord le report des detectés
    names_dataframe=['Hit_Id','Replicon_name','Position','Sequence_length','Gene','Reference_system','Predicted_system','System_Id','System_status','Gene_status','i-evalue','Score','Profile_coverage','Sequence_coverage','Begin_match','End_match']
    report_df = pd.read_table(report_file_detected, names=names_dataframe)

    # on charge aussi le report modifié des detecté
    report_df_modif = pd.read_table(report_detected_modif)

    # On charge le dictionnaire des protéines que je veux
    #dict_protein = set_dict_protein(protein_def)
    #protein_wanted = [item for sublist in list(dict_protein.values()) for item in sublist]

    # On donne les distances maximums pour chaque système
    distance_max = DISTANCE_DICT

    report_df_detected = fill_reportdf_detected(report_df, protein_def, report_df_modif, info_tab, info_folder, distance_max,names_dataframe)

    #Maintenant les verfiés
    report_df_validated = pd.read_table(report_file_validated)
    df_found = read_systems_found(file_found)
    df_found = df_found[df_found.System_status == "V"].reset_index(drop=True)

    report_df_validated = fill_reportdf_validated(report_df_validated, df_found, info_folder)

    print("----------------------")
    print("|Create report merged")
    print("----------------------")

    # On merge les deux report et on remet les colonnes dans le bon ordre
    report_df_merge = report_df_validated.append(report_df_detected, ignore_index=True)

    report_df_merge = report_df_merge[report_df_detected.columns]

    #report_df_merge.to_excel(os.path.join(info_folder,"report_modif_annotation_full.xlsx"),index=False, na_rep=".")
    report_df_merge.to_csv(os.path.join(info_folder,"report_modif_annotation_full.report"),index=False, sep="\t",na_rep=".")

    print()
    print("Done!")
    print()

    return report_df_merge
