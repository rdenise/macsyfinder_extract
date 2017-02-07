# Macsyfinder Sequence Extractor

Script to extract information of macsyfinder .report

Dependencies :
--------------

- Python 3.5 (The Python 2.7 version is no longer updated)
   - Numpy 1.11.2
   - Biopython 1.68
   - Matplotlib 1.5.3
   - Rpy2 2.8.3

Format function definition file
-------------------------------
The function definition file need to be a tabulate separate file with these columns:  

Function name **[TAB]** system name_protein name 1 **[TAB]** system name_protein name 2 **...**  

The comments begin with **#**  

Format annotation file
----------------------
The annotation file need to be a tabulate separate file with these columns:  

species id **[TAB]** taxonomic NCBI id **[TAB]** full_name_of_the_species **[TAB]** kingdom **[TAB]** phylum  **[TAB]** lineage **[TAB]** NCBI ids of the replicon within the genome

Comments in the annotation table begin with **##**  

Format file wanted phylums
--------------------------
This file have two columns each line have one with the kingdom and one phylum.

kingdom **[TAB]** phylum

Comments in the annotation table begin with **#**  
