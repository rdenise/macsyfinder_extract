# Macsyfinder Sequence Extractor

Script to extract information of macsyfinder .report

Dependencies :
--------------

- Python 3.5 (The Python 2.7 version is no longer updated)
   - Numpy 1.12
   - Biopython 1.68
   - Matplotlib 2.0.0
   - Pandas 0.19.2

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
This file have two columns, one with the kingdom and one phylum.

kingdom **[TAB]** phylum

Comments in the annotation table begin with **#**  

Format file distance systems
--------------------------
This file have two columns, one with the system and one maximal distance between two genes (as set in the .xml of the macsyfinder model).

system **[TAB]** distance

Comments in the annotation table begin with **#**  
