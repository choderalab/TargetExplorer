Selection of PDB constructs suitable for expression
===================================================

Description
-----------

These directories contain data generated using the script select-PDB-constructs.py.

The aim is to select construct sequences which are likely to express well, by ranking PDB construct sequences based on parameters such as the expression host, alignment score (against the wild-type sequence), and the likely authenticity of the construct sequence (since PDB files are frequently misannotated).

Manifest
--------

* ../scripts/select-PDB-constructs.py
    * selects PDB constructs based on various criteria and outputs data useful for initiating expression tests.
* PDB\_constructs-data.txt:
    * text file output by the select-PDB-constructs.py script; displays various results for each target protein, including the top-ranked PDB chain and the target\_score from the database.
* PDB\_constructs.xlsx:
    * spreadsheet containing DNA sequences of chosen constructs, along with various related data, to be used to initiate expression testing.
* manual\_exceptions.yaml:
    * list of manual exceptions to be made for construct choices - this is read in by the select-PDB-constructs.py scripts and used to 
* alignments/:
    * HTML files output by the select-PDB-constructs.py scripts; for each target protein domain which has PDB entries with the desired expression\_system annotation, an HTML file displays the sequence alignment of those PDB entries against the UniProt canonical isoform sequence, and a sequence from a plasmid library; the order of the PDB sequences is sorted based on a number of parameters such as those above.

