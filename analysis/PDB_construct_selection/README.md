Selection of PDB constructs suitable for expression
===================================================

Description
-----------

These directories contain data generated using the script select-PDB-constructs.py and related scripts.

The idea here is to select sequences likely to express well, by ranking PDB construct sequences based on parameters such as a given expression host, alignment score (against the wild-type sequence), and the likely authenticity of the construct sequence (since PDB files are frequently misannotated).

Manifest
--------

* PDB\_construct\_selections.txt:
    * text file output by the select-PDB-constructs.py scripts; displays various results for each target protein, including the top-ranked PDB chain.
* alignments/:
    * HTML files output by the select-PDB-constructs.py scripts; for each target protein which has PDB entries with the desired expression\_system annotation, an HTML file displays the sequence alignment of those PDB entries against the UniProt canonical isoform sequence, and a sequence from a plasmid library; the order of the PDB sequences is sorted based on a number of parameters such as those above.
* notes-first-pass/:
    * notes written up after looking at the references in detail for some example target proteins/PDBs
* notes-expr-tags/:
    * notes written up after running the select-PDB-constructs-expr-tags.py script

