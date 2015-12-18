# knoweng_nih
KNOWENG Project.
Author: Edward Huang

>>> python top_pathways_lincs.py 3
or
>>> python top_pathways_lincs.py 4
Reads in the LINCS data for level 3 or 4, and then finds the top pathway, drug
and cell-line pairs by using Fisher's hypergeometric test. The top genes for
each drug are found by taking z-scores the absolute values of which are higher
than 2.5.