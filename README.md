# coca-thesis
This repository contains the materials used during the writing of the "Analysing Dependency Structure of Coordination Using a Constituency Parser and Dependency Length MinimizationPython" bachelor thesis. It contains Python code used to extract the coordination data and the R code used for the statistical analysis.

- coord_gov_BP.py - finds the governor of the coordination. Used by find_coord_data.py
- extract_official.py - use to extract the lengths of conjuncts and produce the final csv table.
- find_coord_data.py - use to find coordination data from a file_raw.csv.
- parseTSV.py - use to parse the TSV file (file.tsv) with the standard BNP model. Results saved to a new csv file (/csv/parse_raw/file_raw.csv).
- parseTSV_large.py - use to parse the TSV file (file.tsv) with the large BNP model. Results saved to a new csv file (/csv/parse_raw_large/file_raw.csv).
- syllables_BP.py - determines syllables length. Used by extract_official.py
- eval.xlsx - evaluation data used to determine accuracy. Not necesary for the Python code 

To install the packages run: pip install -r requirements.txt
