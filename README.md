# coca-thesis
This repository contains the materials used during the writing of the "Analysing Dependency Structure of Coordination Using a Constituency Parser and Dependency Length MinimizationPython" bachelor thesis. It contains Python code used to extract the coordination data and the R code used for the statistical analysis.

- bnp_models_download.py - use to easily download the BNP models.
- COCA_analysis.R - use for statistical analysis
- coord_gov_BP.py - finds the governor of the coordination. Used by find_coord_data.py
- extract_official.py - use to extract the lengths of conjuncts and produce the final csv table.
- find_coord_data.py - use to find coordination data from a file_raw.csv.
- find_coord_data_improved.py - improved version of find_coord_data.py. It fixes the issues pointed out in the limitation section of the thesis.
- parse_sentence.py - use to parse a single sentence with the standard BNP model and retrive its coordination data. Additionally, a syntax tree is drawn.
- parseTSV.py - use to parse the TSV file (file.tsv) with the standard BNP model. Results saved to a new csv file (/csv/parse_raw/file_raw.csv).
- parseTSV_large.py - use to parse the TSV file (file.tsv) with the large BNP model. Results saved to a new csv file (/csv/parse_raw_large/file_raw.csv).
- syllables_BP.py - determines syllables length. Used by extract_official.py
- eval.xlsx - evaluation data used to determine accuracy. Not necesary for the Python code 

The code was writen on Python version 3.9.13., and this is the recomended version.
To install the packages run: 
  ```
  pip install -r requirements.txt
  ```
To install the spacy large English model run:
  ```
  python -m spacy download en_core_web_lg
  ```
To install the benepar model run:
  ```
  bnp_models_download.py
  ```
