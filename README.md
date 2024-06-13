# coca-thesis
This repository contains the materials used during the writing of the "Analysing Dependency Structure of Coordination Using a Constituency Parser and Dependency Length MinimizationPython" bachelor thesis. It contains Python code used to extract the coordination data and the R code used for the statistical analysis.

## List of files ##
- **bnp_models_download.py** - use to easily download the BNP models.
- **COCA_analysis.R** - R code used for statistical analysis
- **coord_gov_BP.py** - finds the governor of the coordination. Used by find_coord_data.py
- **extract_official.py** - use to extract the lengths of conjuncts and produce the final csv table. Takes as input the results of find_coord_data.py or find_coord_data_improved.py.
- **find_coord_data.py** - use to find coordination data from a file_raw.csv, provided by parseTSV.py or parseTSV_large.py.
- **find_coord_data_improved.py** - improved version of find_coord_data.py. It fixes the issues pointed out in the limitation section of the thesis.
- **parse_sentence.py** - use to parse a single sentence with the standard BNP model and retrive its coordination data. Additionally, a syntax tree is drawn.
- **parseTSV.py** - use to parse the TSV file (file.tsv) with the standard BNP model. Results are saved to a new csv file in a specific folder path (/csv/parse_raw/file_raw.csv).
- **parseTSV_large.py** - use to parse the TSV file (file.tsv) with the large BNP model. Results are saved to a new csv file in a specific folder path (/csv/parse_raw_large/file_raw.csv).
- **syllables_BP.py** - determines syllables length. Used by extract_official.py
- **eval.xlsx** - evaluation data used to determine accuracy. Not necesary for the Python code 

## Step by step guide to run the BNP code ##
The code was writen on Python version 3.9.13., and this is the recomended version.
To use the code you need to:
1. Download Python version 3.9.13 from https://www.python.org/downloads/release/python-3913/ (other versions may or may not work).
2. In the installer, check "Add Python 3.9 to PATH" and make sure the "Install Now" option includes pip; if not, add it in "Customize installation".
3. After installation, verify that Python was installed correctly by running the following command in the command prompt (CMD):
   ```
   python --version
   ```
4. Download the code from GitHub at https://github.com/KattGaii/coca-thesis and save it in your folder of choice.
5. Navigate to the folder with the code in the command prompt:
   ```
   cd path_to_code
   ```
6. Activate the virtual environment. The terminal prompt should change to show (venv) before the folder path:
   ```
   python -m venv
   venv\Scripts\activate
   ```
7. Install the required packages:
   ```
   pip install -r requirements.txt
   ```
8. Download the spaCy language model:
   ```
   python -m spacy download en_core_web_lg
   ```
9. Download the BNP models:
   ```
   bnp_models_download.py
   ```
10. You can test the setup by parsing a single sentence using:
    ```
    parse_sentence.py This is an example sentence.
    ```
11. To parse a TSV file, use the parseTSV.py or parseTSV_large, depending on the model. It's recommended to navigate to the folder containing the TSV files and run the command from there:
    ```
    cd path_to_files
    path_to_code/parseTSV_large filename1 filename2 filename3 ...
    ```
12. To deactivate the virtual enviroment use:
   ```
   path_to_code\venv\Scripts\deactivate
   ```

