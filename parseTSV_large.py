#!/usr/bin/env python
"""
parseTSV.py
Load the tsv file(s) passed as command-line argument(s). Parse each
sentence in the file using the Berkeley parser. Write each sentence
parse tree to a csv file. If the input is named 'sample.txt' (or
'/data/sample.txt', or the input is in any number of nested
subdirectories), the result is written to '/csv/parse_raw/sample_raw.csv'.
This code has been heavily influenced by, and parts of it adapted from Julie Kallini's senior thesis project
source: https://github.com/jkallini/PrincetonThesis

"""
import csv
from sys import stderr
import spacy

import pandas as pd

import argparse
import os

from tqdm import tqdm


def get_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Generate parse tree for each line of the cleaned input file(s).')
    parser.add_argument('input_files', nargs='+', type=str,
                        help='path to input file(s)')
    return parser.parse_args()


# Main function.
if __name__ == "__main__":

    args = get_args()

    # Load spacy model
    print("Loading spaCy's large English model...")
    nlp = spacy.load("en_core_web_lg")

    # Disable sentence segmentation
    # - this seems not to work. That is why in find_coord_data.py the sentences are filtered again
    nlp.select_pipes(disable="senter")

    # Integrate with benepar
    print("Integrating spaCy model with Benepar...")

    # Choose the benepar model used for parsing: benepar_en3 or benepar_en3_large
    nlp.add_pipe('benepar', config={'model': 'benepar_en3_large'})

    print()

    i = 1
    tot = str(len(args.input_files))

    input_files = reversed(args.input_files)

    for file in input_files:

        print("(" + str(i) + "/" + tot + ")")
        print("Beginning parse of " + file
              + "! If the input file is large, this may take a few hours...")
        sent_list = pd.read_csv(file, delimiter=r'\s*\t\s*', quoting=csv.QUOTE_NONE, quotechar="", lineterminator="\n")
        data = []

        # Parsing the sentences
        for index, row in tqdm(sent_list.iterrows(), total=len(sent_list.index)):
            coca_sen = row["SENT"]
            coca_sen_id = row["0"]
            coca_type = row["TYPE"]
            coca_year = row["YEAR"]
            coca_file = row["FILE"]
            coca_para_id = row["PARA ID"]
            try:
                doc = nlp(coca_sen)
                for sent in doc.sents:
                    data.append([coca_sen_id,
                                 coca_type,
                                 coca_year,
                                 coca_file,
                                 coca_para_id,
                                 coca_sen,
                                 # Access the parse tree
                                 sent._.parse_string])

            except Exception as e:
                print(str(e), file=stderr)

        columns = ['Sentence Id',
                   'COCA Type',
                   'COCA Year',
                   'COCA File',
                   'COCA Paragraph Id',
                   'COCA Sentence Text',
                   'Sentence Parse Tree']
        df = pd.DataFrame(data, columns=columns)

        # Saving the results to a directory /csv/parse_raw
        dest_name = os.path.splitext(os.path.basename(file))[-2]
        dest_dir = 'csv/parse_raw_large'

        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        # Adding "_raw" suffix to the file and saving it to CSV
        dest = dest_dir + '/' + dest_name + '_raw.csv'
        df.to_csv(dest, index=False)
        print("All done! The result is stored in " +
              dest_dir + '/' + dest_name + '.csv.\n')
        i = i + 1
