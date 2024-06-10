"""
extract_official.py
Takes working data provided by find_coord_data.py and calculates different lengths of conjuncts.

"""

import argparse
import ast
import os
import re
import nltk
from syllables_BP import count_word

import pandas as pd
from tqdm import tqdm
from nltk import Tree


def get_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Get coordination stats from csv input file(s) containing parsed sentences.')
    parser.add_argument('input_files', nargs='+', type=str,
                        help='path to input csv file(s)')
    return parser.parse_args()


def char_count(text):
    """
    Returns the length of string
    :param text: a string representing the conjunct for analysis
    :return: number of syllables or -9999 if error occurs.
    """
    if text:
        return len(text)
    return -9999


def syllable_count(text):
    """
    Returns the number of syllables in the string using the syllables.py
    Created by Michał Woźniak for the Przepiórkowski and Woźniak (2023)
    :param text: a string representing the conjunct for analysis
    :return: number of syllables or -9999 if error occurs.
    """
    if text:
        return count_word(text)
    return -9999


def token_count(conj_pos, parse_tree):
    """
    Returns the number of tokens
    :param conj_pos: position of the conjunct in the tree
    :param parse_tree: sentence parse tree
    :return: number of tokens or -9999 if error occurs.
    """
    if parse_tree and conj_pos:
        tree = Tree.fromstring(parse_tree)
        subtree = tree[conj_pos]
        return len(subtree.leaves())
    return -9999


def word_count(text):
    """
    Returns the number of words (no. of tokens) in the string
    :param text: a string representing the conjunct for analysis
    :return: number of words or -9999 if error occurs.
    """
    if text:
        return len(text.split())
    return -9999


def transform_to_latex_tree(input_str):
    """
    Transforms the nltk format tree into a format compatible with LaTeX TikZ trees.
    This enables using LaTeX or 3rd party websites to visualise the trees
    :param input_str: input string representing the nltk style tree
    :return: string representing the tree in TikZ style
    """
    def escape_special_characters_for_latex(input_text):
        # Define a dictionary of special characters and their LaTeX escape sequences
        special_characters = {
            '&': r'\&',
            '$': r'\$',
            '%': r'\%',
            '#': r'\#',
            '_': r'\_',
            '{': r'\{',
            '}': r'\}',
            '~': r'\textasciitilde',
            '^': r'\textasciicircum',
            '\\': r'\textbackslash',
            '<': r'\textless',
            '>': r'\textgreater',
            '|': r'\textbar',
        }

        # Regular expression pattern to match any of the special characters
        pattern = re.compile('|'.join(re.escape(char) for char in special_characters.keys()))

        # Function to replace the matched special characters with their LaTeX escape sequences
        def replace_special_character(match):
            return special_characters[match.group()]

        # Escape the special characters in the input text
        latex_output = pattern.sub(replace_special_character, input_text)

        return latex_output

    nltk_tree = escape_special_characters_for_latex(input_str)
    try:
        # Tokenize the input string and convert it to a tree
        tree = Tree.fromstring(nltk_tree)

        def convert_to_latex_tree(tree):
            result = ''
            if isinstance(tree, nltk.Tree):
                # result += f"[.{tree.label()} "
                result += f"[{tree.label()} "
                for child in tree:
                    result += convert_to_latex_tree(child)
                result += '] '
            else:
                result += f"{tree} "
            return result

        latex_tree = convert_to_latex_tree(tree)

        return latex_tree
    except Exception as e:
        return f"Error: {str(e)}"


# Main function
if __name__ == '__main__':
    args = get_args()

    i = 1
    tot = str(len(args.input_files))

    for file in args.input_files:

        print("(" + str(i) + "/" + tot + ")")
        print("Extracting official data from " + file + "...")

        sents = pd.read_csv(file)
        data = []

        # Read required variables
        for index, row in tqdm(sents.iterrows(), total=len(sents.index)):
            governor_type = row['Governor Type']
            governor = str(row['Governor'])
            governor_side = row['Governor Side']
            conjunction = row['Conjunction']
            con1_cat = row['1st Conjunct Category']
            conjunct1 = str(row['1st Conjunct Text'])
            con2_cat = row['2nd Conjunct Category']
            conjunct2 = str(row['2nd Conjunct Text'])
            sentence = str(row['Sentence Text'])
            parse_tree = row['Sentence Parse Tree']
            con1_pos = ast.literal_eval(row['1st Conjunct Position'])
            con2_pos = ast.literal_eval(row['2nd Conjunct Position'])
            sen_id = str(row['COCA Paragraph Id']) + "-" + str(row['Sentence Id'])
            coca_type = row['COCA Type']
            coca_file = row['COCA File']

            # Extracting length information from left conjunct
            conj1_word_count = word_count(conjunct1)
            conj1_token_count = token_count(con1_pos, parse_tree)
            conj1_syllable_count = syllable_count(conjunct1)
            conj1_char_count = char_count(conjunct1)

            # Extracting length information from right conjunct
            conj2_word_count = word_count(conjunct2)
            conj2_token_count = token_count(con2_pos, parse_tree)
            conj2_syllable_count = syllable_count(conjunct2)
            conj2_char_count = char_count(conjunct2)

            parse_tree_latex = transform_to_latex_tree(parse_tree[1:-1])

            data.append([
                # Governor data
                governor_side,
                governor,
                governor_type,
                # Conjunction data
                conjunction,
                # Left conjunct
                conjunct1,
                con1_cat,
                conj1_word_count,
                conj1_token_count,
                conj1_syllable_count,
                conj1_char_count,
                # Right conjunct
                conjunct2,
                con2_cat,
                conj2_word_count,
                conj2_token_count,
                conj2_syllable_count,
                conj2_char_count,
                # COCA Sentence data
                sentence,
                parse_tree_latex,
                sen_id,
                coca_type,
                coca_file])

        new_columns = ['governor.position',
                       'governor.word',
                       'governor.tag',

                       'conjunction.word',

                       'L.conjunct',
                       'L.head.tag',
                       'L.words',
                       'L.tokens',
                       'L.syllables',
                       'L.chars',

                       'R.conjunct',
                       'R.head.tag',
                       'R.words',
                       'R.tokens',
                       'R.syllables',
                       'R.chars',

                       'sentence',
                       'sentence.tree',
                       'sent.id',
                       'genre',
                       'converted.from.file']
        df = pd.DataFrame(data, columns=new_columns)

        destination_dir = 'results'

        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)

        dest_name = os.path.splitext(file)[-2]
        dest_name = dest_name.replace("_ccp", "_res")
        destination = destination_dir + '/' + dest_name + ".csv"
        df.to_csv(destination, index=False)

        print("Computing final results successful!")
        print("File saved to ", destination)
