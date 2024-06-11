"""
find_coord_data.py
Goes through parse trees provided by parseTSV.py and searches for binary coordinations in them.
For each coordination found, searches for the governor of given coordination.
Saves results to *name*_ccp.csv
This code has been inspired by Julie Kallini's senior thesis project
source: https://github.com/jkallini/PrincetonThesis

"""
import argparse
import os

import pandas as pd
from tqdm import tqdm
from nltk import ParentedTree
import coord_gov_BP as coord_BP
import re

# Punctuation labels to ignore
punctuation = ["''", "``", "'", ",", ".", ":", "HYPH", "-LRB-", "-RRB-"]
other_labels = ['$']


def get_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Get coordination stats from csv input file(s) containing parsed sentences.')
    parser.add_argument('input_files', nargs='+', type=str,
                        help='path to input csv file(s)')
    return parser.parse_args()


def get_tree_text(tree):
    """
    Returns the leaves of the tree joined with a white space.
    Mainly used for evaluating coordinators.
    :param tree: the tree to be transformed into a string
    :return: tree leaves joined by whitespaces
    """
    return " ".join(tree.leaves())


def remove_whitespace(text):
    """
    Removes unnecessary whitespaces between words and punctuation marks (ex tokens of the parser)
    :param text: text to be cleared of unnecessary whitespaces.
    :return: cleared text
    """
    # Define regular expression pattern to find unnecessary whitespaces
    contractions = r"\s+('re|n't|'m|'d|'s|'ll|'ve|'RE|N'T|'M|'D|'S|'LL|'VE)(\s*)"
    cannot = r"(can)\s+(not)"
    gonna = r"(gon)\s+(na)"
    punctuations = r"(\w|\)|\]|})\s+([.,?!;:])(\s?)\s*"
    quotes_double = r'(")(\s*)([^"]*[^"\s])(\s*)(")'
    quotes_single = r"(')(\s*)([^']*[^'\s])(\s*)(')"
    parentheses_R = r'(\()(\s*)([^\)]*[^\)\s])(\s*)(\))'
    parentheses_S = r'(\[)(\s*)([^\]]*[^\]\s])(\s*)(\])'
    parentheses_C = r'(\{)(\s*)([^\}]*[^\}\s])(\s*)(\})'

    # Substitute the pattern with the punctuations mark followed by a single whitespace
    corrected_text = re.sub(contractions, r'\1\2', text)  # "\1\2" <-- "\1 "
    corrected_text = re.sub(cannot, r'\1\2', corrected_text)
    corrected_text = re.sub(gonna, r'\1\2', corrected_text)
    corrected_text = re.sub(punctuations, r'\1\2\3', corrected_text)
    corrected_text = re.sub(quotes_double, r'\1\3\5', corrected_text)
    corrected_text = re.sub(quotes_single, r'\1\3\5', corrected_text)
    corrected_text = re.sub(parentheses_R, r'\1\3\5', corrected_text)
    corrected_text = re.sub(parentheses_S, r'\1\3\5', corrected_text)
    corrected_text = re.sub(parentheses_C,  r'\1\3\5', corrected_text)
    return corrected_text


def unescape_parentheses(leaf):
    """
    Changes the sequence of escaped parentheses used in the BNP tree to standard parentheses
    :param leaf: current element of the tree being changed
    :return: a string without escaped parentheses
    """
    leaf = re.sub(r'-RRB-', ')', leaf)
    leaf = re.sub(r'-LRB-', '(', leaf)
    leaf = re.sub(r'-RCB-', '}', leaf)
    leaf = re.sub(r'-LCB-', '{', leaf)
    leaf = re.sub(r'-RSB-', ']', leaf)
    leaf = re.sub(r'-LSB-', '[', leaf)
    return leaf


def tree2sentence_text(leaves, sentence):
    """
    Escapes special characters and joins words with spaces
    :param leaves: the parts of the sentence as they appear in the tree
    :param sentence: the actual sentence to which comparison is made
    :return: string corresponding to the actual sentence
    """
    words = []
    for leaf in leaves:
        words.append(unescape_parentheses(leaf))
    escaped_words = [re.escape(word) for word in words]
    regex_pattern = r'\s*'.join(escaped_words)
    match = re.search(regex_pattern, sentence)
    output = match.group()
    return output


def remove_nested_coordination(phrases, tree):
    """
    Removes coordinations that can be possibly multiple conjunct coordinations.
    It removes the coordination only if the nested coordinator is the same as the original one.
    :param phrases: list of all coordinations and their data
    :param tree: the entire sentence tree
    :return: cleared list of coordinations and their data
    """
    output = phrases
    to_be_removed = []
    for coord in phrases:
        coord_pos = coord[4][0]  # position of coord phrase in tree
        if coord_pos[:-1]:  # check if possible parent
            coord_parent = tree[coord_pos[:-1]]
            coord_lvl_pos = coord_pos[-1]
            for child in coord_parent:
                if child.label() == "CC" and get_tree_text(child).lower() == coord[1]:
                    rel_pos = child.treeposition()[-1] - coord_lvl_pos
                    # Only CC which are direct neighbours (excluding punctuation) of coord are considered nested
                    # coordination
                    if abs(rel_pos) > 2:
                        continue
                    if rel_pos == 2 and (coord_parent[coord_lvl_pos + 1].label() not in punctuation and
                                         coord_parent[coord_lvl_pos + 1].label() not in other_labels):
                        continue
                    if rel_pos == -2 and (coord_parent[coord_lvl_pos - 1].label() not in punctuation and
                                          coord_parent[coord_lvl_pos - 1].label() not in other_labels):
                        continue

                    # Remove coord from phrases if not already removed
                    if coord not in to_be_removed:
                        to_be_removed.append(coord)
                    # Remove coord from phrases if it is in phrases, and it is a parent of coord
                    remove = [c for c in phrases if c[4][0] == coord_pos[:-1]]
                    if remove and remove[0] not in to_be_removed:
                        to_be_removed.append(remove[0])
    for el in to_be_removed:
        output.remove(el)
    return output


def get_coord_data(tree, sent):
    """
    Searches a parse tree looking for coordination.
    :param tree: sentence parse tree
    :param sent: sentence text
    :return: coordination data: 1st conjunct, coordinator, 2nd conjunct, governor,
    coordination position, tree leaves, original sentence
    """
    phrases = []
    for s in tree.subtrees(
            lambda t: t.label() == "CC"):
        coordination = s.parent()

        cc_list = []
        cc_pos = 0
        conj_list = []

        for child in list(coordination):
            # Node is a potential coordinator
            if child.label() == 'CC':
                cc_list.append(get_tree_text(child).lower())
                cc_pos = child.treeposition()
            # Node is a punctuation mark
            elif child.label() in punctuation or child.label() in other_labels:
                continue
            # Node is a potential conjunct
            else:
                conj_list.append([child.label(), tree2sentence_text(child.leaves(), sent), child.treeposition()])

        # If the coordination is not binary skip the phrase
        if len(conj_list) != 2:
            continue

        # If there is only one coordinator, and it is before the 1st conjunct skip the phrase
        if len(cc_list) == 1 and cc_pos < conj_list[0][2]:
            continue

        # Search for the governor
        governor = coord_BP.get_governor(coordination, tree)

        phrases.append((conj_list[0], " ".join(cc_list), conj_list[1], governor,
                        (coordination.treeposition(), tree2sentence_text(coordination.leaves(), sent))))

    results = remove_nested_coordination(phrases, tree)
    return results


# Main funciton
if __name__ == '__main__':
    args = get_args()

    i = 1
    tot = str(len(args.input_files))

    for file in args.input_files:

        print("(" + str(i) + "/" + tot + ")")
        print("Gathering coordination data from " + file + "...")

        sents = pd.read_csv(file)

        # Remove the sentences that were split again, due to BNP error, e.g. "et al. [...]" forcing 2 new sentences.
        sents.drop_duplicates(subset=["COCA Sentence Text"], inplace=True)
        data = []

        for index, row in tqdm(sents.iterrows(), total=len(sents.index)):
            coca_type = row["COCA Type"]
            coca_year = row["COCA Year"]
            coca_file = row["COCA File"]
            coca_para_id = row["COCA Paragraph Id"]
            parse_tree = row["Sentence Parse Tree"]
            sen_id = row["Sentence Id"]
            sent = remove_whitespace(row["COCA Sentence Text"])

            # The extra pair of parenthesis is added soley to find the governor. Later it is removed
            parse_tree = "(" + str(parse_tree) + ")"
            tree = ParentedTree.fromstring(parse_tree)

            for coord in get_coord_data(tree, sent):
                category1 = coord[0][0]
                conjunct1 = coord[0][1]
                pos1 = coord[0][2]
                conjunction = coord[1]
                category2 = coord[2][0]
                conjunct2 = coord[2][1]
                pos2 = coord[2][2]
                coord_text = coord[4][1]
                governor_type = coord[3][0]
                governor_position = coord[3][2]
                gov_pos = "0"
                if governor_position != '-':
                    if governor_position < pos1:
                        gov_pos = "L"
                    elif governor_position > pos2:
                        gov_pos = "R"
                if governor_position != '-':
                    governor = tree2sentence_text(coord[3][1], sent)
                else:
                    governor = "-"

                data.append([governor_type,
                             governor,
                             governor_position,
                             gov_pos,
                             conjunction,
                             category1,
                             conjunct1,
                             pos1,
                             category2,
                             conjunct2,
                             pos2,
                             coord_text,
                             sent,
                             parse_tree,
                             0,
                             sen_id,
                             coca_type,
                             coca_year,
                             coca_file,
                             coca_para_id])
        columns = ['Governor Type',
                   'Governor',
                   'Governor Position',
                   'Governor Side',
                   'Conjunction',
                   '1st Conjunct Category',
                   '1st Conjunct Text',
                   '1st Conjunct Position',
                   '2nd Conjunct Category',
                   '2nd Conjunct Text',
                   '2nd Conjunct Position',
                   'Coordination',
                   'Sentence Text',
                   'Sentence Parse Tree',
                   'Coord. Id',
                   'Sentence Id',
                   'COCA Type',
                   'COCA Year',
                   'COCA File',
                   'COCA Paragraph Id']

        df = pd.DataFrame(data, columns=columns)
        columns_to_check = ['1st Conjunct Category', '1st Conjunct Text',
                            '2nd Conjunct Category', '2nd Conjunct Text',
                            'Conjunction', 'Sentence Parse Tree', 'Sentence Text']
        df.drop_duplicates(subset=columns_to_check, inplace=True)
        df.reset_index(inplace=True, drop=True)
        df["Coord. Id"] = df.index + 1

        # Save to /ccp/file_ccp.csv
        dest_name = os.path.splitext(os.path.basename(file))[-2]
        dest_name = dest_name.replace("_raw", "_ccp") 
        dest_dir = 'ccp'

        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        dest = dest_dir + '/' + dest_name + '.csv'
        df.to_csv(dest, index=False)

        print("All done! The result is stored in " + dest + ".\n")
        i = i + 1
