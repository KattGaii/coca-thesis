"""
Changed:
    The method works by adding the new set of parenthesis around the tree to create a new empty label inspired by
    PTB& tree structure.
    For no governor a touple ('-', '-', '-') is returned
This code was originally created by Michał Woźniak and was modified for the purpouse of using the BNP.
"""
import sys


def get_parent(tree, t):
    return t[tree.treeposition()[:-1]]


def get_governor(tree, t):
    """
    Finds the governor and returns its position, value. If it doesn't exist returns ("-", "-", "-").
    :param tree: a subtree of t containing coordination. It's the parent tree of coordination found in get_simple_coordinaion
    :param t: the entire parse tree (nltk tree) of a sentence with added parenthesis at the ends
    :return: a touple
    """
    par = get_parent(tree, t)
    coord_pos = tree.treeposition()[-1]
    head = get_head(par, t, coord_pos)[0]
    return head


def get_head(tree, t, coord_pos, is_1st_level=True):
    """
    Function finding the governor of a coordination.
    :param tree: a subtree containing coordination
    :param t: a whole tree containing 'tree'
    :param coord_pos: position in tree
    :param is_1st_level: if governor is at the level of coordination
    :return:
    """
    # Defined priority groups for 5 categories of phrases
    head_dict = {'V': [['VB', 'VBD', 'VBG', 'VBN', 'VBP', 'VBZ'], ['MD'], ['TO']],
                 'P': [['IN', 'TO'], ['VBG'], ['NP', 'NP-ADV']],
                 'ADJP': [['JJ', 'JJR', 'JJS']],
                 'ADVP': [['ADV', 'RB', 'RBR', 'RBS']],
                 'S': [['TO', 'IN'], ['VP']]  # CHANGE: deleted 'WH' from 1st group
                 }
    # Possible labels of heads of NPs
    n_list = ['NP', 'NN', 'NNS', 'NNP', 'NNPS', 'NML', 'NAC', 'NX']

    head = []

    label = tree.label()

    if label:
        # 4 structures of NPs are taken into consideration
        if label[0] == 'N':
            # Governor immediately after coordination
            if len(tree) > coord_pos + 1:
                p_head = tree[coord_pos + 1]
                if p_head.label() in n_list and tree[coord_pos].label()[0] in ['N', 'U', 'A', 'Q', 'P', 'S']:
                    head.append((p_head.label(), p_head.leaves(), p_head.treeposition()))

            # Governor immediately before coordination
            if coord_pos > 0:
                p_head = tree[coord_pos - 1]
                # CHANGE: added case for S
                if p_head.label() in n_list and (tree[coord_pos].label()[0] in ['N', 'U', 'A', 'Q', 'P', 'V', 'S']):
                    head.append((p_head.label(), p_head.leaves(), p_head.treeposition()))

            # Governor + one of the punctuations below + coordination
            if coord_pos > 1:
                p_head = tree[coord_pos - 2]
                # CHANGE: added case for S
                if tree[coord_pos - 1].label() in [":", ",", ";", "“"] and p_head.label() in n_list and \
                        (tree[coord_pos].label()[0] in ['N', 'U', 'A', 'Q', 'P', 'V', 'S']):
                    head.append((p_head.label(), p_head.leaves(), p_head.treeposition()))

            # Coordination + one of the phrases below + governor
            if len(tree) > coord_pos + 2:
                p_head = tree[coord_pos + 2]
                if p_head.label() in n_list and tree[coord_pos].label()[0] in ['N', 'U', 'A', 'Q', 'P', 'V', 'S'] and \
                        tree[coord_pos + 1].label() in ["JJ", "NP", "PP"]:
                    head.append((p_head.label(), p_head.leaves(), p_head.treeposition()))

        else:
            # For S first two characters of the label are checked because there are many subcategories beggining with "WH"
            if label[0] == 'S':
                for group in head_dict['S']:
                    for child in tree:
                        if child.label()[:2] in group and child.treeposition()[-1] != coord_pos:
                            head.append((child.label(), child.leaves(), child.treeposition()))
            else:
                if label[0] in head_dict.keys():
                    for group in head_dict[label[0]]:
                        for child in tree:
                            if child.label() in group and child.treeposition()[-1] != coord_pos:
                                head.append((child.label(), child.leaves(), child.treeposition()))
                else:
                    if label[:4] in head_dict.keys():
                        for group in head_dict[label[:4]]:
                            for child in tree:
                                if child.label() in group and child.treeposition()[-1] != coord_pos:
                                    head.append((child.label(), child.leaves(), child.treeposition()))

    # No governor found
    if head == []:
        # end of tree - governor doesn't exist
        if tree.label() == '':
            head.append(("-", "-", "-"))

        # governor can be higher up the tree
        else:
            # The function recursively checks higher nodes of the tree if the governor hasn't been found
            head = get_head(get_parent(tree, t), t, tree.treeposition()[-1], False)

    return head


# function transforming a tree into a sentence which can be written in LaTeX file
def latex_sentence(tree):
    leaves = tree.leaves()
    pos = tree.treepositions('leaves')
    sen = ""
    punct = [",", ".", ";", ":", "...", "'s", "n't", "'re", "''", ")", "'"]
    slash = ["$", "%", "&", "#"]

    for i, leaf in enumerate(leaves):
        lab = tree[pos[i][:-1]].label()
        if lab[0] == "-" and lab[-1] == "-" and lab != "-LRB-" and lab != "-RRB-":
            continue
        if lab == "-LRB-":
            leaf = "("
        if lab == "-RRB-":
            leaf = ")"
        if leaf in punct:
            sen += leaf
        else:
            if leaf == "%":
                sen += "\\%"
            else:
                word = ""
                for char in leaf:
                    if char in slash:
                        word += "\\"
                    word += char
                if sen != "":
                    if sen[-1] == "`" or sen[-1] == "(":
                        sen = sen + word
                    else:
                        sen = sen + " " + word
                else:
                    sen = word

    return sen


# function transforming a tree into a sentence
def sentence(tree):
    leaves = tree.leaves()
    pos = tree.treepositions('leaves')
    sen = ""
    punct = [",", ".", ";", ":", "...", "'s", "n't", "'re", "''", "'", ")", "%"]
    for i, leaf in enumerate(leaves):
        lab = tree[pos[i][:-1]].label()
        if lab[0] == "-" and lab[-1] == "-" and lab != "-LRB-" and lab != "-RRB-":
            continue
        if lab == "-LRB-":
            leaf = "("
        if lab == "-RRB-":
            leaf = ")"
        if leaf in punct:
            sen += leaf
        else:
            if sen != "":
                if sen[-1] == "`" or sen[-1] == "(":
                    sen = sen + leaf
                else:
                    sen = sen + " " + leaf
            else:
                sen = leaf

    return sen


if __name__ == '__main__':
    print("Paste parse tree string to evaluate it:")
    from nltk import ParentedTree
    from find_coord_data import get_coord_data

    line = "((S (NP (NP (DT The) (NNS students)) (CC and) (NP (DT the) (NN audience))) (VP (VBP break) (PP (IN into) (NP (NN applause)))) (. .)))"
    tree = ParentedTree.fromstring(line)
    coord = get_coord_data(tree, sentence(tree))
    for c in coord:
        print("conj 1 cat:\t", c[0][0])
        print("conj 1:\t", c[0][1])
        print("conj 2 cat:\t", c[2][0])
        print("conj 2:\t", c[2][1])
        print("conjunction:\t", c[1])
        print("governor type:\t", c[3][0])
        print("governor:\t", c[3][1])
    tree.draw()
    do_input = False
    while do_input:
        line = sys.stdin.readline().strip()
        if not line:
            break
        tree = ParentedTree.fromstring(line)
        coord = get_coord_data(tree, sentence(tree))
        for c in coord:
            print("conj 1 cat:\t", c[0][0])
            print("conj 1:\t", c[0][1])
            print("conj 2 cat:\t", c[2][0])
            print("conj 2:\t", c[2][1])
            print("conjunction:\t", c[1])
            print("governor type:\t", c[3][0])
            print("governor:\t", c[3][1])
        tree.draw()
