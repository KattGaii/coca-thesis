import sys
import spacy
from nltk import ParentedTree
from find_coord_data import get_coord_data
from benepar import BeneparComponent, NonConstituentException


if __name__ == '__main__':
    # Load spacy model
    print("Loading spaCy's large English model...")
    nlp = spacy.load("en_core_web_lg")
    # Disable sentence segmentation
    nlp.disable_pipes("senter")

    # Integrate with benepar
    print("Integrating spaCy model with Benepar...")
    # print("You may ignore any messages about TensorFlow not being optimized.")
    nlp.add_pipe('benepar', config={'model': 'benepar_en3'})
    print()
    print("\nEnter your sentences and phrases below.")
    print("======================================================================================================================")
    while True:
        line = sys.stdin.readline().strip()
        if not line:
            break
        doc = nlp(line)
        for sent in doc.sents:
            parse_tree = '(' + sent._.parse_string + ')'
            print(parse_tree)
            tree = ParentedTree.fromstring(parse_tree)
            coord = get_coord_data(tree, line)
            for c in coord:
                print("conj 1 cat:\t", c[0][0])
                print("conj 1:\t", c[0][1])
                print("conj 2 cat:\t", c[2][0])
                print("conj 2:\t", c[2][1])
                print("conjunction:\t", c[1])
                print("governor type:\t", c[3][0])
                print("governor:\t", c[3][1])
            tree.draw()
        print("\nEnter your sentences and phrases below.")
        print("======================================================================================================================")

