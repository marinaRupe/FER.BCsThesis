from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer

STAPHILOCCOCUS = "alignments/alignments_staphi"
SALMONELLA = "alignments/alignments_salmonella"
KLEBSIELLA = "alignments/alignments_klebsiella"

STAPHILOCCOCUS_NEW = "alignments/new/alignments_staphi"
SALMONELLA_NEW = "alignments/new/alignments_salmonella"
KLEBSIELLA_NEW = "alignments/new/alignments_klebsiella"

STAPHILOCCOCUS_TEST = "alignments/test/outs/staphi"
SALMONELLA_TEST = "alignments/test/outs/salmonella"
KLEBSIELLA_TEST = "alignments/test/outs/klebsiella"

INPUT_FILE = SALMONELLA + "_x1.sam"

#DatabaseReducer().generate()
EM = EMAlgorithm().start(INPUT_FILE)


