from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer

INPUT_FILE = "alignments/alignments_staphi_x10.sam"

# dbReducer = DatabaseReducer()
EM = EMAlgorithm().start(INPUT_FILE)


