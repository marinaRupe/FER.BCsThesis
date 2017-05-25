from EMAlgorithm import EMAlgorithm
from alignments import AlignmentsFiles as aligFile

INPUT_FILE = aligFile.KLEBSIELLA_x1

# DatabaseReducer().generate()
EM = EMAlgorithm().start(INPUT_FILE)
