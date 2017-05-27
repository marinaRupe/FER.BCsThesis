from pathlib import Path
from os import system
from sys import argv, exit
from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer
from res.ResourceFiles import REDUCED_DB_FILE

try:
    INPUT_FILE = str(argv[1])
except:
    print("Missing argument: input file")
    exit(1)

reducedDbPath = Path(REDUCED_DB_FILE)
if not reducedDbPath.is_file():
    DatabaseReducer().generate()

OUTPUT_FILE_NAME = ''.join(INPUT_FILE.split('/')[-1].split('.')[:-1]) + ".sam"
ALIGNMENTS_FILE = "alignments/out/" + OUTPUT_FILE_NAME

system("graphmap align -r " + "./" + REDUCED_DB_FILE + " -d " + INPUT_FILE + " -o " + "./" + ALIGNMENTS_FILE)
EM = EMAlgorithm().start(ALIGNMENTS_FILE)
