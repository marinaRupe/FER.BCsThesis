from pathlib import Path
import sys
from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer
from alignments import AlignmentsFiles as aligFile
from res.ResourceFiles import REDUCED_DB_FILE

if len(sys.argv) > 1:
    INPUT_FILE = str(sys.argv[1])
else:
    INPUT_FILE = aligFile.KLEBSIELLA_x1

reducedDb = Path(REDUCED_DB_FILE)
if not reducedDb.is_file():
    DatabaseReducer().generate()

ALIGNMENTS_FILE = ''.join(INPUT_FILE.split('.')[:-1]) + ".sam"

#graphmap align -r REDUCED_DB_FILE -d INPUT_FILE -o ALIGNMENTS_FILE

EM = EMAlgorithm().start(ALIGNMENTS_FILE)
