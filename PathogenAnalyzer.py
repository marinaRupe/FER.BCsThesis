from pathlib import Path
from os import system
from sys import argv, exit
from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer
from res.ResourceFiles import REDUCED_DB_FILE


def main():
    try:
        INPUT_FILE = str(argv[1])
    except:
        print("Missing argument: input file")
        exit(1)

    if not Path(REDUCED_DB_FILE).is_file():
        DatabaseReducer().generate()

    OUTPUT_FILE_NAME = ''.join(INPUT_FILE.split('/')[-1].split('.')[:-1]) + ".sam"
    ALIGNMENTS_FILE = "alignments/out/" + OUTPUT_FILE_NAME

    if not Path(ALIGNMENTS_FILE).is_file():
        system("graphmap align -r " + "./" + REDUCED_DB_FILE + " -d " + INPUT_FILE + " -o " + "./" + ALIGNMENTS_FILE)

    EMAlgorithm().start(ALIGNMENTS_FILE)


if __name__ == "__main__":
    main()
