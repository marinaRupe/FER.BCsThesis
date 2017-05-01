from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer

STAPHILOCCOCUS = "alignments/alignments_staphi"
SALMONELLA = "alignments/alignments_salmonella"
KLEBSIELLA = "alignments/alignments_klebsiella"

STAPHILOCCOCUS_NEW = "alignments/new/alignments_staphi"
SALMONELLA_NEW = "alignments/new/alignments_salmonella"
KLEBSIELLA_NEW = "alignments/new/alignments_klebsiella"

STAPHILOCCOCUS2 = "alignments/Ivan/outs/Staphi_x1.out"
SALMONELLA2 = "alignments/Ivan/outs/salmonella_x1.out"
KLEBSIELLA2 = "alignments/Ivan/outs/klebsiella_x1.fastq.out"

INPUT_FILE = KLEBSIELLA + "_x1.sam"

#DatabaseReducer().generate()
EM = EMAlgorithm().start(INPUT_FILE)


