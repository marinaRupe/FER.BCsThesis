import math
from TaxonomyTree import TaxonomyTree
from threading import Thread
from multiprocessing import Process


class EMAlgorithm:
    def __init__(self):
        self.taxTree = TaxonomyTree()
        self.reads, self.genomes = [], []
        self.pi_list, self.delta_list = [], []
        self.a_list, self.b_list = [], []
        self.q_list = []
        self.y_list = []
        self.groups = {}       # key - TI
        self.parentTIs = {}    # key - TI (get parent TI by organism TI)

    def start(self, alignmentsFile):
        taxTreeThread = Thread(target=self.taxTree.build)
        taxTreeThread.start()

        # First substep
        alignments = self.calculateInitialParameters(alignmentsFile)
        result = self.getResult()
        taxTreeThread.join()
        print("\nFirst result:\n")
        self.printResult(result)
        bestTIs = self.getBestTIsPerGroup(result)

        # Second substep
        self.calculateInitialParameters(alignmentsFile, alignments, bestTIs)
        result = self.getResult()
        print("\nFinal result:\n")
        self.printResult(result)

    def calculateInitialParameters(self, alignmentsFile, alignments={}, bestTIs=[]):
        IS_SECOND_STEP = bool(alignments)
        self.reads, self.a_list, self.b_list, self.y_list = [], [], [], []
        q_dict, map_freq, unique, non_unique = {}, {}, {}, {}
        genomes = set()
        maxScore = 0

        if not IS_SECOND_STEP:
            print("\nSetting the initial parameters...")
            with open(alignmentsFile) as alignFile:
                for line in alignFile:
                    if not line.startswith('@'):
                        fields = line.strip().split("\t")
                        QNAME = fields[0]
                        RNAME = fields[2]  # reference marker gene

                        if RNAME != "*":
                            # calculate score
                            CIGAR = fields[5]
                            matches = n = 0
                            num = ''
                            for i in range(len(CIGAR)):
                                c = CIGAR[i]
                                if c.isdigit():
                                    num += c
                                else:
                                    if c == 'M':
                                        matches += int(num)
                                    n += int(num)
                                    num = ''
                            score = n / matches

                            TIs = RNAME.split("|")[3].split(",")
                            alignments[QNAME] = TIs, score
        else:
            print("Resetting the parameters...")

        for read in alignments:
            self.reads.append(read)
            TIs, score = alignments[read][0], alignments[read][1]

            if score > maxScore:
                maxScore = score

            if IS_SECOND_STEP:
                TIsLeft = []
                for TI in TIs:
                    if TI in bestTIs:
                        TIsLeft.append(TI)
                TIs = TIsLeft

            for TI in TIs:
                genomes.add(TI)
                q_dict[read, TI] = score
                map_freq[TI] = map_freq.get(TI, 0) + 1

                if len(TIs) == 1:
                    unique[TI] = unique.get(TI, 0) + 1
                else:
                    non_unique[TI] = non_unique.get(TI, 0) + 1

            if len(TIs) == 1:
                self.y_list.append(1)
            else:
                self.y_list.append(0)

        self.genomes = list(genomes)
        self.q_list = []
        for i in range(len(self.reads)):
            q_list_for_read = []
            for j in range(len(self.genomes)):
                q_list_for_read.append(math.exp(q_dict.get((self.reads[i], self.genomes[j]), 0) / maxScore))
            self.q_list.append(q_list_for_read)

        self.pi_list, self.delta_list = [], []
        pi0 = delta0 = 1.0 / len(self.genomes)

        for i in range(len(self.genomes)):
            freq = map_freq[self.genomes[i]]
            unique_reads = unique.get(self.genomes[i], 0)
            non_unique_reads = non_unique.get(self.genomes[i], 0)

            self.a_list.append(freq + unique_reads)
            self.b_list.append(freq + non_unique_reads)

            self.pi_list.append(pi0)
            self.delta_list.append(delta0)

        return alignments

    def EStep(self):
        # expectation of parameters
        h = [[0 for j in range(len(self.genomes))] for i in range(len(self.reads))]
        h_sum = 0
        for i in range(len(self.reads)):
            for j in range(len(self.genomes)):
                h[i][j] = self.pi_list[j] * pow(self.delta_list[j], 1 - self.y_list[i]) * self.q_list[i][j]
                h_sum += h[i][j]

        h = [[h[i][j] / h_sum for j in range(len(self.genomes))] for i in range(len(self.reads))]

        return h, h_sum

    def MStep(self, h, N):
        # maximization of parameters
        pi_list = list()
        delta_list = list()
        a_sum = sum(self.a_list)
        b_sum = sum(self.b_list)

        for j in range(len(self.genomes)):
            h_j_sum_by_reads = 0
            h_j_with_y_sum_by_reads = 0
            y_sum = 0

            for i in range(len(self.reads)):
                h_j_sum_by_reads += h[i][j]
                h_j_with_y_sum_by_reads += h[i][j] * (1 - self.y_list[i])
                y_sum += 1 - self.y_list[i]

            pi = self.calculatePi(h_j_sum_by_reads, self.a_list[j], a_sum, N)
            pi_list.append(pi)

            delta = self.calculateDelta(h_j_with_y_sum_by_reads, y_sum, self.b_list[j], b_sum)
            delta_list.append(delta)

        return pi_list, delta_list

    def calculateLogLikelihood(self):
        log_likelihood = 0

        for i in range(len(self.reads)):
            inner_sum = 0
            for j in range(len(self.genomes)):
                inner_sum += self.pi_list[j] * (pow(self.delta_list[j], 1 - self.y_list[i]) * self.q_list[i][j])
            log_likelihood += math.log(inner_sum)

        return log_likelihood

    def getBestTIsPerGroup(self, result):
        print("\nGetting the best TIs per group...")

        for TI in self.genomes:
            genome = [r[0] for r in result if r[1] == TI][0], TI
            parentTI = self.taxTree.taxNodes[TI].parent.taxId
            groupGenome = self.groups.get(parentTI, (0, None))

            if genome[0] > groupGenome[0]:
                self.groups[parentTI] = genome
                self.parentTIs[TI] = parentTI

        bestTIs = []
        for group in self.groups:
            bestTIs.append(self.groups[group][1])

        return bestTIs

    def getResult(self):
        EPSILON = pow(10, -8)
        finished = False
        log_likelihood = None

        while not finished:
            h, N = self.EStep()
            new_pi_list, new_delta_list = self.MStep(h, N)
            new_log_likelihood = self.calculateLogLikelihood()

            convergency_of_log_likelihood = (log_likelihood is not None) and (abs(new_log_likelihood - log_likelihood) < EPSILON)
            log_likelihood = new_log_likelihood

            # check if algorithm converges
            finished = False
            for i in range(len(self.pi_list)):
                conv_pi = abs(new_pi_list[i] - self.pi_list[i]) < EPSILON
                conv_delta = abs(new_delta_list[i] - self.delta_list[i]) < EPSILON
                convergency_of_parameters = conv_pi and conv_delta

                if convergency_of_parameters or convergency_of_log_likelihood:
                    finished = True
                    break

            if not finished:
                self.pi_list, self.delta_list = new_pi_list, new_delta_list

        sum_pi = sum(self.pi_list)
        solution = []
        for i in range(len(self.genomes)):
            solution.append((self.pi_list[i] / sum_pi, self.genomes[i]))

        return sorted(solution, reverse=True)

    def printResult(self, result):
        N = 5
        NO_NAME = "(no name found)"
        names = [NO_NAME for i in range(N)]
        TIs = [genome[1] for genome in result[:N]]

        for TI in TIs:
            if names[TIs.index(TI)] == NO_NAME:
                names[TIs.index(TI)] = self.taxTree.taxonomyNames.get(TI, None)

        for i in range(N):
            print("{}. {}".format(i + 1, names[i]))
            print("     {:10}  {:>.8}".format(result[i][1], result[i][0]))

        return TIs

    @staticmethod
    def calculatePi(h_j_sum_by_R, a_j, a_sum, N):
        return (h_j_sum_by_R + a_j) / (N + a_sum)

    @staticmethod
    def calculateDelta(h_j_with_y_sum_by_R, y_sum, b_j, b_sum):
        return (h_j_with_y_sum_by_R + b_j) / (y_sum + b_sum)
