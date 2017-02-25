import math


class EMAlgorithm:
    def __init__(self):
        self.reads = list()
        self.genomes = list()
        self.pi_list = list()
        self.delta_list = list()
        self.q_list = list()
        self.a_list = list()
        self.b_list = list()
        self.y_list = list()

    def start(self, alignmentsFile):
        self.calculateInitialParameters(alignmentsFile)
        result = self.getResult()
        self.printResult(result)

    def calculateInitialParameters(self, alignmentsFile):
        self.reads = list()
        self.a_list = list()
        self.b_list = list()
        self.y_list = list()
        genomes = set()
        alignments = dict()
        q_dict = dict()
        max_MAPQ = 0
        map_freq = dict()
        sum_map_freq = 0

        with open(alignmentsFile) as alignFile:
            for line in alignFile:
                fields = line.strip().split("\t")
                QNAME = fields[0]
                RNAME = fields[2]  # reference marker gene
                MAPQ = int(fields[4])  # mapping quality (if 255 - not avaliable)

                # if read is paired
                if RNAME != "*":
                    TIs = RNAME.split("|")[3].split(",")

                    alignments[QNAME] = list(TIs)
                    self.reads.append(QNAME)

                    for TI in TIs:
                        genomes.add(TI)
                        map_freq[TI] = map_freq.get(TI, 0) + 1
                        sum_map_freq += 1

                    if len(TIs) == 1:
                        self.y_list.append(1)
                    else:
                        self.y_list.append(0)

                    for TI in TIs:
                        q_dict[QNAME, TI] = MAPQ

                    if MAPQ > max_MAPQ:
                        max_MAPQ = MAPQ

                else:
                    alignments[QNAME] = list()

        self.genomes = list(genomes)
        self.q_list = list()

        for i in range(len(self.reads)):
            q_list_for_read = list()
            for j in range(len(self.genomes)):
                q_list_for_read.append(math.exp(q_dict.get((self.reads[i], self.genomes[j]), 0) / max_MAPQ))
            self.q_list.append(q_list_for_read)

        # TODO:
        self.pi_list = list()
        self.delta_list = list()
        avg_map_freq = sum_map_freq / len(map_freq)
        for i in range(len(self.genomes)):
            freq = map_freq[self.genomes[i]]
            self.a_list.append(freq)

            if freq < avg_map_freq:
                self.b_list.append(avg_map_freq - freq)
            else:
                self.b_list.append(0)

            self.pi_list.append(100)
            self.delta_list.append(100)

    def getResult(self):
        finished = False
        while not finished:
            print("Current solution:" + str(self.pi_list[:20]))
            h = self.EStep()
            new_pi_list, new_delta_list = self.MStep(h)

            # check if algorithm converges
            finished = False
            for i in range(len(self.pi_list)):
                conv_pi = abs(new_pi_list[i] - self.pi_list[i]) < pow(10, -10)
                conv_delta = abs(new_delta_list[i] - self.delta_list[i]) < pow(10, -10)
                convergency_of_parameters = conv_pi and conv_delta

                if convergency_of_parameters:
                    finished = True
                    break

            if not finished:
                self.pi_list = new_pi_list
                self.delta_list = new_delta_list

        sum_pi = sum(self.pi_list)
        solution = list()
        for i in range(len(self.genomes)):
            solution.append((self.pi_list[i] / sum_pi, self.genomes[i]))

        return sorted(solution)

    @staticmethod
    def printResult(result):
        names = ["" for i in range(5)]
        TIs = [genome[1] for genome in result[:5]]

        print("\nFinal result:\n")
        with open("names.dmp") as namesFile:
            for line in namesFile:
                taxName = line.split('|')
                TI = taxName[0].strip()
                name = taxName[1].strip()

                if TI in TIs:
                    if not names[TIs.index(TI)]:
                        names[TIs.index(TI)] = name

        for i in range(5):
            print("{}. {}".format(i + 1, names[i]))
            print("     {:10}  {:>.8}".format(result[i][1], result[i][0]))

    def EStep(self):
        # expectation of parameters
        h = [[0 for j in range(len(self.genomes))] for i in range(len(self.reads))]
        h_sum = 0
        for i in range(len(self.reads)):
            for j in range(len(self.genomes)):
                h[i][j] = self.pi_list[j] * pow(self.delta_list[j], 1 - self.y_list[i]) * self.q_list[i][j]
                h_sum += h[i][j]

        h = [[h[i][j] / h_sum for j in range(len(self.genomes))] for i in range(len(self.reads))]

        return h

    def MStep(self, h):
        # maximization of parameters
        pi_list = list()
        delta_list = list()
        a_sum = sum(self.a_list)
        b_sum = sum(self.b_list)

        for j in range(len(self.genomes)):
            h_j_sum_by_reads = 0
            h_j_with_y_sum_by_reads = 0
            y_sum = 0  # TODO: izdvoji

            for i in range(len(self.reads)):
                h_j_sum_by_reads += h[i][j]
                h_j_with_y_sum_by_reads += h[i][j] * (1 - self.y_list[i])
                y_sum += 1 - self.y_list[i]

            pi = self.calculatePi(h_j_sum_by_reads, self.a_list[j], a_sum)
            pi_list.append(pi)

            delta = self.calculateDelta(h_j_with_y_sum_by_reads, y_sum, self.b_list[j], b_sum)
            delta_list.append(delta)

        return pi_list, delta_list

    def calculateLogLikelihood(self):
        for i in range(len(self.reads)):
            for j in range(len(self.genomes)):
                return self.pi_list[j]  # TODO

    def calculatePi(self, h_j_sum_by_R, a_j, a_sum):
        N = len(self.reads)
        return (h_j_sum_by_R + a_j) / (N + a_sum)

    @staticmethod
    def calculateDelta(h_j_with_y_sum_by_R, y_sum, b_j, b_sum):
        return (h_j_with_y_sum_by_R + b_j) / (y_sum + b_sum)
