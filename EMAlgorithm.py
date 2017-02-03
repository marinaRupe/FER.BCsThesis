
class EMAlgorithm:
    def __init__(self, dataset, initialThetaList, options):
        self.dataset = dataset
        self.options = options
        self.occurrences = self.calculateOccurrences()
        self.currentSolutions = initialThetaList
        self.newSolutions = initialThetaList
        self.noOfThetas = len(self.currentSolutions)  # genomesCount

    def getSolution(self):
        finished = False

        while not finished:
            print("Current solution:" + str(self.currentSolutions))

            h = self.EStep()
            newSolution = self.MStep(h)

            # check if algorithm converges
            finished = True
            for i in range(0, len(self.currentSolutions)):
                if (newSolution[i] - self.currentSolutions[i]) > 0.000001:
                    finished = False
                    break
            if not finished:
                self.currentSolutions = newSolution

        return self.currentSolutions

    def EStep(self):
        # expectation of parameters
        sumOfExpectedNumbers = [[0, 0], [0, 0]]  # for every coin

        for i in range(0, len(self.dataset)):
            noOfHeads = self.occurrences[i]['H']
            noOfTails = self.occurrences[i]['T']

            tempProbabilities = list()
            mixtureDensity = 0
            for j in range(0, self.noOfThetas):
                tempProbability = pow(self.currentSolutions[j], noOfHeads) * pow(1 - self.currentSolutions[j], noOfTails)
                tempProbabilities.append(tempProbability)
                mixtureDensity += tempProbability

            probabilities = list()
            for tempProbability in tempProbabilities:
                probabilities.append(tempProbability / mixtureDensity)

            for j in range(0, len(probabilities)):
                expectedNOOfHeads = noOfHeads * probabilities[j]
                expectedNOOfTails = noOfTails * probabilities[j]

                # expectedNumbers = expectedNOOfHeads, expectedNOOfTails
                sumOfExpectedNumbers[j][0] += expectedNOOfHeads
                sumOfExpectedNumbers[j][1] += expectedNOOfTails

        return sumOfExpectedNumbers

    def MStep(self, h):
        # maximization of parameters
        newSolution = list()

        for hi in h:
            head = hi[0]
            tail = hi[1]

            newParameter = head / (head + tail)
            newSolution.append(newParameter)

        return newSolution

    def calculateOccurrences(self):
        # calculate number of occurrences for every entry
        occurrences = list()
        for x in self.dataset:
            occurrencesForX = dict()

            for option in self.options:
                occurrencesForX[option] = 0

            for xi in x:
                occurrencesForX[xi] += 1
            occurrences.append(occurrencesForX)

        return occurrences


