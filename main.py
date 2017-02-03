from EMAlgorithm import EMAlgorithm
from DatabaseReducer import DatabaseReducer


def loadInitialData():
    # data
    x1 = "HTTTHHTHTH"
    x2 = "HHHHTHHHHH"
    x3 = "HTHHHHHTHH"
    x4 = "HTHTTTHHTT"
    x5 = "THHHTHHHTH"

    dataset = [x1, x2, x3, x4, x5]

    # initial guess of the parameters
    thetaA0 = 0.60
    thetaB0 = 0.50

    thetaList = [thetaA0, thetaB0]

    options = ['H', 'T']

    return dataset, thetaList, options


def printResult(result):
    print("\nFinal result: ")
    for theta in result:
        print(round(theta, 2))

dataset, initialThetaList, options = loadInitialData()
dbReducer = DatabaseReducer()

# EM = EMAlgorithm(dataset, initialThetaList, options)
# printResult(EM.getSolution())


