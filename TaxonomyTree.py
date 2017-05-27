from re import sub
from res.ResourceFiles import NAMES_FILE, NODES_FILE, NODES_STATS_FILE
from TaxonomyTreeNode import TaxonomyTreeNode


class TaxonomyTree:
    def __init__(self, databaseMode=False):
        self.databaseMode = databaseMode
        self.taxonomyNames = {}     # key = TI
        self.taxIDFromName = {}     # key = taxName
        self.taxNodes = {}          # key = taxName

    def build(self):
        printInfo = self.databaseMode
        self.parseTaxonomyNamesFile(NAMES_FILE, printInfo)
        self.parseTaxonomyNodesFile(NODES_FILE, printInfo)

    def parseTaxonomyNamesFile(self, taxonomyNamesFile, printInfo):
        if printInfo:
            print("Preparing taxonomy names...")

        SCIENTIFIC_NAME = "scientific name"
        with open(taxonomyNamesFile) as namesFile:
            for line in namesFile:
                taxName = line.split('|')
                TI = taxName[0].strip()

                if self.databaseMode:
                    name = sub('[^0-9a-zA-Z]+', '_', taxName[1].strip())
                else:
                    name = taxName[1].strip()

                self.taxIDFromName[name] = TI

                if taxName[3].strip() == SCIENTIFIC_NAME and self.taxonomyNames.get(TI, None) is None:
                    self.taxonomyNames[TI] = name
        if printInfo:
            print("---Done.")

    def parseTaxonomyNodesFile(self, taxonomyNodesFile, printInfo):
        if printInfo:
            print("Preparing taxonomy nodes...")

        with open(taxonomyNodesFile) as nodesFile:
            for line in nodesFile:
                node = line.split('|')
                TI, parentTaxId, rank = node[0].strip(), node[1].strip(), node[2].strip()
                taxName = self.taxonomyNames[TI]
                self.addToTaxonomyTree(TI, parentTaxId, taxName, rank)
        if printInfo:
            print("---Done.")

    def addToTaxonomyTree(self, TI, parentTI, taxName, rank):
        if parentTI in self.taxNodes:
            parentNode = self.taxNodes[parentTI]
        else:
            parentNode = TaxonomyTreeNode(parentTI, None, self.taxonomyNames.get(parentTI, None), None)

        if TI in self.taxNodes:
            taxNode = self.taxNodes[TI]
            if taxNode.parent is None:
                taxNode.parent = parentNode
            if taxNode.rank is None or taxNode.rank == "no rank":
                taxNode.rank = rank
        else:
            taxNode = TaxonomyTreeNode(TI, parentNode, taxName, rank)

        parentNode.addChild(taxNode)
        self.taxNodes[parentTI] = parentNode
        self.taxNodes[TI] = taxNode

    def taxIdHasName(self, taxId):
        name = self.taxonomyNames.get(taxId, None)
        return name is not None

    # print statistics
    def printNodesStatistic(self):
        print("\t\tNames count: " + str(len(self.taxonomyNames)))
        print("\t\tNodes count: " + str(len(self.taxNodes)))
        print("\t\tNames without nodes: " + str(len(self.taxonomyNames) - len(self.taxNodes)))

    def saveTaxNodes(self):
        ROOT_TI = '1'
        with open(NODES_STATS_FILE, 'w') as ns:
            node = self.taxNodes[ROOT_TI]
            ns.write(node.taxId + " " + node.name + "\n")
            new_list = node.children
            counter = space_count = 1
            while new_list:
                old_list = new_list
                new_list = []
                for child in old_list:
                    space = " " * space_count
                    ns.write(space + child.taxId + " " + child.name + " (" + child.rank + ")" + "\n")
                    new_list += child.children
                    counter += 1
                space_count += 1
            print("\t\tChildren counter (from organism with TI 1): ", counter)
