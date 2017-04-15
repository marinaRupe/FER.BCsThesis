import re
import itertools
from TaxonomyTree import TaxonomyTreeNode

NAMES_FILE = "res/names.dmp"
NODES_FILE = "res/nodes.dmp"
STRAINS_ASSEMBLY_FILE = "res/assembly_summary_genbank_and_refseq.txt"
MARKERS_FILE = "res/markers_info.txt"
NOT_PAIRED_CLADES_FILE = "res/notPairedClades.txt"
CODING_SEQUENCES_FILE = "res/markeri.out"
REDUCED_DB_FILE = "reducedDatabase/reducedDatabase.fa"


class DatabaseReducer:
    def __init__(self):
        self.taxonomyNames = dict()     # key = TI
        self.taxonomyRanks = dict()     # key = rank
        self.taxonomyNodes = dict()     # key = taxName
        self.strainAssemblies = dict()  # key = assembly
        self.strainNames = dict()       # key = taxName
        self.markers = dict()           # key = GI, position

        self.parseTaxonomyNamesFile(NAMES_FILE)
        self.parseStrainsAssemblyFile(STRAINS_ASSEMBLY_FILE)

        self.parseTaxonomyNodesFile(NODES_FILE)
        self.buildTaxonomyTree()

        self.parseMarkersFile(MARKERS_FILE, NOT_PAIRED_CLADES_FILE)
        self.pairMarkers(CODING_SEQUENCES_FILE, REDUCED_DB_FILE)

    def parseTaxonomyNamesFile(self, taxonomyNamesFile):
        print("Preparing taxonomy names...")

        with open(taxonomyNamesFile) as namesFile:
            for line in namesFile:
                taxName = line.split('|')
                taxId = taxName[0].strip()
                name = re.sub('[^0-9a-zA-Z]+', '_', taxName[1].strip())
                self.taxonomyNames[taxId] = name

        print("---Done.")

    def parseTaxonomyNodesFile(self, taxonomyNodesFile):
        print("Preparing taxonomy nodes...")

        with open(taxonomyNodesFile) as nodesFile:
            for line in nodesFile:
                node = line.split('|')
                taxId, parentTaxId, rank = node[0].strip(), node[1].strip(), node[2].strip()
                taxName = self.taxonomyNames[taxId]

                if rank not in self.taxonomyRanks:
                    self.taxonomyRanks[rank] = [(taxId, parentTaxId, taxName)]
                else:
                    self.taxonomyRanks[rank].append((taxId, parentTaxId, taxName))

        print("---Done.")

    def parseStrainsAssemblyFile(self, assemblyFile):
        print("Preparing strains assemblies...")

        with open(assemblyFile) as strainsAssemblyFile:
            for line in strainsAssemblyFile:
                if not line.startswith('#'):
                    data = line.split('\t')
                    assembly = data[0].strip().split('.')[0]
                    gbrs_paired_asm = data[17].strip().split('.')[0]
                    speciesTI, organismName = data[6].strip(), data[7].strip()

                    self.strainAssemblies[assembly] = speciesTI
                    self.strainAssemblies[gbrs_paired_asm] = speciesTI
                    self.strainNames[organismName] = speciesTI

        # manually add not existing assemblies
        self.strainAssemblies['GCF_000219255'] = '1491'
        self.strainAssemblies['GCF_000295795'] = '1225174'
        self.strainAssemblies['GCF_000335395'] = '1262663'
        print("---Done.")

    def buildTaxonomyTree(self):
        print("Building taxonomy tree...")

        for rank in self.taxonomyRanks:
            taxes = self.taxonomyRanks[rank]

            for tax in taxes:
                taxId, parentTaxId, taxName = tax
                parentName = self.taxonomyNames.get(parentTaxId, None)

                # if there's no name for the parent, there is no parent node
                if parentName is None:
                    parentTaxNode = None
                else:
                    parentTaxNode = self.taxonomyNodes.get(parentName, None)

                    # if there's no node for the parent - create it
                    if parentTaxNode is None:
                        parentTaxNode = TaxonomyTreeNode(parentTaxId, None, parentName, None)
                        self.taxonomyNodes[parentName] = parentTaxNode

                taxNode = TaxonomyTreeNode(taxId, parentTaxNode, taxName, rank)

                node = self.taxonomyNodes.get(taxName, "Not existing node")
                if node == "Not existing node":
                    self.taxonomyNodes[taxName] = taxNode
                else:
                    for child in node.children:
                        taxNode.addChild(child)
                    self.taxonomyNodes[taxName] = taxNode

                if parentTaxNode is not None:
                    parentTaxNode.addChild(taxNode)

        self.printNodesStatistic()
        # self.taxonomyNames.clear()
        self.taxonomyRanks.clear()
        print("---Done.")

    def parseMarkersFile(self, markersFile, notPairedCladesFile):
        print("Preparing markers...")

        markersWithNoTIs = 0
        with open(markersFile, 'r') as f1, open(notPairedCladesFile, 'w') as f2:
            for line in f1:
                marker = line.strip().split("\t")
                giInfo = marker[0]

                if giInfo.startswith("gi|"):
                    GI = giInfo.split('|')[1]
                    position = giInfo.split('|')[-1].replace(':', '')
                    info = eval(marker[1])
                    ext = info["ext"]
                    clade = info["clade"].strip()
                    rank, cladeName = clade.split('__')[0], clade.split('__')[1]
                    taxNode = self.taxonomyNodes.get(cladeName, None)

                    if taxNode is not None:
                        # if rank is not species, make clade to species level
                        if rank != 's':
                            speciesTIs = self.getSpecies(taxNode) + self.getTIsFromExt(ext)
                        else:
                            speciesTIs = [taxNode.taxId] + self.getTIsFromExt(ext)

                        self.markers[GI, position] = speciesTIs

                    else:
                        # if clade name is strain assembly
                        if rank == 't':
                            taxId = self.strainAssemblies.get(cladeName, None)
                            if taxId is not None:
                                speciesTIs = [taxId] + self.getTIsFromExt(ext)
                                self.markers[GI, position] = speciesTIs

                            else:
                                f2.write(cladeName + "\n")
                                speciesTIs = self.getTIsFromExt(ext)

                                if len(speciesTIs) > 0:
                                    self.markers[GI, position] = speciesTIs
                                else:
                                    markersWithNoTIs += 1

                        else:
                            # try with strain names
                            taxId = self.strainNames.get(cladeName, None)
                            if taxId is not None:
                                speciesTIs = [taxId] + self.getTIsFromExt(ext)
                                self.markers[GI, position] = speciesTIs

                            else:
                                # TODO: not paired clades
                                f2.write(cladeName + "\n")
                                speciesTIs = self.getTIsFromExt(ext)

                                if len(speciesTIs) > 0:
                                    self.markers[GI, position] = speciesTIs
                                else:
                                    markersWithNoTIs += 1

        self.printMarkersStatistics(markersWithNoTIs)
        self.strainAssemblies.clear()
        self.strainNames.clear()
        self.taxonomyNodes.clear()
        print("---Done.")

    def pairMarkers(self, codingSequencesFile, reducedDatabase):
        print("Pairing markers...")
        paired, notPaired = 0, 0

        with open(codingSequencesFile, 'r') as f, open(reducedDatabase, 'w') as db:
            for line1, line2 in itertools.zip_longest(*[f] * 2):
                marker = line1.strip().split("\t")
                GI = marker[0].split("|")[1]
                position = marker[0].split("|")[4].replace(':', '')

                TIs = self.markers.get((GI, position), None)
                self.markers.pop((GI, position), None)

                # check if there is a name for a taxID

                if (TIs is not None) and len(TIs) != 0:
                    namedTIs = []
                    for taxId in TIs:
                        if self.taxIdHasName(taxId):
                            namedTIs.append(taxId)
                    TIs = namedTIs

                    paired += 1
                    db.write(">gi|" + GI + "|ti|" + ",".join(TIs) + "\n")
                    db.write(line2.strip() + "\n")
                else:
                    notPaired += 1

        self.printPairingMarkersStatistics(paired, notPaired)
        self.taxonomyNames.clear()
        self.markers.clear()
        print("---Done.")

    def getTIsFromExt(self, ext):
        TIs = set()
        for assembly in ext:
            TI = self.strainAssemblies.get(assembly, None)
            if TI is not None:
                TIs.add(TI)
        return list(TIs)

    @staticmethod
    def getSpecies(taxNode):
        # find children on species level
        notSpeciesNodes, speciesNodes = [taxNode], []
        while len(notSpeciesNodes) > 0:
            newNotSpeciesNodes = []
            for node in notSpeciesNodes:
                if node.hasChildren:
                    for child in node.children:
                        if child.rank == "species":
                            speciesNodes.append(child.taxId)
                        else:
                            newNotSpeciesNodes.append(child)
            notSpeciesNodes = newNotSpeciesNodes

        return speciesNodes

    def taxIdHasName(self, taxId):
        name = self.taxonomyNames.get(taxId, None)
        return name is not None

    # print statistics
    def printNodesStatistic(self):
        print("\t\tNames count: " + str(len(self.taxonomyNames)))
        print("\t\tNodes count: " + str(len(self.taxonomyNodes)))
        print("\t\tNames without nodes: " + str(len(self.taxonomyNames) - len(self.taxonomyNodes)))

    def printMarkersStatistics(self, markersWithNoTIsCount):
        if markersWithNoTIsCount > 0:
            print("\t\tMarkers with no TIs: " + str(markersWithNoTIsCount))
        print("\t\tRemaining markers: " + str(len(self.markers)))

    def printPairingMarkersStatistics(self, pairedCount, notPairedCount):
        print("\t\tPaired: " + str(pairedCount))
        print("\t\tNot paired sequences: " + str(notPairedCount))
        print("\t\tNot paired markers: " + str(len(self.markers)))

