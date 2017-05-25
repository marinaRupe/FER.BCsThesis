import itertools
import re

from TaxonomyTree import TaxonomyTreeNode
from res.ResourceFiles import NAMES_FILE, NODES_FILE, STRAINS_ASSEMBLY_FILE, MARKERS_FILE, NOT_PAIRED_CLADES_FILE,\
    CODING_SEQUENCES_FILE, REDUCED_DB_FILE, NODES_STATS_FILE


class DatabaseReducer:
    def __init__(self):
        self.taxonomyNames = {}     # key = TI
        self.taxIDFromName = {}     # key = taxName
        self.taxNodes = {}          # key = taxName
        self.strainAssemblies = {}  # key = assembly
        self.strainNames = {}       # key = taxName
        self.markers = {}           # key = GI, position

    def generate(self):
        self.parseTaxonomyNamesFile(NAMES_FILE)
        self.parseStrainsAssemblyFile(STRAINS_ASSEMBLY_FILE)
        self.parseTaxonomyNodesFile(NODES_FILE)
        self.parseMarkersFile(MARKERS_FILE, NOT_PAIRED_CLADES_FILE)
        self.pairMarkers(CODING_SEQUENCES_FILE, REDUCED_DB_FILE)

    def parseTaxonomyNamesFile(self, taxonomyNamesFile):
        print("Preparing taxonomy names...")

        with open(taxonomyNamesFile) as namesFile:
            for line in namesFile:
                taxName = line.split('|')
                TI = taxName[0].strip()
                name = re.sub('[^0-9a-zA-Z]+', '_', taxName[1].strip())
                self.taxonomyNames[TI] = name
                self.taxIDFromName[name] = TI

        print("---Done.")

    def parseTaxonomyNodesFile(self, taxonomyNodesFile):
        print("Preparing taxonomy nodes...")

        with open(taxonomyNodesFile) as nodesFile:
            for line in nodesFile:
                node = line.split('|')
                TI, parentTaxId, rank = node[0].strip(), node[1].strip(), node[2].strip()
                taxName = self.taxonomyNames[TI]
                self.addToTaxonomyTree(TI, parentTaxId, taxName, rank)

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

        print("---Done.")

    def addToTaxonomyTree(self, TI, parentTI, taxName, rank):
        if parentTI in self.taxNodes:
            parentNode = self.taxNodes[parentTI]
        else:
            parentNode = TaxonomyTreeNode(parentTI, None, self.taxonomyNames.get(parentTI, None), None)

        if TI in self.taxNodes:
            taxNode = self.taxNodes[TI]
            if taxNode.rank is None or taxNode.rank == "no rank":
                taxNode.rank = rank
        else:
            taxNode = TaxonomyTreeNode(TI, parentNode, taxName, rank)

        parentNode.addChild(taxNode)
        self.taxNodes[parentTI] = parentNode
        self.taxNodes[TI] = taxNode

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
                    clade, ext = info["clade"].strip(), info["ext"]
                    rank, cladeName = clade.split('__')[0], clade.split('__')[1]

                    cladeTI = self.taxIDFromName.get(cladeName, "not found")
                    taxNode = self.taxNodes.get(cladeTI, None)

                    if taxNode is not None:
                        # if rank is not species, make clade to species level
                        if rank != 's':
                            speciesTIs = self.getSpecies(taxNode) | self.getTIsFromExt(ext)
                        else:
                            speciesTIs = self.getTIsFromExt(ext) | {taxNode.taxId}

                        self.markers[GI, position] = speciesTIs

                    else:
                        # if clade name is strain assembly
                        if rank == 't':
                            TI = self.strainAssemblies.get(cladeName, None)
                            if TI is not None:
                                speciesTIs = self.getTIsFromExt(ext) | {TI}
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
                            TI = self.strainNames.get(cladeName, None)
                            if TI is not None:
                                speciesTIs = self.getTIsFromExt(ext) | {TI}
                                self.markers[GI, position] = speciesTIs

                            else:
                                f2.write(cladeName + "\n")
                                speciesTIs = self.getTIsFromExt(ext)

                                if len(speciesTIs) > 0:
                                    self.markers[GI, position] = speciesTIs
                                else:
                                    markersWithNoTIs += 1

        self.printMarkersStatistics(markersWithNoTIs)
        self.printNodesStatistic()

        self.strainAssemblies.clear()
        self.strainNames.clear()
        self.taxNodes.clear()
        self.taxIDFromName.clear()
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

                if (TIs is not None) and len(TIs) != 0:
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
        return TIs

    @staticmethod
    def getAllChildNodes(taxNode):
        SPECIES_RANK = "species"
        new_child_nodes = taxNode.children
        TIs = set()

        while new_child_nodes:
            child_nodes = new_child_nodes
            new_child_nodes = []
            for child in child_nodes:
                if child.rank == SPECIES_RANK:
                    TIs.add(child.taxId)
                new_child_nodes += child.children

        return TIs

    @staticmethod
    def getSpecies(taxNode):
        SPECIES_RANK = "species"
        GENUS_RANK = "genus"

        # find children on species level
        if taxNode.rank == GENUS_RANK:
            return DatabaseReducer.getAllChildNodes(taxNode)

        speciesTIs = set()
        notSpeciesNodes = [taxNode]
        while len(notSpeciesNodes) > 0:
            newNotSpeciesNodes = []
            for node in notSpeciesNodes:
                if node.rank == GENUS_RANK:
                    speciesTIs.update(DatabaseReducer.getAllChildNodes(node))
                elif node.hasChildren:
                    for child in node.children:
                        if child.rank == SPECIES_RANK:
                            speciesTIs.add(child.taxId)
                        else:
                            newNotSpeciesNodes.append(child)
            notSpeciesNodes = newNotSpeciesNodes

        return speciesTIs

    def taxIdHasName(self, taxId):
        name = self.taxonomyNames.get(taxId, None)
        return name is not None

    # print statistics
    def printNodesStatistic(self):
        print("\t\tNames count: " + str(len(self.taxonomyNames)))
        print("\t\tNodes count: " + str(len(self.taxNodes)))
        print("\t\tNames without nodes: " + str(len(self.taxonomyNames) - len(self.taxNodes)))
        # self.saveTaxNodes()

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

    def printMarkersStatistics(self, markersWithNoTIsCount):
        if markersWithNoTIsCount > 0:
            print("\t\tMarkers with no TIs: " + str(markersWithNoTIsCount))
        print("\t\tRemaining markers: " + str(len(self.markers)))

    def printPairingMarkersStatistics(self, pairedCount, notPairedCount):
        print("\t\tPaired: " + str(pairedCount))
        print("\t\tNot paired sequences: " + str(notPairedCount))
        print("\t\tNot paired markers: " + str(len(self.markers)))
