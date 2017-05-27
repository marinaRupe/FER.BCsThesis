from itertools import zip_longest
from TaxonomyTree import TaxonomyTree
from res.ResourceFiles import STRAINS_ASSEMBLY_FILE, MARKERS_FILE, NOT_PAIRED_CLADES_FILE,\
    CODING_SEQUENCES_FILE, REDUCED_DB_FILE


class DatabaseReducer:
    def __init__(self):
        self.taxTree = TaxonomyTree(databaseMode=True)
        self.strainAssemblies = {}  # key = assembly
        self.strainTIByName = {}    # key = taxName
        self.markers = {}           # key = GI, position

    def generate(self):
        self.taxTree.build()
        self.parseStrainsAssemblyFile(STRAINS_ASSEMBLY_FILE)
        self.parseMarkersFile(MARKERS_FILE, NOT_PAIRED_CLADES_FILE)
        self.pairMarkers(CODING_SEQUENCES_FILE, REDUCED_DB_FILE)

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
                    self.strainTIByName[organismName] = speciesTI

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
                    clade, ext = info["clade"].strip(), info["ext"]
                    rank, cladeName = clade.split('__')[0], clade.split('__')[1]

                    cladeTI = self.taxTree.taxIDFromName.get(cladeName, "not found")
                    taxNode = self.taxTree.taxNodes.get(cladeTI, None)

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
                            TI = self.strainTIByName.get(cladeName, None)
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
        self.taxTree.printNodesStatistic()

        self.strainAssemblies.clear()
        self.strainTIByName.clear()
        self.taxTree.taxNodes.clear()
        self.taxTree.taxIDFromName.clear()
        print("---Done.")

    def pairMarkers(self, codingSequencesFile, reducedDatabase):
        print("Pairing markers...")
        paired, notPaired = 0, 0

        with open(codingSequencesFile, 'r') as f, open(reducedDatabase, 'w') as db:
            for line1, line2 in zip_longest(*[f] * 2):
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
        self.taxTree.taxonomyNames.clear()
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

    def printMarkersStatistics(self, markersWithNoTIsCount):
        if markersWithNoTIsCount > 0:
            print("\t\tMarkers with no TIs: " + str(markersWithNoTIsCount))
        print("\t\tRemaining markers: " + str(len(self.markers)))

    def printPairingMarkersStatistics(self, pairedCount, notPairedCount):
        print("\t\tPaired: " + str(pairedCount))
        print("\t\tNot paired sequences: " + str(notPairedCount))
        print("\t\tNot paired markers: " + str(len(self.markers)))
