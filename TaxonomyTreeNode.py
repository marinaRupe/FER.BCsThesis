class TaxonomyTreeNode:
    def __init__(self, taxId, parent, name, rank):
        self.taxId = taxId
        self.name = name
        self.rank = rank
        self.parent = parent
        self.children = []

    def addChild(self, child):
        self.children.append(child)

    def setChildren(self, children):
        self.children = children

    def hasChildren(self):
        return len(self.children) > 0
