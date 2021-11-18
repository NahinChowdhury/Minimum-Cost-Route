class Route:
    def __init__(self):
        self.POIs = []
        self.items = []
        self.cost = 0


    def addPOI(self, poi):
        self.POIs.append(poi)

        
    def addItem(self, item):
        self.items.append(item)
        
    def incrementCost(self, cost):
        self.cost += cost