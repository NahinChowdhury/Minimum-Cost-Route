class Route:
    def __init__(self):
        self.POIs = []
        self.items = []
        self.costs = []
        self.latestPOI = None
        self.latestItem = None
        self.latestCost = None

        self.totalCost = 0
        self.n = []

    def __str__(self):
        return "poi: "+ str(self.POIs) + ", items: " + str(self.items) + ", costs: " + str(self.costs) + ", totalCost: " + str(self.totalCost) + ", n: " + str(self.n)
    
    def addPOI(self, poi):
        self.POIs.append(poi)
        self.latestPOI = poi

        self.n.append(0)

    def addItem(self, item):
        self.items.append(item)
        self.latestItem = item

    def addCost(self, cost):
        self.costs.append(cost)
        self.incrementTotalCost(cost) # incrementTotalCost updates latestCost for me
        
    def addFinalCost(self, cost):
        self.totalCost += cost

    def incrementTotalCost(self, cost):
        self.totalCost += cost
        self.latestCost = cost

    def removeLatestPOI(self):
        self.POIs.pop()
        self.items.pop()
        cost = self.costs.pop()
        self.totalCost -= cost
        self.n.pop()

        if(len(self.POIs) == 0):
            self.latestPOI = None
            self.latestItem = None
            self.latestCost = None
            return

        self.latestPOI = self.POIs[-1]
        self.latestItem = self.items[-1]
        self.latestCost = self.costs[-1]
        return

    def updateN(self, index, value = 1):
        # if len(self.POIs) - 1 < index:
        #     self.n.append(0)
        #     return

        if (index - len(self.POIs) - 1) >= 2:
            print("There is a gap of more than 2: {}".format(index- len(self.POIs)))


        self.n[index] += value
        return 