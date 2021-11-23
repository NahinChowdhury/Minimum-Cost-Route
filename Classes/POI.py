class POI:
    def __init__(self, id, long, lat, eID, name = ""):
        self.id = id
        self.long = long 
        self.lat = lat
        self.name = name
        self.eID = eID
        self.useful = 1 # 1 if we use it
        self.items = [0 for item in range(1000)] # should have N indices and items that can be bought here should have values greater than 0
        self.n1 = None # closest edge first node
        self.n2 = None # closest edge second node
        self.closestNode = None # closest edge

    def setUseless(self):
        self.useful = 0

    # print the id, name, lat, lon
    def __str__(self):
        return str(self.id) + " " + self.name + " " + str(self.lat) + " " + str(self.lon)
    