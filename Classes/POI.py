class POI:
    def __init__(self, id, lat, lon, name):
        self.id = id
        self.lat = lat
        self.lon = lon 
        self.name = name
        self.useful = 1 # 1 if we use it
        self.items = [0 for item in range(1000)] # should have N indices and items that can be bought here should have values greater than 0


    def setUseless(self):
        self.useful = 0

    # print the id, name, lat, lon
    def __str__(self):
        return str(self.id) + " " + self.name + " " + str(self.lat) + " " + str(self.lon)
    