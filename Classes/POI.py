class POI:
    def __init__(self, id, lat, lon):
        self.id = id
        self.lat = lat
        self.lon = lon 
        self.useful = 1 # 1 if we use it
        self.items = [] # should have N indices and items that can be bought here should have values greater than 0


    def setUseless(self):
        self.useful = 0

        
    