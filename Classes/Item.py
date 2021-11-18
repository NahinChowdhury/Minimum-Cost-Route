class Item():
    def __init__(self, id, name=""):
        self.id = id
        self.name = name
        self.POIs = []


    def __str__(self):
        return "{}\n=====\n{}\nValue: {}\n".format(self.name, self.id, self.POIs)

    def __repr__(self):
        return self.name
