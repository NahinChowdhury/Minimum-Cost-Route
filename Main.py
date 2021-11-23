from Classes.Item import Item
from Classes.POI import POI
from Classes.Route import Route

import networkx as nx
import osmnx as ox
import json
import math
import time
from shapely.geometry import Point

"""
Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks. Computers, Environment and Urban Systems 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004
"""
def main():
    # Create a route
    # route = Route()
    # # Create a POI
    # poi = POI(1, 1.11, 3.22)
    # # Create an item
    # item = Item(1, "Apple")
    # # Add the POI to the route
    # route.addPOI(poi)
    # # Add the item to the POI
    # poi.setUseless()

    # route.POIs[0].id = 15
    # # Print the route
    # print(route.POIs[0].id)


    #G = ox.graph_from_place("Edmonton, CA", network_type="drive")
    # fig, ax = ox.plot_graph(G)
    #print(type(G))


    # print(G.nodes)
    # print("Below is items")
    # nodes = dict(G.nodes.items())
    # edges = dict(G.edges.items())
    # print(edges)

    # filters the original node we got from the graph
    # filteredNodes = []
    # for key in nodes.keys():
    #     res = [key, nodes[key]['y'], nodes[key]['x']]
    #     filteredNodes.append(res)

    # filters the original node we got from the graph
    # print(edges[(29851773, 277619420, 0)])
    # filteredEdges = []
    # for key in edges.keys():
    #     print(key)
    #     res = [key[0], key[1], edges[key]['length']]
    #     filteredEdges.append(res)

    # # print(filteredNodes)
    # print(filteredEdges)

    # print("Below is nodes")
    # gdf_nodes, gdf_edges = ox.graph_to_gdfs(G)
    # gdf_nodes.head()
    
    # gdf_edges = gdf_edges.to_dict('list')
    # print(type(gdf_edges))
    # print(gdf_edges)

    # print(len(gdf_edges["osmid"]))
    # print(len(gdf_edges["lanes"]))
    # print(len(gdf_edges["ref"]))
    # print(len(gdf_edges["highway"]))
    # print(len(gdf_edges["maxspeed"]))
    # print(len(gdf_edges["oneway"]))
    # print(len(gdf_edges["length"]))
    # print(len(gdf_edges["name"]))


    # gdf_nodes = gdf_nodes.to_dict('list')
    # print(type(gdf_nodes))
    # print(gdf_nodes)


    # print("Below is edges")
    # print(type(gdf_edges))

    # place = "Bunker Hill, Los Angeles, California"
    # tags = {"building": True}
    # gdf = ox.geometries_from_place(place, tags)
    # gdf.shape
    # fig, ax = ox.plot_footprints(gdf, figsize=(3, 3))
    # print(gdf)

    # used to create a list of nodes(id, y, x) in edmonton
    # with open('edges.txt', 'w') as f:
    #     for edge in filteredEdges:
    #         f.write("{} {} {}\n".format(edge[0], edge[1], edge[2]))

    # f = open('./Data/stores.json',)
    # stores = json.load(f)
    # print(stores["elements"][0]["tags"]["name"])
    # f.close()

    # with open('stores.txt', 'w') as f:
    #     for i in range(len(stores["elements"])):
    #         name = stores["elements"][i]['tags']['name'] # one shop's name was missing, so I put Unknown there
    #         f.write("{} {} {} {}\n".format(stores["elements"][i]['id'], stores["elements"][i]['lat'], stores["elements"][i]['lon'], name))

    # f.close()

        # create a 2d array calculating the distance between each pois
    # sp = [[0 for x in range(len(pois))] for y in range(len(pois))]
    # for i in range(len(pois)):
    #     for j in range(len(pois)):
    #         if i != j:
    #             sp[i][j] = calculateDistance(pois[i], pois[j])

    # point = (pois[0].lat, pois[0].lon)

    # u, v, k, edge_geom, dist = ox.distance.get_nearest_edge(G, point, return_geom=True, return_dist=True)

    # # create shapely point geometry object as (x, y), that is (lng, lat)
    # point_geom = Point(reversed(point))

    # # use shapely to find the point along the edge that is closest to the reference point
    # nearest_point_on_edge = edge_geom.interpolate(edge_geom.project(point_geom))
    # nearest_point_on_edge.coords[0]
    # # dist = ox.distance.shortest_path(G, pois[0].id, pois[1].id, weight='length', cpus=1)
    # print(nearest_point_on_edge.coords[0])
    # print(nearest_point_on_edge)
    # print(u)
    # print(v)

    travelWeight = 0.2
    costWeight = 1 - travelWeight

    graph = nx.Graph()
    
    edges = {}
    readGraph(graph, edges, "./datasets/Amsterdam/roadnetwork/RoadVerticesAMS.txt", "./datasets/Amsterdam/roadnetwork/RoadEdgesAMS.txt")
    # print(edges[107365])
    # print("graph completed")

    pois = []
    readPOIs(pois, graph, edges, "./datasets/Amsterdam/poi/originals/PoiAMS50.txt")

    items = []
    readItems(items, "./datasets/Amsterdam/poi/originals/Item_Cost/StoresPerItemAMS_Cost_INward.txt")
    # print(len(pois[0].items))

    readStoreItems(pois, "./datasets/Amsterdam/poi/originals/Item_Cost/ItemsPerStoreAMS_Cost_INward.txt")

    sp = [[0 for j in range(51)] for i in range(51)] # sp has distance between each POI
    # print("poi len: {}".format(len(pois)))
    readSP(sp, len(pois), "./datasets/Amsterdam/poi/originals/ShortestPathPoi50.txt")
    #print(sp)
    # print(pois[0].closestNode)

    startNode = graph.nodes[8020]
    # calculate distance between each POI and startNode
    startToPOIs = []
    calculateLocationToPOIsDist(graph, startNode, pois, startToPOIs, sp, 1)
    
    endNode = graph.nodes[88896]
    # calculate distance between each POI and endNode
    endToPOIs = []
    calculateLocationToPOIsDist(graph, endNode, pois, endToPOIs, sp, travelWeight)

    # have a 3d array for start to each POI and buying items
    startToPOIItemArray = [[[None, None] for i in range(1000)] for i in range(len(pois))]
    createStart3DArray(startToPOIItemArray, startToPOIs, pois, travelWeight, costWeight) # indices with value null means that the item cannot be bought from that POI
    startToPOIItemArray0_sorted = (sorted(startToPOIItemArray[50], key=lambda x: float('inf') if x[1] is None else x[1])) # this is how we sort rows in a 2d array
    # need to pass a single row to sort the row above

    # maybe I can sort every single row of items so that our algorithm can find the cheapest item to buy faster
    


    # sort startToPOIItemArray by the cost of the item to be bought
    # startToPOIItemArray = sorted(startToPOIItemArray, key=lambda x: x[1])
    # print(startToPOIItemArray0_sorted)
    # print(startToPOIItemArray[1])
    # print(startToPOIItemArray[1])

    # create a 3D array of 51 POIs with 51 POIs with 999 items with default value of None
    currentToNextItemArray = [[[[None, None] for j in range(1000)] for i in range(len(pois))] for k in range(len(pois))]
    create3DArray(currentToNextItemArray, sp, pois, travelWeight, costWeight)
    # print(len(currentToNextItemArray))
    # print(currentToNextItemArray[0])
    # print(len(currentToNextItemArray[0][0]))

    # our last POIS to ending location is just travelWeight times the travel distance
    # i do that by passing the travel weight while calculating the distance

    # now we can implement the algorithm






def readGraph(graph, edges, vertexLocation, edgeLocation):
    with open(vertexLocation, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()
            graph.add_node(int(line[0]), long = float(line[1]), lat = float(line[2]), poi = 0)
        # print(graph.nodes[0])

    f.close()

    with open(edgeLocation, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()
            graph.add_edge(int(line[1]), int(line[2]))
            edges[int(line[0])] = [int(line[1]), int(line[2])]

    # print(graph.edges[64061, 64106])
    f.close()

def readSP(sp, len, location):

    row = 0
    with open(location, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()
            for col in range(len):
                sp[row][col] = float(line[col])
            row += 1
            if row == len:
                break

    f.close()

def readPOIs(pois, graph, edges, location):
    i = 0
    # we should only take 51 stores because items are stored in 51 stores only
    with open(location, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()

            id = int(line[0])
            long = float(line[1])
            lat = float(line[2])
            eID = int(line[3])

            poi = POI(id, long, lat, eID)
            poi.n1 = edges[eID][0]
            poi.n2 = edges[eID][1]

            distN1 = calculateDistance(graph.nodes[poi.n1], poi)
            distN2 = calculateDistance(graph.nodes[poi.n2], poi)

            if distN1 < distN2:
                poi.closestNode = poi.n1
            else:
                poi.closestNode = poi.n2

            pois.append(poi)

    f.close()

    # print(len(pois))

    # for i in pois:
    #     # print(str(i.id) + " " + i.name + " " + str(i.lat) + " " + str(i.lon))
    #     print(i.name)

def calculateDistance(node1, node2):
    graphNodeLat = node1['lat']
    graphNodeLong = node1['long']
    poiNodeLat = node2.lat
    poiNodeLong = node2.long
    return math.sqrt(math.pow(graphNodeLat - poiNodeLat, 2) + math.pow(graphNodeLong - poiNodeLong, 2))


def calculateLocationToPOIsDist(graph, startNode, pois, startToPOIs, sp, weight=1):
    for poi in pois:
        dist = ox.distance.great_circle_vec(startNode['lat'], startNode['long'], poi.lat, poi.long, earth_radius=6371009)
        dist *= weight 
        startToPOIs.append(dist)
        # print(dist)
    # print(len(startToPOIs))
    # print()


def createStart3DArray(startToPOIItemArray, startToPOIs, pois, travelWeight, costWeight):

    for i in range(len(startToPOIItemArray)):
        # now i am looking at each poi's item block
        dist = startToPOIs[i] # i have the dist from start to poi[i]
        
        for j in range(len(pois[i].items)):
            # now i am looking at each item in poi[i]
            if pois[i].items[j] == 0:
                continue
            else:
                cost = pois[i].items[j]
                startToPOIItemArray[i][j][0] = j
                startToPOIItemArray[i][j][1] = costWeight*cost + travelWeight*dist 

def create3DArray(currentToNextItemArray, sp, pois, travelWeight, costWeight):
    # currentToNextItem arr is currently all None

    for current in range(len(currentToNextItemArray)):

        for next in range(len(currentToNextItemArray)):
            dist = sp[current][next]
            for item in range(len(pois[next].items)):
                if pois[next].items[item] == 0:
                    continue
                else:
                    cost = pois[next].items[item]
                    currentToNextItemArray[current][next][item][0] = item
                    currentToNextItemArray[current][next][item][1] = costWeight*cost + travelWeight*dist


def readItems(items, location):
    with open(location, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()
            item = Item(line[0])
            item.POIs += line[1:]
            item.POIs = [int(poi) for poi in item.POIs]
            items.append(item)
    f.close()
    # print(items[0].POIs)
    # print(items[10].POIs)
    # print(items[20].POIs)
    # print(items[30].POIs)

def readStoreItems(pois, location):
    with open(location, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()
            poi = int(line[0])
            item = int(line[1])
            pois[poi].items[item] = float(line[2])

    f.close()

    # for i in pois[0].items:
    #     print(i)




main()

