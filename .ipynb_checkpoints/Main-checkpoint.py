from Classes.Item import Item
from Classes.POI import POI
from Classes.Route import Route
from DataStructures.PriorityQueue import PriorityQueue

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
    # # Print the routes
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
    preCompiutationStartTime = time.time()
    travelWeight = 1
    costWeight = 1 - travelWeight

    itemsToBuy = [1,2]#,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

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
    calculateLocationToPOIsDist(startNode, pois, startToPOIs, 1)
    
    endNode = graph.nodes[88896]
    # calculate distance between each POI and endNode
    endToPOIs = []
    calculateLocationToPOIsDist(endNode, pois, endToPOIs, travelWeight)

    # have a 3d array for start to each POI and buying items
    startToPOIItemArray = [[[None, None] for i in range(1000)] for i in range(len(pois))]
    createStart3DArray(startToPOIItemArray, startToPOIs, pois, itemsToBuy, travelWeight, costWeight) # indices with value null means that the item cannot be bought from that POI
    # startToPOIItemArray0_sorted = (sorted(startToPOIItemArray[0], key=lambda x: float('inf') if x[1] is None else x[1])) # this is how we sort rows in a 2d array
    # need to pass a single row to sort the row above
    # print(startToPOIItemArray[1])

    # maybe I can sort every single row of items so that our algorithm can find the cheapest item to buy faster
    

    # sort startToPOIItemArray by the cost of the item to be bought
    # startToPOIItemArray = sorted(startToPOIItemArray, key=lambda x: x[1])
    # print(startToPOIItemArray0_sorted)
    # print(startToPOIItemArray[1])
    # print(startToPOIItemArray[1])

    # create a 3D array of 51 POIs with 51 POIs with 999 items with default value of None
    currentToNextItemArray = [[[[None, None] for j in range(1000)] for i in range(len(pois))] for k in range(len(pois))]
    create3DArray(currentToNextItemArray, sp, pois, itemsToBuy, travelWeight, costWeight)
    # print(len(currentToNextItemArray))
    # print(currentToNextItemArray[0])
    # print(len(currentToNextItemArray[0][0]))
    # currentToNextItemArray0_12_sorted = (sorted(currentToNextItemArray[0][12], key=lambda x: float('inf') if x[1] is None else x[1])) # this is how we sort rows in a 2d array
    # print(currentToNextItemArray0_12_sorted)

    preCompiutationEndTime = time.time()
    preCompiutationDuration = preCompiutationEndTime - preCompiutationStartTime
    print("preCompiutationDuration: {}".format(preCompiutationDuration))
    # our last POIS to ending location is just travelWeight times the travel distance
    # i do that by passing the travel weight while calculating the distance

    # now we can implement the algorithm

    firstAlgStartTime = time.time()
    route = findBestRoute(startToPOIItemArray.copy(), currentToNextItemArray.copy(), endToPOIs.copy(), itemsToBuy.copy())
    firstAlgEndTime = time.time()
    firstAlgDuration = firstAlgEndTime - firstAlgStartTime
    print(route)
    print("firstAlgDuration: {}".format(firstAlgDuration))


    print()

    secondAlgStartTime = time.time()
    route = findBestRouteContinued(startToPOIItemArray.copy(), currentToNextItemArray.copy(), endToPOIs.copy(), itemsToBuy.copy())
    secondAlgEndTime = time.time()
    secondAlgDuration = secondAlgEndTime - secondAlgStartTime
    print(route)
    print("secondAlgDuration: {}".format(secondAlgDuration))





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


def calculateDistance(node1, node2):
    graphNodeLat = node1['lat']
    graphNodeLong = node1['long']
    poiNodeLat = node2.lat
    poiNodeLong = node2.long
    return math.sqrt(math.pow(graphNodeLat - poiNodeLat, 2) + math.pow(graphNodeLong - poiNodeLong, 2))

def calculateLocationToPOIsDist(startNode, pois, startToPOIs, weight=1):
    for poi in pois:
        dist = ox.distance.great_circle_vec(startNode['lat'], startNode['long'], poi.lat, poi.long, earth_radius=6371009)
        dist *= weight 
        startToPOIs.append(dist)
        # print(dist)
    # print(len(startToPOIs))
    # print()

def createStart3DArray(startToPOIItemArray, startToPOIs, pois, itemsToBuy, travelWeight, costWeight):

    for i in range(len(startToPOIItemArray)):
        # now i am looking at each poi's item block
        dist = startToPOIs[i] # i have the dist from start to poi[i]
        
        for j in range(len(pois[i].items)):
            # now i am looking at each item in poi[i]
            if j not in itemsToBuy or pois[i].items[j] == 0:
                continue
            else:
                cost = pois[i].items[j]
                startToPOIItemArray[i][j][0] = j
                startToPOIItemArray[i][j][1] = costWeight*cost + travelWeight*dist 

def create3DArray(currentToNextItemArray, sp, pois, itemsToBuy, travelWeight, costWeight):
    # currentToNextItem arr is currently all None

    for current in range(len(currentToNextItemArray)):

        for next in range(len(currentToNextItemArray)):
            dist = sp[current][next]
            for item in range(len(pois[next].items)):
                if item not in itemsToBuy or pois[next].items[item] == 0:
                    continue
                else:
                    cost = pois[next].items[item]
                    currentToNextItemArray[current][next][item][0] = item
                    currentToNextItemArray[current][next][item][1] = costWeight*cost + travelWeight*dist


def findBestRoute(startToPOIItemArray, currentToNextItemArray, endToPOIs, itemsToBuy):
    potentialRoute = None
    PQ = PriorityQueue()

    startRoute = Route()
    startRouteInfo = MinCost(startRoute, startToPOIItemArray.copy())
    # print(startRouteInfo)

    startRoute.addPOI(startRouteInfo[0])
    startRoute.addItem(startRouteInfo[1])
    startRoute.addCost(startRouteInfo[2]) # no need to error check here because we know it is not None

    # print(startRoute)

    PQ.insert(startRoute)

    # startRoute1 = Route()
    # startRoute1.addPOI(1)
    # startRoute1.addItem(1)
    # startRoute1.addCost(1)

    # PQ.insert(startRoute1)
    # print("poi len: {}".format(len(startRoute.POIs)))

    count = 0
    print(count)
    while not PQ.isEmpty():
        r = PQ.pop()
        
        if SatisfyList(r, itemsToBuy):
            # print("satisfied")
            # get the latestPOI
            # for debugging
            if len(r.items) != len(itemsToBuy):
                print("error")
                for i in range(len(r.items)):
                    print("poi: {}, item: {}, weight: {}".format(r.pois[i], r.items[i], r.costs[i]))
                
                return

            latestPOI = r.latestPOI
            finalWeight = endToPOIs[latestPOI]
            r.addFinalCost(finalWeight)
            potentialRoute = r
            return potentialRoute

        nextR = nextItemRoute(r, currentToNextItemArray.copy())

        # special case for then we want to get the next best POI from starting point
        nextMinCostR = duplicateRoute(r)
        # print("poi len: {}".format(len(r.POIs)))
        if len(r.POIs) == 1:
            # print("in")
            nextMinCostR.removeLatestPOI()
            # nextMinCostR.updateN(0)
            # print(nextMinCostR.n)
            n = r.n[0] + 1

            startRouteInfo = MinCost(nextMinCostR, startToPOIItemArray.copy(), n)

            if startRouteInfo == None or None in startRouteInfo:
                continue
            # print(startRouteInfo)
            # print(startRouteInfo)
            
            nextMinCostR.addPOI(startRouteInfo[0])
            nextMinCostR.n[0] = n
            nextMinCostR.addItem(startRouteInfo[1])
            nextMinCostR.addCost(startRouteInfo[2])
            # print(nextMinCostR.POIs)

        else:
            nextMinCostR = nextMinCostRoute(r, currentToNextItemArray.copy())

        if nextR != None:
            # print("nextR below")
            # print(nextR)
            PQ.insert(nextR)
        if nextMinCostR != None:
            # print("nextMinCost below")
            # print(nextMinCostR)
            PQ.insert(nextMinCostR)
        
        
        count += 1
        # print(count)

    print("count for first alg: {}".format(count))
    return potentialRoute



def findBestRouteContinued(startToPOIItemArray, currentToNextItemArray, endToPOIs, itemsToBuy):
    potentialRoute = None
    PQ = PriorityQueue()

    startRoute = Route()
    startRouteInfo = MinCost(startRoute, startToPOIItemArray.copy())
    # print(startRouteInfo)

    startRoute.addPOI(startRouteInfo[0])
    startRoute.addItem(startRouteInfo[1])
    startRoute.addCost(startRouteInfo[2])

    # print(startRoute)

    PQ.insert(startRoute)

    weightUpperBound = float('inf')

    # startRoute1 = Route()
    # startRoute1.addPOI(1)
    # startRoute1.addItem(1)
    # startRoute1.addCost(1)

    # PQ.insert(startRoute1)
    # print("poi len: {}".format(len(startRoute.POIs)))

    count = 0
    pruned = 0
    while not PQ.isEmpty():
        r = PQ.pop()
        print(r.POIs)
        if r.totalCost > weightUpperBound:
            pruned += 1
            continue
        
        if SatisfyList(r, itemsToBuy):
            # print("satisfied")
            # get the latestPOI
            # for debugging
            if len(r.items) != len(itemsToBuy):
                print("error")
                for i in range(len(r.items)):
                    print("poi: {}, item: {}, weight: {}".format(r.pois[i], r.items[i], r.costs[i]))
                
                return

            latestPOI = r.latestPOI
            finalWeight = endToPOIs[latestPOI]
            r.addFinalCost(finalWeight)

            if potentialRoute == None:
                potentialRoute = r
                weightUpperBound = r.totalCost
            else:
                if r.totalCost < potentialRoute.totalCost:
                    print("accepted below")
                    potentialRoute = r
                    weightUpperBound = r.totalCost
            print("Upperbound cost: " + str(potentialRoute.totalCost) + ", potential cost: {}".format(r.totalCost) )
            continue
            # return potentialRoute

        nextR = nextItemRoute(r, currentToNextItemArray.copy())

        nextMinCostR = duplicateRoute(r)
        
        # print("poi len: {}".format(len(r.POIs)))
        if len(r.POIs) == 1:
            # print("in")
            nextMinCostR.removeLatestPOI()

            # nextMinCostR.updateN(0)
            n = r.n[0] + 1
            startRouteInfo = MinCost(nextMinCostR, startToPOIItemArray.copy(), n)
            # print(startRouteInfo)
            if startRouteInfo == None or None in startRouteInfo:
                continue
            
            nextMinCostR.addPOI(startRouteInfo[0])
            nextMinCostR.n[0] = n
            nextMinCostR.addItem(startRouteInfo[1])
            nextMinCostR.addCost(startRouteInfo[2])

        else:
            nextMinCostR = nextMinCostRoute(r, currentToNextItemArray.copy())

        if nextR != None and nextR.totalCost < weightUpperBound:
            # print("nextR below")
            # print(nextR)
            PQ.insert(nextR)
        if nextMinCostR != None and nextMinCostR.totalCost < weightUpperBound:
            # print("nextMinCost below")
            # print(nextMinCostR)
            PQ.insert(nextMinCostR)
        
        # print(count)
        count += 1

    print("count for second alg: {}".format(count))
    print("pruned: {}".format(pruned))
    return potentialRoute














def SatisfyList(route, itemsToBuy):
    # if there is one item missing, return false
    for i in itemsToBuy:
        if i not in route.items:
            return False
    return True



def nextItemRoute(route, currentToNextItemArray):
    newRoute = duplicateRoute(route)
    
    currPOI = route.latestPOI
    nextRouteInfo = getMinCost(newRoute, currentToNextItemArray[currPOI].copy())

    if nextRouteInfo == None or None in nextRouteInfo:
        return None

    newRoute.addPOI(nextRouteInfo[0])
    newRoute.addItem(nextRouteInfo[1])
    newRoute.addCost(nextRouteInfo[2])

    # newRoute.updateN(len(newRoute.POIs) - 1)

    return newRoute


def nextMinCostRoute(route, currentToNextItemArray):
    

    newRoute = duplicateRoute(route)
    
    newRoute.removeLatestPOI()
    n = route.n[len(newRoute.POIs) - 1] + 1

    newRoute.n[len(newRoute.POIs) - 1] = n

    currPOI = newRoute.latestPOI
    nextRouteInfo = getMinCost(newRoute, currentToNextItemArray[currPOI].copy())
    if nextRouteInfo == None or None in nextRouteInfo:
        return None

    newRoute.addPOI(nextRouteInfo[0])
    newRoute.addItem(nextRouteInfo[1])
    newRoute.addCost(nextRouteInfo[2])

    

    return newRoute



def MinCost(route, startToPOIItemArray, n = 0):

    if len(route.n) == 0:
        pass
    else:
        n = route.n[len(route.POIs) - 1]

    for nVal in route.n:
        if nVal >= len(startToPOIItemArray):
            return None

    # get the sorted weights for each row(or POI)
    # create a new array with poi number, item number, and weight
    # sort the array by weight
    # get the first item number and poi number and weight

    sortedItemsPerPOI = []
    for i in range(len(startToPOIItemArray)):
        sortedItemsPerPOI.append(sorted(startToPOIItemArray[i], key=lambda x: float('inf') if x[1] is None else x[1]))

    # print(startToPOIItemArray[0])

    # for i in sortedItemsPerPOI:
        # print(i[0])

    # print("\n")
    sortedPOIs = []
    for i in range(len(sortedItemsPerPOI)):
        poiItemWeightList = []
        poiItemWeightList.append(i) 
        poiItemWeightList.append(sortedItemsPerPOI[i][0][0]) 
        poiItemWeightList.append(sortedItemsPerPOI[i][0][1]) # i 0 1 here because i want the best weight from each poi
        sortedPOIs.append(poiItemWeightList)

    sortedPOIs = sorted(sortedPOIs, key=lambda x: float('inf') if x[2] is None else x[2])
    # print(sortedPOIs)

    if n == 0:
        # print(sortedPOIs)
        pass
    return sortedPOIs[n]

def getMinCost(route, array):

    n = route.n[- 1]
    
    for nVal in route.n:
        if nVal >= len(array):
            return None

    sortedItemsPerPOI = []
    for i in range(len(array)):
        sortedItemsPerPOI.append(sorted(array[i], key=lambda x: float('inf') if x[1] is None else x[1]))


    # I have a list of sorted items
    # i should filter each row such that eahc row only has the items i need to buy
    # then i should sort the rows by weight
    # then i should get the first item and poi number and weight

    sortedPOIs = []
    for i in range(len(sortedItemsPerPOI)):
        poiItemWeightList = []
        poiItemWeightList.append(i)

        sortedItemsPerPOI = filterArray(sortedItemsPerPOI.copy(), route, i) # filters every
        poiItemWeightList.append(sortedItemsPerPOI[i][0][0]) # new filtered array's index 0 should have the lowest weight item which i havent bough yet
        poiItemWeightList.append(sortedItemsPerPOI[i][0][1])
        sortedPOIs.append(poiItemWeightList)

    
    sortedPOIs = sorted(sortedPOIs, key=lambda x: float('inf') if x[2] is None else x[2])
    # print(sortedPOIs)

    return sortedPOIs[n]



def filterArray(sortedItemsPerPOI, route, i):
    
    j = 0
    length = len(sortedItemsPerPOI[i])
    while j <= length - 1:
        if sortedItemsPerPOI[i][j][1] == None: # if the item is not bought from the poi as it has cost of None
            return sortedItemsPerPOI
        if sortedItemsPerPOI[i][j][0] in route.items:
            sortedItemsPerPOI[i].pop(j)
            length = len(sortedItemsPerPOI[i])
        else:
            j += 1
    
    return sortedItemsPerPOI


def duplicateRoute(route):
    newRoute = Route()
    newRoute.POIs = [] 
    for poi in route.POIs:
        newRoute.addPOI(poi)

    newRoute.items = []
    for item in route.items:
        newRoute.addItem(item)

    newRoute.costs = []
    for cost in route.costs:
        newRoute.addCost(cost)

    newRoute.totalCost = route.totalCost

    i = 0
    for n in route.n:
        newRoute.n[i] = n


    return newRoute


main()

