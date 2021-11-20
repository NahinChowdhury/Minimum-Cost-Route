from Classes.Item import Item
from Classes.POI import POI
from Classes.Route import Route

import networkx as nx
import osmnx as ox
import json
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


    G = ox.graph_from_place("Edmonton, CA", network_type="drive")
    # fig, ax = ox.plot_graph(G)
    print(type(G))

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
    pois = []
    readPOIs(pois, "./Data/stores.txt")

    items = []
    readItems(items, "./Data/storesPerItem.txt")
    # print(len(pois[0].items))

    readStoreItems(pois, "./Data/itemsPerStore.txt")


    # create a 2d array calculating the distance between each pois
    # sp = [[0 for x in range(len(pois))] for y in range(len(pois))]
    # for i in range(len(pois)):
    #     for j in range(len(pois)):
    #         if i != j:
    #             sp[i][j] = calculateDistance(pois[i], pois[j])

    point = (pois[0].lat, pois[0].lon)

    u, v, k, edge_geom, dist = ox.distance.get_nearest_edge(G, point, return_geom=True, return_dist=True)

    # create shapely point geometry object as (x, y), that is (lng, lat)
    point_geom = Point(reversed(point))

    # use shapely to find the point along the edge that is closest to the reference point
    nearest_point_on_edge = edge_geom.interpolate(edge_geom.project(point_geom))
    nearest_point_on_edge.coords[0]
    # dist = ox.distance.shortest_path(G, pois[0].id, pois[1].id, weight='length', cpus=1)
    print(nearest_point_on_edge.coords[0])

def readPOIs(pois, location):
    # we should only take 51 stores because items are stored in 51 stores only
    with open(location, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.split()
            poi = POI(int(line[0]), float(line[1]), float(line[2]), " ".join(line[3:]))
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




main()

