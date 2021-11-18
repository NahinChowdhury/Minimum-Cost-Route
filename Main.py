from Classes.Item import Item
from Classes.POI import POI
from Classes.Route import Route

import networkx as nx
import osmnx as ox

"""
Boeing, G. 2017. OSMnx: New Methods for Acquiring, Constructing, Analyzing, and Visualizing Complex Street Networks. Computers, Environment and Urban Systems 65, 126-139. doi:10.1016/j.compenvurbsys.2017.05.004
"""
def main():
    # Create a route
    route = Route()
    # Create a POI
    poi = POI(1, 1.11, 3.22)
    # Create an item
    item = Item(1, "Apple")
    # Add the POI to the route
    route.addPOI(poi)
    # Add the item to the POI
    poi.setUseless()

    route.POIs[0].id = 15
    # Print the route
    print(route.POIs[0].id)


    G = ox.graph_from_place("Edmonton, CA", network_type="drive")
    fig, ax = ox.plot_graph(G)
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


    # f.close()


main()