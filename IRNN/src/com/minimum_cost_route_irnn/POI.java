package com.minimum_cost_route_irnn;

import org.graphstream.graph.Node;

public class POI {

    int ID;
    int inactive;
    double longitude;
    double latitude;
    String eId;
    double distance;
    double distance2;
    Node n1;
    Node n2;
    int  near;
    double[]items;
    void setvalues(int id, double lon, double lat, String e, double d )
    {
        ID= id;
        longitude = lon;
        latitude = lat;
        eId = e;
        distance = d;
        inactive=0;
    };

}
