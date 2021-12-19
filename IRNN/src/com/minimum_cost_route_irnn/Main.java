package com.minimum_cost_route_irnn;

import net.sf.geographiclib.Geodesic;
import net.sf.geographiclib.GeodesicData;
import org.graphstream.algorithm.Dijkstra;
import org.graphstream.graph.Edge;
import org.graphstream.graph.EdgeRejectedException;
import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.SingleGraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Random;

import static java.lang.Integer.parseInt;

public class Main {

    public static void readGraph(Graph graph){

        try (BufferedReader br = new BufferedReader(new FileReader("./datasets/Amsterdam/roadnetwork/RoadVerticesAMS.txt"))) {
            String line;

            while ((line = br.readLine()) != null) {

                String[] splited = line.split("\\s+"); 							// splits based on white spaces
                graph.addNode(splited[0]);										// add node
                Node n = graph.getNode(splited[0]);								// get the node as object
                n.setAttribute("long", Double.parseDouble(splited[1]));			// set longitude attribute
                n.setAttribute("lat", Double.parseDouble(splited[2]));			// set latitude attribute
                n.setAttribute("poi", 0);										// set the node as a normal node
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        try (BufferedReader br = new BufferedReader(new FileReader("./datasets/Amsterdam/roadnetwork/RoadEdgesAMS.txt"))) {
            String line;

            while ((line = br.readLine()) != null) {

                String[] splited = line.split("\\s+");							// split based on white spaces
                try
                {
                    graph.addEdge(splited[0], splited[1], splited[2]);    	// adds the edge

                    Node n1 = graph.getNode(splited[1]);					// gets the nodes from both ends of the edge
                    Node n2 = graph.getNode(splited[2]);

                    // gets the road network distance between the edge
                    GeodesicData g = Geodesic.WGS84.Inverse(n1.getAttribute("lat"), n1.getAttribute("long"), n2.getAttribute("lat"), n2.getAttribute("long"));
                    //System.out.println(g.s12);
                    Edge e1 = graph.getEdge(splited[0]);					// get the edge as object
                    e1.setAttribute("weight", g.s12);						// set the road network distance as the weight of the edge
                    //System.out.println(e1.getId());
                }
                catch (EdgeRejectedException e) {

                }

            }
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public static void readItems(ArrayList<Item> items, String location){
        //"./datasets/Amsterdam/poi/originals/StoresPerItemAMS.txt"

        try (BufferedReader br = new BufferedReader(new FileReader(location))) {

            String line;
            while ((line = br.readLine()) != null) {

                String[] splited = line.split("\\s+");							// split based on white spaces
                Item itm = new Item();
                itm.id = parseInt(splited[0]);
                //System.out.println(line);
                for(int i = 1; i<splited.length;i++)
                {
                    int shp = parseInt(splited[i]);
                    if(!itm.pois.contains(shp))
                    {
                        itm.pois.add(shp);
                    }
                }
                items.add(itm);
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static void readPoi(ArrayList<POI> P, int total, Graph graph, String location){
        // Reads the POIs  "./datasets/Amsterdam/poi/originals/GasStationsAMS.txt"

        try (BufferedReader br = new BufferedReader(new FileReader(location))) {

            String line;
            while ((line = br.readLine()) != null) {

                String[] splited = line.split("\\s+");							// split based on white spaces
                POI nPoi = new POI();
                nPoi.items = new double[total];
                nPoi.setvalues(parseInt(splited[0]), Double.parseDouble(splited[1]), Double.parseDouble(splited[2]), splited[3], Double.parseDouble(splited[4]));
                //System.out.println(splited[3]);

                Edge e = graph.getEdge(nPoi.eId);
                //System.out.println(splited[3]);
                //System.out.println(e.getId());
                //System.out.println(e.getNode1().getId());
                nPoi.n1 = e.getNode0();
                nPoi.n2 = e.getNode1();
                P.add(nPoi);
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static void readStoreItems(ArrayList<POI> pois, String location){
        try (BufferedReader br = new BufferedReader(new FileReader(location))) {

            String line;
            while ((line = br.readLine()) != null) {

                String[] splited = line.split("\\s+");							// split based on white spaces

                int poi = parseInt(splited[0]);
                int  itm = parseInt(splited[1]);
                //System.out.println(line);
                if(itm>1000)continue;
                //System.out.println(splited[2] + ": "+ splited[2].length());
                pois.get(poi).items[itm] = Double.parseDouble(splited[2]);

            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    public static void readSP(double [][] sp, int len, String location){
        int i=0;
//		int i=50;
        try (BufferedReader br = new BufferedReader(new FileReader(location))) {

            String line;
            while ((line = br.readLine()) != null) {
//				line = line.replaceAll("\r", "");
                String[] splited = line.split("\\s+"); // split based on white spaces
//				System.out.println(splited[49]);

                for (int j=0;j<len;j++)
                {
                    sp[i][j]= Double.parseDouble(splited[j]);
                }
                i++;

            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static void calculateLocationToPOIsDist(Node node, ArrayList<POI> pois, double [] locationToPOIs, Graph graph, Float weight){

        if(weight == null) weight = 1.0f;

        for (int i=0;i<pois.size();i++)
        {
            POI p = pois.get(i);
            Edge e = graph.getEdge(p.eId);
            p.n1 = e.getNode0();
            p.n2 = e.getNode1();
            Node n = p.n1;
            double s1 = n.getAttribute("long") ;
            s1 = s1 - p.longitude;
            s1 = s1*s1;

            double s2 = n.getAttribute("lat") ;
            s2 = s2 - p.latitude;
            s2 = s2*s2;
            double S1 = Math.sqrt(s1+s2);

            n = p.n2;
            s1 = n.getAttribute("long") ;
            s1 = s1 - p.longitude;
            s1 = s1*s1;

            s2 = n.getAttribute("lat") ;
            s2 = s2 - p.latitude;
            s2 = s2*s2;
            double S2 = Math.sqrt(s1+s2);
            if (S1<S2) p.near = 1;
            else p.near = 2;
            p.distance = S1;
            p.distance2 = S2;
//			p.distance2 = graph.getEdge(p.eId).getAttribute("weight");
//			p.distance2 = p.distance2 - p.distance;
        }

        Dijkstra dijkstra = new Dijkstra(Dijkstra.Element.EDGE ,"result", "weight");
        dijkstra.init(graph);
        dijkstra.setSource(node);
        dijkstra.compute();
        for(int j=0;j<pois.size();j++)
        {
            POI p1 = pois.get(j);
            //System.out.println(p1.n1.getId());
            double d1 = dijkstra.getPathLength(p1.n1)+p1.distance;
            double d2 = dijkstra.getPathLength(p1.n2)+p1.distance2;
//			if(p1.near==1)
//			{
//				d1  = d1 + p1.distance;
//				d2  = d2 + p1.distance2;
//			}
//			else {
//				d1  = d1 + p1.distance2;
//				d2  = d2 + p1.distance;
//			}
//			if (d1<d2) wareDist[j]=	Math.round((d1/maxd) * 100.0) / 100.0;
//			else wareDist[j]=Math.round((d2/maxd) * 100.0) / 100.0;
//			if (d1<d2) wareDist[j]=	precision(d1/maxd);
            if (d1<d2) locationToPOIs[j]=	precision(d1) * weight;
//			else wareDist[j]= precision(d2/maxd);
            else locationToPOIs[j]= precision(d2) * weight;

        }
    }

    public static void createStart3DArray(double [][][] startToPOIItemArray, double [] startToPOIs, ArrayList<POI> pois, ArrayList<Integer> itemsToBuy, Float travelWeight, Float costWeight){
        for (int i=0;i<startToPOIItemArray.length;i++)
        {
            double dist = startToPOIs[i];
            double itemCost = pois.get(i).items[itemsToBuy.get(0)];
//            for (int j=0;j<pois.get(i).items.length;j++)
//            {
//                if( !itemsToBuy.contains(j) || pois.get(i).items[j]== (double) 0) continue;
                if(itemCost <= 0) continue;
                else{
//                    double cost = pois.get(i).items[j];
                    startToPOIItemArray[i][0][0] = i;
                    startToPOIItemArray[i][0][1] = costWeight*itemCost + travelWeight*dist;
                }
//            }
        }
        return;
    }

    public static void create3DArray(double [][][][] currentToNextItemArray, double [][] sp, ArrayList<POI> pois, ArrayList<Integer> itemsToBuy, Float travelWeight, Float costWeight){
        for (int current=0; current <currentToNextItemArray.length;current++)
        {
            for (int next=0;next<currentToNextItemArray[current].length;next++)
            {
                double dist = sp[current][next];
                double itemCost = pois.get(next).items[itemsToBuy.get(0)];

//                for (int item=0;item<currentToNextItemArray[current][next].length;item++)
//                {
//                    if( !itemsToBuy.contains(item) || pois.get(next).items[item]== (double) 0) continue;
                if(itemCost <= 0) continue;
                else{
//                        double cost = pois.get(next).items[item];
                        currentToNextItemArray[current][next][0][0] = next;
                        currentToNextItemArray[current][next][0][1] = costWeight*itemCost + travelWeight*dist;
                    }
//                }
            }
        }
    }

    public static void createAll3DArray(double [][][][] currentToNextItemArray, double [][] sp, ArrayList<POI> pois, ArrayList<Integer> allItemsToBuy, Float travelWeight, Float costWeight){
        for (int current=0; current <currentToNextItemArray.length;current++)
        {
            for (int next=0;next<currentToNextItemArray[current].length;next++)
            {
                double dist = sp[current][next];
                for (int item=0;item<currentToNextItemArray[current][next].length;item++)
                {
                    if( !allItemsToBuy.contains(item) || pois.get(next).items[item]== (double) 0) continue;
                    else{
                        double cost = pois.get(next).items[item];
                        currentToNextItemArray[current][next][item][0] = item;
                        currentToNextItemArray[current][next][item][1] = costWeight*cost + travelWeight*dist;
                    }
                }
            }
        }
    }

    public static ArrayList<Integer> getAvailablePOIs(ArrayList<POI>pois, ArrayList<Integer>itemsToBuy){

        ArrayList<Integer> availablePOIs = new ArrayList<Integer>();
        for (int i=0;i<pois.size();i++)
        {
            if(pois.get(i).items[itemsToBuy.get(0)]>0)
                availablePOIs.add(i);
        }
        return availablePOIs;
    }

    public static Route computeNewRoute(Route route, double [][][][] allCurrentToNextItemArray, double [][] sortedStartToPOIItemArray, double [][][] sortedCurrentToNextItemArray, double [] endToPOIs, ArrayList<Integer> availablePOIs, ArrayList<Integer> itemsToBuy){

        Route finalRoute = new Route(route);
        // itemToBuy is given
        // availablePOIs is given

        BigDecimal minCost = BigDecimal.valueOf(Double.MAX_VALUE);//max value acts as infinity
        BigDecimal minCost2 = BigDecimal.valueOf(Double.MAX_VALUE);//max value acts as infinity
        Integer candidateNN = null;
        Integer prevPOI = null;
        Integer nextPOI = null;

        // should update candidateNN and prevPOI and nextPOI by considering going from start location to availblePOIs to first poi
        // then run the loop for all pois in the route

        for(int j = 0; j < availablePOIs.size(); j++){

            double firstCost = sortedStartToPOIItemArray[availablePOIs.get(j)][1];
            int itemToBuy = route.items.get(0);
            double secondCost = allCurrentToNextItemArray[availablePOIs.get(j)][route.pois.get(0)][itemToBuy][1];
            double combinedNewCost = firstCost + secondCost;

            // if item cannot be bought from any of the stores, then the combined cost is infinity
            if( ( minCost.floatValue() + minCost2.floatValue() ) > combinedNewCost){
                minCost = BigDecimal.valueOf(firstCost);
                minCost2 = BigDecimal.valueOf(secondCost);
                candidateNN = availablePOIs.get(j);
                prevPOI = -1;
                nextPOI = 0;
            }
        }

        for(int i = 0; i < route.pois.size(); i++){
            double potentialMinCost = Double.POSITIVE_INFINITY;
            double potentialMinCost2 = Double.POSITIVE_INFINITY;
            Integer potentialCandidateNN = null;
            Integer potentialPrevPOI = null;
            Integer potentialNextPOI = null;



            if(i == route.pois.size() - 1) {

                for (int j = 0; j < availablePOIs.size(); j++) {
                    double firstCost = allCurrentToNextItemArray[i][availablePOIs.get(j)][itemsToBuy.get(0)][1];
                    double secondCost = endToPOIs[availablePOIs.get(j)];
                    double combinedNewCost = firstCost + secondCost;

                    if (( potentialMinCost + potentialMinCost2 ) > combinedNewCost) {
                        potentialMinCost = firstCost;
                        potentialMinCost2 = secondCost;
                        potentialCandidateNN = availablePOIs.get(j);
                        potentialPrevPOI = i;
                        potentialNextPOI = -1; // -1 represents end location
                    }

                }
            }else if(route.pois.get(i) == route.pois.get(i+1)){

                continue;

            }else{

                for(int j = 0; j < availablePOIs.size(); j++){

                    double firstCost = sortedCurrentToNextItemArray[route.pois.get(i)][availablePOIs.get(j)][1];
                    int itemToBuy = route.items.get(i+1);
                    double secondCost = allCurrentToNextItemArray[availablePOIs.get(j)][route.pois.get(i+1)][itemToBuy][1];
                    double combinedNewCost = firstCost + secondCost;

                    // if item cannot be bought from any of the stores, then the combined cost is infinity
                    if( ( potentialMinCost + potentialMinCost2 ) > combinedNewCost){
                        potentialMinCost = firstCost;
                        potentialMinCost2 = secondCost;
                        potentialCandidateNN = availablePOIs.get(j);
                        potentialPrevPOI = i;
                        potentialNextPOI = i+1;
                    }
                }
            }

            if( ( minCost.floatValue() + minCost2.floatValue() ) > potentialMinCost){
                minCost = BigDecimal.valueOf(potentialMinCost);
                minCost2 = BigDecimal.valueOf(potentialMinCost2);
                candidateNN = potentialCandidateNN;
                prevPOI = potentialPrevPOI;
                nextPOI = potentialNextPOI;
            }

        }

        if(candidateNN != null){
            int changed = 0;
            for(int i = 0; i < route.pois.size()-1; i++){
                if(prevPOI == i && nextPOI == i+1){
                    finalRoute.totalCost -= finalRoute.costs.get(i);
                    finalRoute.pois.add(i+1, candidateNN);
                    finalRoute.costs.add(i+1, minCost.floatValue() );
                    finalRoute.costs.set(i+2, minCost2.floatValue());

                    finalRoute.items.add(i+1, itemsToBuy.get(0));
                    finalRoute.totalCost += minCost.floatValue() + minCost2.floatValue();
                    changed = 1;
                    System.out.println("----------------------------------------------------------------------------------------------------------------------");
                    System.out.println("----------------------------------------------------------------------------------------------------------------------");
                    break;
                }
            }
            if(changed == 0 && prevPOI == route.pois.size()-1 && nextPOI == -1){
                finalRoute.totalCost -= finalRoute.costs.get( finalRoute.costs.size() - 1 );
                finalRoute.pois.add(candidateNN);
                finalRoute.costs.set( finalRoute.costs.size() - 1 , minCost.floatValue() );
                finalRoute.items.add(itemsToBuy.get(0));

                finalRoute.totalCost += minCost.floatValue() + minCost2.floatValue();

            }else if(changed == 0 && prevPOI == -1 && nextPOI == 0){
                finalRoute.totalCost -= finalRoute.costs.get(0);
                finalRoute.pois.add(0, candidateNN);
                finalRoute.costs.add(0, minCost.floatValue() );
                finalRoute.costs.set(1, minCost2.floatValue() );
                finalRoute.items.add(0, itemsToBuy.get(0));

                finalRoute.totalCost += minCost.floatValue() + minCost2.floatValue();
            }
        }else{
            System.out.println("No candidate found");
        }

        return finalRoute;
    }

    public static Comparator<Route> idComparator = new Comparator<Route>(){

        @Override
        public int compare(Route c1, Route c2) {
            return (int) (c1.totalCost - c2.totalCost);
        }
    };

    public static double precision(double val) {
        DecimalFormat df = new DecimalFormat("#.000");
        return Double.valueOf(df.format(val));
    }



    public static void main(String[] args) {
	// write your code here

        Float travelWeight = 1.0f; // the weight that the user gives to travel
        Float costWeight = 1 - travelWeight; // the weight that the user gives to cost

        ArrayList<Integer> itemsToBuy = new ArrayList<Integer>();

//        for (int i = 0;i<1;i++){
        Random rand = new Random(); //instance of random class
        int upperBound = 999;
        itemsToBuy.add(rand.nextInt(upperBound)); // creating a shopping list
//        }
        // pois: [43, 43, 38, 38, 27]
        int[] givenPois = new int[]{43, 43, 38, 38, 27};
        ArrayList<Integer> allItemsToBuy = new ArrayList<Integer>();
        allItemsToBuy.add(0);
        allItemsToBuy.add(1);
        allItemsToBuy.add(2);
        allItemsToBuy.add(3);
        allItemsToBuy.add(4);

        Float[] givenCosts = new Float[]{150.0f, 10.0f, 15300.0f, 0.0f, 1000.0f};




        // items: [0, 4, 2, 3, 1]
        Route route = new Route();
        for (int i = 0;i<givenPois.length;i++){
            route.pois.add(givenPois[i]);
        }
        for (int i = 0;i<allItemsToBuy.size();i++){
            route.items.add(allItemsToBuy.get(i));
        }
        for (int i = 0;i<givenCosts.length;i++){
            route.costs.add(givenCosts[i]);
            route.totalCost += givenCosts[i];
        }
        route.latestPoi = route.pois.get(route.pois.size()-1);
        route.latestItem = route.items.get(route.items.size()-1);
        route.latestCost = givenCosts[givenCosts.length-1];

        long preComputationStartTime = System.currentTimeMillis();
        ArrayList<POI> pois = new ArrayList<POI>();
        ArrayList<Item> items = new ArrayList<Item>();

        Graph graph = new SingleGraph("Amsterdam");
        readGraph(graph); // reading graph from available data

        readItems(items,"./datasets/Amsterdam/poi/originals/Item_Cost/StoresPerItemAMS_Cost_INward.txt");

        readPoi(pois, items.size(), graph, "./datasets/Amsterdam/poi/originals/PoiAMS50.txt");										// Reads all the POI information from input file

        readStoreItems(pois, "./datasets/Amsterdam/poi/originals/Item_Cost/ItemsPerStoreAMS_Cost_INward.txt");

//        pickPoi(pois, 20); // pickPois should set poi to inactive and then make sure sp of other pois is set to infinity

//        setnewItem(pois, items);

        double [][] sp = new double [pois.size()][pois.size()] ;
        readSP(sp, pois.size(), "./datasets/Amsterdam/poi/originals/ShortestPathPoi50.txt");


        Node startNode = graph.getNode(8020);
        double [] startToPOIs = new double [pois.size()];
        calculateLocationToPOIsDist(startNode, pois, startToPOIs, graph, 1.0f); // calculate the distance from start node to all POIs. Weight is zero here because we use multiple the weight with the distance while computing routes

        Node endNode = graph.getNode(88896);
        double [] endToPOIs = new double [pois.size()];
        calculateLocationToPOIsDist(endNode, pois, endToPOIs, graph, travelWeight); // calculate the distance from end node to all POIs. Weight is travelWeight here because we don't compute the new travel distance while computing the route

        double [][][] startToPOIItemArray = new double [pois.size()][1][2]; // start to [poi][item][0] is the itemID, start to [poi][item][1] is the combined weight of the cost of item and travel distance

        // set startToPOIItemArray elements to infinity. Can't use null because of multiple type cast conflict. Infinity works when we want to sort by lowest cost
        for (int i=0;i<pois.size();i++){
//            for (int j=0;j<items.size();j++){
                startToPOIItemArray[i][0][0] = Double.POSITIVE_INFINITY;
                startToPOIItemArray[i][0][1] = Double.POSITIVE_INFINITY;
//            }
        }
//        for (int iter = 0; iter < 100; iter++) {

//            itemsToBuy.set(0, rand.nextInt(upperBound));
        createStart3DArray(startToPOIItemArray, startToPOIs, pois, itemsToBuy, travelWeight, costWeight); // create 3D array containing the combined weight of cost and travel dist from start loc to all pois and items


        double [][][][] currentToNextItemArray = new double [pois.size()][pois.size()][1][2]; // current poi to [next poi][item][0] is the itemID, current poi to [next poi][item][1] is the combined weight of the cost of item and travel distance
        double [][][][] allCurrentToNextItemArray = new double [pois.size()][pois.size()][items.size()][2]; // current poi to [next poi][item][0] is the itemID, current poi to [next poi][item][1] is the combined weight of the cost of item and travel distance

        // set currentToNextItemArray elements to infinity. Can't use null because of multiple type cast conflict. Infinity works when we want to sort by lowest cost
        for (int i=0;i<pois.size();i++){
            for (int j=0;j<pois.size();j++){
//                for (int k=0;k<items.size();k++){
                    currentToNextItemArray[i][j][0][0] = Double.POSITIVE_INFINITY;
                    currentToNextItemArray[i][j][0][1] = Double.POSITIVE_INFINITY;
//                }
            }

        }

        for (int i=0;i<pois.size();i++){
            for (int j=0;j<pois.size();j++){
                for (int k=0;k<items.size();k++){
                    allCurrentToNextItemArray[i][j][k][0] = Double.POSITIVE_INFINITY;
                    allCurrentToNextItemArray[i][j][k][1] = Double.POSITIVE_INFINITY;
                }
            }

        }

        create3DArray(currentToNextItemArray, sp, pois, itemsToBuy, travelWeight, costWeight); // techincally creating a 4D array from one POi to another POi to buy item X. [current][next][item] has 2 values: [0] is the itemID, [1] is the combined weight of the cost of item and travel distance
        createAll3DArray(allCurrentToNextItemArray, sp, pois, allItemsToBuy, travelWeight, costWeight); // techincally creating a 4D array from one POi to another POi to buy item X. [current][next][item] has 2 values: [0] is the itemID, [1] is the combined weight of the cost of item and travel distance


        double [][] sortedStartToPOIItemArray = new double [startToPOIItemArray.length][2];

        for(int i=0;i<startToPOIItemArray.length;i++){

            sortedStartToPOIItemArray[i][0] = i;
            sortedStartToPOIItemArray[i][1] = startToPOIItemArray[i][0][1];
        }
//        System.out.println(1);
        java.util.Arrays.sort(sortedStartToPOIItemArray, new java.util.Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                return Double.compare(a[1], b[1]);
            }
        });


//        for(int i = 0; i < startToPOIItemArray.length; i++){
//            java.util.Arrays.sort(startToPOIItemArray[i], new java.util.Comparator<double[]>() {
//                public int compare(double[] a, double[] b) {
//                    return Double.compare(a[1], b[1]);
//                }
//            });
//        }

        // sorting the currentToNextItemArray by weight
        for(int i=0;i<currentToNextItemArray.length;i++){
            for (int j=0;j<currentToNextItemArray[i].length;j++) {
                java.util.Arrays.sort(currentToNextItemArray[i][j], new java.util.Comparator<double[]>() {
                    public int compare(double[] a, double[] b) {
                        return Double.compare(a[1], b[1]);
                    }
                });
            }
        }

        double [][][] sortedCurrentToNextItemArray = new double [currentToNextItemArray.length][currentToNextItemArray.length][2]; // we all know that only one item is being bought. so we go from one poi to another and get the cost of buying the item to ge bought.
        // if we have multiple items to buy, then we should change this implementation

        for(int i=0;i<currentToNextItemArray.length;i++){
            for (int j=0;j<currentToNextItemArray[i].length;j++) {
                sortedCurrentToNextItemArray[i][j][0] = j;
                sortedCurrentToNextItemArray[i][j][1] = currentToNextItemArray[i][j][0][1];
            }
        }

        for(int i=0;i<sortedCurrentToNextItemArray.length;i++){
            java.util.Arrays.sort(sortedCurrentToNextItemArray[i], new java.util.Comparator<double[]>() {
                public int compare(double[] a, double[] b) {
                        return Double.compare(a[1], b[1]);
                    }
            });
        }

        for(int i=0;i<allCurrentToNextItemArray.length;i++){
            for (int j=0;j<allCurrentToNextItemArray[i].length;j++) {
                java.util.Arrays.sort(allCurrentToNextItemArray[i][j], new java.util.Comparator<double[]>() {
                    public int compare(double[] a, double[] b) {
                        return Double.compare(a[1], b[1]);
                    }
                });
            }
        }

        ArrayList<Integer> availablePOIs = new ArrayList<>();
        availablePOIs = getAvailablePOIs(pois, itemsToBuy);


        long preComputationEndTime = System.currentTimeMillis();
        double preComputationDuration = (double) (preComputationEndTime - preComputationStartTime)/ (double) 1000;

        // before the for loop, the upperbound is used as the start to some poi smallest weight

        // go over every single poi
        // in the 3d array, for each POI, find the next POI with the lowest combined weight of cost and travel distance
        // The poi with the lowest

//        double [] upperBound = new double[3];
//        upperBound[0] = Float.POSITIVE_INFINITY; // last poi from the sorted array
//        upperBound[1] = sortedStartToPOIItemArray[0][0]; // next poi number
//        upperBound[2] = sortedStartToPOIItemArray[0][1]; // cost
//
//        for(int i = 0; i < route.pois.size(); i++){
//            int poi = route.pois.get(i);
//            if(sortedCurrentToNextItemArray[poi][0][0] != Float.POSITIVE_INFINITY &&
//                    sortedCurrentToNextItemArray[poi][0][1] != Float.POSITIVE_INFINITY &&
//                    sortedCurrentToNextItemArray[poi][0][1] < upperBound[2]){
//                upperBound[0] = poi;
//                upperBound[1] = sortedCurrentToNextItemArray[poi][0][0];
//                upperBound[2] = sortedCurrentToNextItemArray[poi][0][1];
//            }
//
//        }
//
//        int lastPOIindex = route.pois.lastIndexOf( (int) upperBound[0]);
//        route.pois.add(lastPOIindex, (int) upperBound[1]);
//        route.items.add(lastPOIindex, itemsToBuy.get(0)); // need to make sure the original route has multiple duplicate pois for different items bought


            Route finalRoute = new Route();
            finalRoute = computeNewRoute(route, allCurrentToNextItemArray, sortedStartToPOIItemArray, sortedCurrentToNextItemArray, endToPOIs, availablePOIs, itemsToBuy);

            System.out.println(route.toString());
            System.out.println(finalRoute.toString());
            System.out.println("\n");
//        }


    }
}
