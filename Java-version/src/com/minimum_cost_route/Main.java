package com.minimum_cost_route;

import net.sf.geographiclib.Geodesic;
import net.sf.geographiclib.GeodesicData;
import org.graphstream.algorithm.Dijkstra;
import org.graphstream.graph.*;
import org.graphstream.graph.implementations.SingleGraph;

import java.io.*;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.*;

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

    public static void pickPoi(ArrayList<POI> pois, int N)
    {
//		int count =0;
        ArrayList<POI> nPois = new ArrayList<POI>();
//		while(count<N)
        // N is 25
        for(int i = 0; i < 51; i++)
        {
			Random rand = new Random();
			int rand_int1 = rand.nextInt(pois.size());
			POI poi = pois.get(rand_int1);
//            POI poi = pois.get(i);
//			System.out.println(poi.ID + " " + poi.latitude + " " + poi.longitude);

            if(i % 2 == 1){

                if(!nPois.contains(poi))
                {
                    nPois.add(poi);
                    //System.out.println("Picked: "+ poi.ID);
//				count++;
                }
            }

        }

//        for(int i=0;i<pois.size();i++)
//        {
//            POI p = pois.get(i);
//            if(!nPois.contains(p))
//            {
////                Arrays.fill(p.items,0.0);
////                p.inactive=1;
//
//            }
//
//        }
        int length = pois.size();
        int i = 0;
        while(i < length){
            POI p = pois.get(i);
            if(!nPois.contains(p)){
                pois.remove(p);
            }else{
                i++;
            }
        }
        return;

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
            for (int j=0;j<pois.get(i).items.length;j++)
            {
                if( !itemsToBuy.contains(j) || pois.get(i).items[j]== (double) 0) continue;
                else{
                    double cost = pois.get(i).items[j];
                    startToPOIItemArray[i][j][0] = j;
                    startToPOIItemArray[i][j][1] = costWeight*cost + travelWeight*dist;
                }
            }
        }
        return;
    }

    public static void create3DArray(double [][][][] currentToNextItemArray, double [][] sp, ArrayList<POI> pois, ArrayList<Integer> itemsToBuy, Float travelWeight, Float costWeight){
        for (int current=0; current <currentToNextItemArray.length;current++)
        {
            for (int next=0;next<currentToNextItemArray[current].length;next++)
            {
                double dist = sp[current][next];
                for (int item=0;item<currentToNextItemArray[current][next].length;item++)
                {
                    if( !itemsToBuy.contains(item) || pois.get(next).items[item]== (double) 0) continue;
                    else{
                        double cost = pois.get(next).items[item];
                        currentToNextItemArray[current][next][item][0] = item;
                        currentToNextItemArray[current][next][item][1] = costWeight*cost + travelWeight*dist;
                    }
                }
            }
        }
    }

    public static BigDecimal [] MinCost(Route route, double [][][] startToPOIItemArray, Integer n){

        if(n == null){
            if(route.n.size() == 0){
                n = 0;
            }
            else{
                n = route.n.get(route.n.size()-1);
                if( n == startToPOIItemArray.length-1) return null;
            }
        }


        for(int i=0;i<route.n.size();i++){
            if(route.n.get(i) >= startToPOIItemArray.length){
                return null;
            }
        }

//        double [][] sortedItemsPerPOI = new double [startToPOIItemArray.length][1000];
        for(int i=0;i<startToPOIItemArray.length;i++){
            java.util.Arrays.sort(startToPOIItemArray[i], new java.util.Comparator<double[]>() {
                public int compare(double[] a, double[] b) {
                    return Double.compare(a[0], b[0]);
                }
            });

        }

        double [][] sortedPOIs = new double [startToPOIItemArray.length][3];

        for(int i=0;i<startToPOIItemArray.length;i++){

            sortedPOIs[i][0] = i;
            sortedPOIs[i][1] = startToPOIItemArray[i][0][0];
            sortedPOIs[i][2] = startToPOIItemArray[i][0][1];
        }
//        System.out.println(1);
        java.util.Arrays.sort(sortedPOIs, new java.util.Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                return Double.compare(a[2], b[2]);
            }
        });

        if(sortedPOIs[n][0] == Float.POSITIVE_INFINITY || sortedPOIs[n][1] == Float.POSITIVE_INFINITY || sortedPOIs[n][2] == Float.POSITIVE_INFINITY){ // in case you picked an item which is not being sold
            return null;
        }

        BigDecimal [] returnValue = new BigDecimal[3];

        returnValue[0] = new BigDecimal(""+sortedPOIs[n][0]);

        returnValue[1] = new BigDecimal(""+sortedPOIs[n][1]);

        returnValue[2] = new BigDecimal(""+sortedPOIs[n][2]);


        return returnValue;

    }

    public static BigDecimal [] getMinCost(Route route, double [][][] currentToNextItemArray){
        Integer n = route.n.get(route.n.size()-1);

        for(int i=0;i<route.n.size();i++){
            if(route.n.get(i) >= currentToNextItemArray.length){
                return null;
            }
        }

        double [][][] sortedItemsPerPOI = new double[currentToNextItemArray.length][currentToNextItemArray[0].length][];
//        ArrayList<ArrayList<ArrayList<ArrayList<Double>>>> sortedItemsPerPOI= new ArrayList<>();

        for (int i = 0; i < currentToNextItemArray.length; i++) {
            ArrayList<ArrayList<ArrayList<Double>>> temp = new ArrayList<ArrayList<ArrayList<Double>>>();
            for (int j = 0; j < currentToNextItemArray[i].length; j++){

                sortedItemsPerPOI[i][j] = Arrays.copyOf(currentToNextItemArray[i][j], currentToNextItemArray[i][j].length);
//                  ArrayList<ArrayList<Double>> temp2 = new ArrayList<ArrayList<Double>>();
//                  for(int k = 0; k < currentToNextItemArray[i][j].length; k++){
//                      ArrayList<Double> temp3 = new ArrayList<Double>();
//                      temp3.add(Double.valueOf(currentToNextItemArray[i][j][k]));
//                      temp2.add(temp3);
//                  }
//                temp.add(temp2);
            }
//            sortedItemsPerPOI.add(temp);


        }
//        ArrayList<Integer> temp = new ArrayList<Integer>();
//        for(int i=0;i<sortedItemsPerPOI.length;i++){
//            java.util.Arrays.sort(sortedItemsPerPOI[i], new java.util.Comparator<double[]>() {
//                public int compare(double[] a, double[] b) {
//                    return Double.compare(a[0], b[0]);
//                }
//            });
//
//        }




        double [][] sortedPOIs = new double [sortedItemsPerPOI.length][3];

        for(int i=0;i<sortedItemsPerPOI.length;i++){

            sortedItemsPerPOI = filterArray(sortedItemsPerPOI, route, i);
        }

        // sorting again because the sorted array got modified once.
        for(int i=0;i<sortedItemsPerPOI.length;i++){
            java.util.Arrays.sort(sortedItemsPerPOI[i], new java.util.Comparator<double[]>() {
                public int compare(double[] a, double[] b) {
                    return Double.compare(a[1], b[1]);
                }
            });

        }



        for(int i=0;i<sortedItemsPerPOI.length;i++){

            sortedPOIs[i][0] = i;

            sortedPOIs[i][1] = sortedItemsPerPOI[i][0][0];
            sortedPOIs[i][2] = sortedItemsPerPOI[i][0][1];
        }



//        System.out.println(1);
        java.util.Arrays.sort(sortedPOIs, new java.util.Comparator<double[]>() {
            public int compare(double[] a, double[] b) {
                return Double.compare(a[2], b[2]);
            }
        });

        if(sortedPOIs[n][0] == Float.POSITIVE_INFINITY || sortedPOIs[n][1] == Float.POSITIVE_INFINITY || sortedPOIs[n][2] == Float.POSITIVE_INFINITY){ // in case you picked an item which is not being sold
            return null;
        }

        BigDecimal [] returnValue = new BigDecimal[3];

        returnValue[0] = new BigDecimal(""+sortedPOIs[n][0]);

        returnValue[1] = new BigDecimal(""+sortedPOIs[n][1]);

        returnValue[2] = new BigDecimal(""+sortedPOIs[n][2]);


        return returnValue;




    }

    public static double [][][] filterArray( double [][][] sortedItemsPerPOI, Route route,Integer i){
        int j = 0;
        int length = sortedItemsPerPOI[i].length;
        while(j < length){

            if(sortedItemsPerPOI[i][j][1] == Float.POSITIVE_INFINITY){
                return sortedItemsPerPOI;
            }


            if( route.items.contains((int) sortedItemsPerPOI[i][j][0])){
                sortedItemsPerPOI[i][j][0] = Float.POSITIVE_INFINITY;
                sortedItemsPerPOI[i][j][1] = Float.POSITIVE_INFINITY;
                j++;
            }else{
                j++;
            }


        }
        return sortedItemsPerPOI;
    }

    public static boolean SatisfyList(Route route, ArrayList<Integer> itemsToBuy){
        for (int i=0;i<itemsToBuy.size();i++){
            if(!route.items.contains(itemsToBuy.get(i))){
                return false;
            }
        }
        return true;
    }

    public static Route nextItemRoute(Route route, double [][][][] currentToNextItemArray){

        Route newRoute = new Route(route);

        Integer currPoi = route.latestPoi;

        BigDecimal [] nextRouteInfo = new BigDecimal[3];
        nextRouteInfo = getMinCost(newRoute, currentToNextItemArray[currPoi].clone());

        if( nextRouteInfo == null || nextRouteInfo[0] == null || nextRouteInfo[1] == null || nextRouteInfo[2] == null){
            return null;
        }

        if( nextRouteInfo[0].floatValue() == Float.POSITIVE_INFINITY || nextRouteInfo[1].floatValue() == Float.POSITIVE_INFINITY || nextRouteInfo[2].floatValue() == Float.POSITIVE_INFINITY){
            return null;
        }

        newRoute.addPOI( nextRouteInfo[0].intValue());
        newRoute.addItem( nextRouteInfo[1].intValue());
        newRoute.addCost( (Float) nextRouteInfo[2].floatValue());

        return newRoute;
    }

    public static Route nextMinCostRoute(Route route, double [][][][] currentToNextItemArray){
        Route newRoute = new Route(route);
        newRoute.removeLatestPOI();
        newRoute.n.set(newRoute.n.size()-1, newRoute.n.get(newRoute.n.size()-1) + 1);
        Integer currPoi = route.latestPoi;

        BigDecimal [] nextRouteInfo = new BigDecimal[3];
        nextRouteInfo = getMinCost(newRoute, currentToNextItemArray[currPoi].clone());


        if( nextRouteInfo == null || nextRouteInfo[0] == null || nextRouteInfo[1] == null || nextRouteInfo[2] == null){
            return null;
        }

        if( nextRouteInfo[0].floatValue() == Float.POSITIVE_INFINITY || nextRouteInfo[1].floatValue() == Float.POSITIVE_INFINITY || nextRouteInfo[2].floatValue() == Float.POSITIVE_INFINITY){
            return null;
        }


        newRoute.addPOI( nextRouteInfo[0].intValue());
        newRoute.addItem( nextRouteInfo[1].intValue());
        newRoute.addCost( (Float) nextRouteInfo[2].floatValue());

        return newRoute;
    }

    public static Route findBestRoute(double [][][] startToPOIItemArray, double [][][][] currentToNextItemArray, double [] endToPOIs, ArrayList<Integer> itemsToBuy){

        Route potentialRoute = null;
        Comparator<Route> C = new Comp();
        PriorityQueue<Route> PQ = new PriorityQueue<Route>(10, C);
        Route startRoute = new Route();
        BigDecimal [] startRouteInfo = new BigDecimal[3];
        startRouteInfo = MinCost(startRoute, startToPOIItemArray.clone(), null); // creating incomplete route with lowest combined cost from starting loc

        startRoute.addPOI( startRouteInfo[0].intValue());
        startRoute.addItem( startRouteInfo[1].intValue());
        startRoute.addCost( (Float) startRouteInfo[2].floatValue());

        PQ.add(startRoute);

        int count = 0;
        while(!PQ.isEmpty()){
            Route r = new Route();
            r = PQ.poll();
            if(count == 72){
                int a = 0;
            }

            if(SatisfyList(r, itemsToBuy)){
                if(r.items.size() != itemsToBuy.size()){
                    System.out.println("Error");
                    for(int i = 0; i < r.items.size(); i++){
                        System.out.print("poi: "+ r.pois.get(i) + ", item: " + r.items.get(i) + ", weight: " + r.costs.get(i));
                    }
                    return null;
                }

                Integer latestPoi = r.latestPoi;
                BigDecimal finalWeight = new BigDecimal( "" + endToPOIs[latestPoi] ); //""+endToPOIs[latestPoi] because BigDecimal creates unnecessary decimal places otherwise
                r.addFinalCost(finalWeight.floatValue());
                potentialRoute = r;
                System.out.println("First Alg Count: " + count);
                return potentialRoute;
            }
            Route nextR = nextItemRoute(r, currentToNextItemArray);

            Route nextMinCostR = new Route(r);

            if(r.pois.size() == 1){
                nextMinCostR.removeLatestPOI();
                Integer n = r.n.get(0) + 1;

                if(n >= currentToNextItemArray.length){ // incase we explore every single possible route from start to a POI
                    if(nextR != null){
                        PQ.add(nextR);
                        count++;
                    }
                    continue;
                }
                BigDecimal [] startRouteNewInfo = new BigDecimal[3];
                startRouteNewInfo = MinCost(nextMinCostR, startToPOIItemArray.clone(), n);

                if( startRouteNewInfo == null || startRouteNewInfo[0] == null || startRouteNewInfo[1] == null || startRouteNewInfo[2] == null){
                    if(nextR != null){
                        PQ.add(nextR);
                        count++;
                    }
                    continue;
                }

                if( startRouteNewInfo[0].floatValue() == Float.POSITIVE_INFINITY || startRouteNewInfo[1].floatValue() == Float.POSITIVE_INFINITY || startRouteNewInfo[2].floatValue() == Float.POSITIVE_INFINITY){
                    if(nextR != null){
                        PQ.add(nextR);
                        count++;
                    }
                    continue;
                }

                nextMinCostR.addPOI( startRouteNewInfo[0].intValue());
                nextMinCostR.n.set(0, n);
                nextMinCostR.addItem( startRouteNewInfo[1].intValue());
                nextMinCostR.addCost( (Float) startRouteNewInfo[2].floatValue());

            } else{
                nextMinCostR = nextMinCostRoute(r, currentToNextItemArray.clone());
            }

            if(nextR != null){
                PQ.add(nextR);
            }
            if(nextMinCostR != null){
                PQ.add(nextMinCostR);
            }
            count++;

        }


        System.out.println("First Alg Count: " + count);
        return potentialRoute;

    }

    public static Route findBestRouteContinued(double [][][] startToPOIItemArray, double [][][][] currentToNextItemArray, double [] endToPOIs, ArrayList<Integer> itemsToBuy){

        Route potentialRoute = null;
        Comparator<Route> C = new Comp();
        PriorityQueue<Route> PQ = new PriorityQueue<Route>(10, C);
        Route startRoute = new Route();
        BigDecimal [] startRouteInfo = new BigDecimal[3];
        startRouteInfo = MinCost(startRoute, startToPOIItemArray.clone(), null);

        startRoute.addPOI( startRouteInfo[0].intValue());
        startRoute.addItem( startRouteInfo[1].intValue());
        startRoute.addCost( (Float) startRouteInfo[2].floatValue());

        PQ.add(startRoute);

        Float weightUpperBound = Float.POSITIVE_INFINITY;
        int count = 0;
        int pruned = 0;
        while(!PQ.isEmpty()){
            Route r = new Route();
            r = PQ.poll();

            if(r.totalCost > weightUpperBound){
                pruned++;
                continue;
            }

            if(SatisfyList(r, itemsToBuy)){
                if(r.items.size() != itemsToBuy.size()){
                    System.out.println("Error");
                    for(int i = 0; i < r.items.size(); i++){
                        System.out.print("poi: "+ r.pois.get(i) + ", item: " + r.items.get(i) + ", weight: " + r.costs.get(i));
                    }
                    return null;
                }

                Integer latestPoi = r.latestPoi;
                BigDecimal finalWeight = new BigDecimal( "" + endToPOIs[latestPoi] ); //""+endToPOIs[latestPoi] because BigDecimal creates unnecessary decimal places otherwise
                r.addFinalCost(finalWeight.floatValue());


                if(potentialRoute == null){
                    potentialRoute = r;
                    weightUpperBound = r.totalCost;
                }else{
                    if(r.totalCost < potentialRoute.totalCost){
                        potentialRoute = r;
                        weightUpperBound = r.totalCost;
                    }
                }
                continue;
//                return potentialRoute;
            }
            Route nextR = nextItemRoute(r, currentToNextItemArray);

            Route nextMinCostR = new Route(r);

            if(r.pois.size() == 1){
                nextMinCostR.removeLatestPOI();
                Integer n = r.n.get(0) + 1;

                if(n >= currentToNextItemArray.length){ // incase we explore every single possible route from start to a POI
                    if(nextR != null){
                        PQ.add(nextR);
                        count++;
                    }
                    continue;
                }

                BigDecimal [] startRouteNewInfo = new BigDecimal[3];
                startRouteNewInfo = MinCost(nextMinCostR, startToPOIItemArray.clone(), n);

                if( startRouteNewInfo == null || startRouteNewInfo[0] == null || startRouteNewInfo[1] == null || startRouteNewInfo[2] == null){
                    if(nextR != null){
                        PQ.add(nextR);
                        count++;
                    }
                    continue;
                }

                if( startRouteNewInfo[0].floatValue() == Float.POSITIVE_INFINITY || startRouteNewInfo[1].floatValue() == Float.POSITIVE_INFINITY || startRouteNewInfo[2].floatValue() == Float.POSITIVE_INFINITY){
                    if(nextR != null){
                        PQ.add(nextR);
                        count++;
                    }
                    continue;
                }

                nextMinCostR.addPOI( startRouteNewInfo[0].intValue());
                nextMinCostR.n.set(0, n);
                nextMinCostR.addItem( startRouteNewInfo[1].intValue());
                nextMinCostR.addCost( (Float) startRouteNewInfo[2].floatValue());

            } else{
                nextMinCostR = nextMinCostRoute(r, currentToNextItemArray.clone());
            }

            if(nextR != null && nextR.totalCost < weightUpperBound){
                PQ.add(nextR);
            }
            if(nextMinCostR != null && nextMinCostR.totalCost < weightUpperBound){
                PQ.add(nextMinCostR);
            }
            count++;

        }


        System.out.println("Second Alg Count: " + count);
        System.out.println("Second Alg Pruned: " + pruned);

        return potentialRoute;

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

    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
	// write your code here

        Float travelWeight = 0.0f; // the weight that the user gives to travel
        Float costWeight = 1 - travelWeight; // the weight that the user gives to cost

        ArrayList<Integer> itemsToBuy = new ArrayList<Integer>();

        for (int i = 0;i<5;i++){
            itemsToBuy.add(i); // creating a shopping list
        }


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

        double [][][] startToPOIItemArray = new double [pois.size()][items.size()][2]; // start to [poi][item][0] is the itemID, start to [poi][item][1] is the combined weight of the cost of item and travel distance

        // set startToPOIItemArray elements to infinity. Can't use null because of multiple type cast conflict. Infinity works when we want to sort by lowest cost
        for (int i=0;i<pois.size();i++){
            for (int j=0;j<items.size();j++){
                startToPOIItemArray[i][j][0] = Double.POSITIVE_INFINITY;
                startToPOIItemArray[i][j][1] = Double.POSITIVE_INFINITY;
            }
        }

        createStart3DArray(startToPOIItemArray, startToPOIs, pois, itemsToBuy, travelWeight, costWeight); // create 3D array containing the combined weight of cost and travel dist from start loc to all pois and items


        double [][][][] currentToNextItemArray = new double [pois.size()][pois.size()][items.size()][2]; // current poi to [next poi][item][0] is the itemID, current poi to [next poi][item][1] is the combined weight of the cost of item and travel distance

        // set currentToNextItemArray elements to infinity. Can't use null because of multiple type cast conflict. Infinity works when we want to sort by lowest cost
        for (int i=0;i<pois.size();i++){
            for (int j=0;j<pois.size();j++){
                for (int k=0;k<items.size();k++){
                    currentToNextItemArray[i][j][k][0] = Double.POSITIVE_INFINITY;
                    currentToNextItemArray[i][j][k][1] = Double.POSITIVE_INFINITY;
                }
            }

        }

        create3DArray(currentToNextItemArray, sp, pois, itemsToBuy, travelWeight, costWeight); // techincally creating a 4D array from one POi to another POi to buy item X. [current][next][item] has 2 values: [0] is the itemID, [1] is the combined weight of the cost of item and travel distance


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

        long preComputationEndTime = System.currentTimeMillis();
        double preComputationDuration = (double) (preComputationEndTime - preComputationStartTime)/ (double) 1000;

        System.out.println("Precomputation time: " + preComputationDuration);

        double firstALgAvgDuration = 0;
        double firstAlgAvgTotalCost = 0;

        double secondALgAvgDuration = 0;
        double secondAlgAvgTotalCost = 0;

        int n = 5; // number of times to run the algorithm
        for (int i=0;i<n;i++) {
            long firstALgStartTime = System.currentTimeMillis();

            Route firstRoute;

            firstRoute = findBestRoute(startToPOIItemArray.clone(), currentToNextItemArray.clone(), endToPOIs.clone(), itemsToBuy); // find the first best route from start to end
            System.out.println("Route: " + firstRoute.toString());


            long firstAlgEndTime = System.currentTimeMillis();
            System.out.println("First algorithm time: " + (double) (firstAlgEndTime - firstALgStartTime) / (double) 1000);

            firstALgAvgDuration += (double) (firstAlgEndTime - firstALgStartTime) / (double) 1000;
            firstAlgAvgTotalCost += firstRoute.totalCost;

//            long secondALgStartTime = System.currentTimeMillis();
//
//            Route secondRoute;
//
//            secondRoute = findBestRouteContinued(startToPOIItemArray.clone(), currentToNextItemArray.clone(), endToPOIs.clone(), itemsToBuy);
//            System.out.println("Route: " + secondRoute.toString());
//
//
//            long secondAlgEndTime = System.currentTimeMillis();
//            System.out.println("First algorithm time: " + (double) (secondAlgEndTime - secondALgStartTime) / (double) 1000);
//            System.out.println("");
//            System.out.println("");
//
//            secondALgAvgDuration += (double) (secondAlgEndTime - secondALgStartTime) / (double) 1000;
//            secondAlgAvgTotalCost += secondRoute.totalCost;

        }

        System.out.println("Average first algorithm time: " + (double) firstALgAvgDuration/ (double) n);
//        System.out.println("Average second algorithm time: " + (double) secondALgAvgDuration/ (double) n);

        System.out.println("Average first algorithm total cost: " + (double) firstAlgAvgTotalCost/ (double) n);
//        System.out.println("Average second algorithm total cost: " + (double) secondAlgAvgTotalCost/ (double) n);

    }
}
