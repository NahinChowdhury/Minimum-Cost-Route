package com.minimum_cost_route_irnn;

import java.util.Comparator;

class Comp implements Comparator<Route> {
    public int compare(Route r1, Route r2)
    {
        if (r1.totalCost>r2.totalCost)return -1;
        else if (r1.totalCost<r2.totalCost) return 1;
        return 0;

    }

}


//class Comp_Dist implements Comparator<Route> {
//    public int compare(Route r1, Route r2)
//    {
//        if (r1.distance<r2.distance)return -1;
//        else if (r1.distance>r2.distance) return 1;
//        else return 0;
//
//    }
//
//}


//class Nei_Comp implements Comparator<Neighbor> {
//    public int compare(Neighbor n1, Neighbor n2)
//    {
//        if (n1.distance<n2.distance)return -1;
//        else if (n1.distance>n2.distance) return 1;
//        return 0;
//
//    }
//}


//class Pair_Comp implements Comparator<IS_Pair> {
//    public int compare(IS_Pair r1, IS_Pair r2)
//    {
//        if (r1.cost<r2.cost)return -1;
//        else if (r1.cost>r2.cost) return 1;
//        return 0;
//
//    }
//
//}
