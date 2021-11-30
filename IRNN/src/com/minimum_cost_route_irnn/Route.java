package com.minimum_cost_route_irnn;

import java.util.ArrayList;


public class Route {
    ArrayList<Integer> pois = new ArrayList<Integer>();
    ArrayList<Integer> items = new ArrayList<Integer>();
    ArrayList<Float> costs = new ArrayList<Float>();
    Integer latestPoi = null;
    Integer latestItem = null;
    Float latestCost =  null;

    Float totalCost = (float) 0.0;
    ArrayList<Integer> n = new ArrayList<Integer>();

    public Route(Route r) {
        for (int i=0;i<r.pois.size();i++)
        {
            this.pois.add(r.pois.get(i));
        }

        for (int i=0;i<r.items.size();i++)
        {
            this.items.add(r.items.get(i));
        }

        for (int i=0;i<r.costs.size();i++)
        {
            this.costs.add(r.costs.get(i));
        }

        for(int i=0;i<r.n.size();i++){
            this.n.add(r.n.get(i));
        }

        this.latestPoi = r.latestPoi;
        this.latestItem = r.latestItem;
        this.latestCost = r.latestCost;
        this.totalCost = r.totalCost;

    }
    public Route(){}

    public String toString() {
        return "poi: "+ this.pois + ", items: " + this.items + ", costs: " + this.costs + ", totalCost: " + this.totalCost + ", n: " + this.n;
    }

    public void addPOI(Integer poi){
        this.pois.add(poi);
        this.latestPoi = poi;
        this.n.add(0);
        return;
    }

    public void addItem(Integer item){
        this.items.add(item);
        this.latestItem = item;
        return;
    }

    public void addCost(Float cost){
        this.costs.add(cost);
        this.incrementTotalCost(cost);
        return;
    }

    public void incrementTotalCost(Float cost){
        this.totalCost += cost;
        this.latestCost = cost;
        return;
    }

    public void addFinalCost(Float cost){
        this.totalCost += cost;
        return;
    }

    public void removeLatestPOI(){
        this.pois.remove(this.pois.size()-1);
        this.items.remove(this.items.size()-1);
        Float cost = this.costs.remove(this.costs.size()-1);
        this.totalCost -= cost;
        this.n.remove(this.n.size()-1);

        if(this.pois.size() == 0){
            this.latestPoi = null;
            this.latestItem = null;
            this.latestCost = null;
            return;
        }

        this.latestPoi = this.pois.get(this.pois.size()-1);
        this.latestItem = this.items.get(this.items.size()-1);
        this.latestCost = this.costs.get(this.costs.size()-1);

        return;
    }

    public void updateN(Integer index, Integer value){
        if( (index - this.pois.size() - 1) >= 2){
            System.out.println("These is a gap of more than 2: " + (index - this.pois.size()));
        }
        this.n.set(index, value);
        return;
    }

}
