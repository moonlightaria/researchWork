public class EdgeData {
    //sum of all color distances with adjacent pixels
    double sumWeight;
    //weight of the edge
    double weight;
    //number of connections to other pixels
    int connections;
    public EdgeData(double weight) {
        this.weight = weight;
    }
    public EdgeData(int Connections, double sumWeight) {
        this.sumWeight = sumWeight;
        connections = Connections;
    }
    public EdgeData() {
        connections =0;
        sumWeight =0;
    }
    //calculates the weight of the node using the color segmentation method if the data for it exists
    public void calcWeight(){
        if (sumWeight > 0) {
            weight = sumWeight / connections;
        }
    }
    //return a new edge data with the info of the arguments merged (prob could be optimized)
    public static EdgeData merge(EdgeData data, EdgeData v) {
        return new EdgeData(data.connections + v.connections, data.sumWeight + v.sumWeight);
    }
    EdgeData copy(){
        EdgeData copy = new EdgeData();
        copy.connections = this.connections;
        copy.weight = this.weight;
        copy.sumWeight = this.sumWeight;
        return copy;
    }
}
