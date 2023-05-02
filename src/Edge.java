import java.util.Comparator;

public class Edge implements Comparator<Edge> {
    Node nodeA;
    Node nodeB;

    EdgeData data;

    public Edge(Node a, Node b, EdgeData edgeData){
        nodeA = a;
        nodeB = b;
        data = edgeData;
        nodeA.sumWeight += data.weight;
        nodeB.sumWeight += data.weight;
    }
    public Node getNodeA() {
        return nodeA;
    }
    public void setNodeA(Node nodeA) {
        this.nodeA = nodeA;
    }
    public Node getNodeB() {
        return nodeB;
    }
    public void setNodeB(Node nodeB) {
        this.nodeB = nodeB;
    }
    public double getWeight() {
        return data.weight;
    }
    public void setWeight(int weight) {
        data.weight = weight;
    }
    @Override
    public String toString(){
        return nodeA.nodeId + ":" + nodeB.nodeId + "(" + data.sumWeight + ") ";
    }
    @Override
    public int hashCode(){
        return nodeA.nodeId + nodeB.nodeId;
    }
    @Override
    public int compare(Edge e1, Edge e2) {
        if (e1.getWeight() < e2.getWeight())
            return -1;
        if (e1.getWeight() > e2.getWeight())
            return 1;
        return 0;
    }
    @Override
    public boolean equals(Object o){
        if (this == o) return true;
        if (o == null) return false;
        if (getClass() != o.getClass()) return false;
        Edge e = (Edge) o;
        if (this.nodeA == e.nodeA){
            return this.nodeB == e.nodeB;
        }else if(this.nodeA == e.nodeB){
            return this.nodeB == e.nodeA;
        }
        return false;
    }
}
