import processing.core.PImage;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Graph {
    //region data
    //enum to specify type of graph, used for file naming
    enum GraphTypes{
        segment,
        shortest,
        widest,
        kruskal,
        roza
    }
    final GraphTypes type;
    //nodes contained in the graph
    ArrayList<Node> nodes;
    //stores the edges, edges are made in pairs (i.e. edge from a to b and edge from b to a), hashset to prevent duplicate edges between nodes
    ArrayList<HashSet<Edge>> edges;
    //number of edges (not pairs), debug field
    int edgeCount = 0;
    //mid-point of the image
    Point midPoint;
    //endregion
    //region constructors
    //constructor to create spanning trees
    public Graph(ArrayList<Node> nodes, Point midPoint, GraphTypes type) {
        this.nodes = nodes;
        this.type = type;
        this.midPoint = midPoint;
        initializeEdges(nodes.size());
    }
    //constructor for the color segmentation
    public Graph(LinkedList<HashSet<Point>> points, HashMap<Point, EdgeData> neighbors, PImage originalImage) {
        type = GraphTypes.segment;
        midPoint = new Point(originalImage.width/2,originalImage.height/2);
        createNodes(points, originalImage);
        long start = System.currentTimeMillis();
        initializeEdges(nodes.size());
        for (Map.Entry<Point, EdgeData> map: neighbors.entrySet()){
            addEdge(map.getKey(), map.getValue());
        }
        System.out.println("Graph: " + edgeCount + " Edges Created: time elapsed: " + (System.currentTimeMillis() - start));
    }
    //constructor for average segmentation
    public Graph(LinkedList<HashSet<Point>> points, HashSet<Point> neighbors, PImage originalImage) {
        type = GraphTypes.segment;
        midPoint = new Point(originalImage.width/2,originalImage.height/2);
        createNodes(points, originalImage);
        long start = System.currentTimeMillis();
        initializeEdges(nodes.size());
        for (Point p: neighbors){
            addEdge(p.x,p.y);
        }
        System.out.println("Graph: " + edgeCount + " Edges Created: time elapsed: " + (System.currentTimeMillis() - start));
    }
    //function to initialize the nodes in the graph based on the points in the segmentation
    private void createNodes(LinkedList<HashSet<Point>> points, PImage originalImage) {
        long start = System.currentTimeMillis();
        nodes = new ArrayList<>(points.size());
        AtomicInteger currentId = new AtomicInteger();
        points.forEach(n -> nodes.add(new Node(n, originalImage, currentId.getAndIncrement())));
        System.out.println("Graph: " + nodes.size() + " Nodes Created: time elapsed: " + (System.currentTimeMillis() - start));
    }
    //initializes the structure for edges
    private void initializeEdges(int size) {
        edges = Stream.generate(HashSet<Edge>::new).limit(size).collect(Collectors.toCollection(ArrayList::new));
    }
    //node id in partial graph doesn't line up to index, use partial id instead to reference nodes in the partial graph, possible better implementation to get around this
    public static Graph partialGraph(HashSet<Node> nodes, Graph basis, Node added) {
        ArrayList<Node> node = new ArrayList<>(nodes);
        node.add(added);
        Graph g = new Graph(node, new Point(-1,-1), GraphTypes.segment);
        for (int i = 0; i < g.nodes.size(); i++) {
            Node n = g.nodes.get(i);
            n.partialId = i;
            for (Edge e : basis.edges.get(n.nodeId)) {
                if (nodes.contains(e.nodeA) && nodes.contains(e.nodeB)) {
                    g.addEdgePartial(e);
                }
            }
        }
        return g;
    }
    //endregion
    //region addEdge
    private void addEdge(Point p, EdgeData data) {
        data.calcWeight();
        edgeCount++;
        edges.get(p.x).add(new Edge(nodes.get(p.x), nodes.get(p.y), data));
        edges.get(p.y).add(new Edge(nodes.get(p.y), nodes.get(p.x), data));
    }
    public void addEdge(int x, int y, EdgeData data) {
        edgeCount++;
        data.calcWeight();
        edges.get(x).add(new Edge(nodes.get(x), nodes.get(y), data));
        edges.get(y).add(new Edge(nodes.get(y), nodes.get(x), data));
    }
    public void addEdge(int x, int y) {
        edgeCount++;
        EdgeData e = new EdgeData(Color.rozaColorDiff(nodes.get(x).avgColor, nodes.get(y).avgColor));
        edges.get(x).add(new Edge(nodes.get(x), nodes.get(y), e));
        edges.get(y).add(new Edge(nodes.get(y), nodes.get(x), e));
    }
    public void addEdge(Node n1, Node n2) {
        addEdge(n1.nodeId, n2.nodeId);
    }
    public void addEdge(Edge e) {
        edgeCount++;
        edges.get(e.nodeA.nodeId).add(e);
        Edge e2 = new Edge(e.nodeB, e.nodeA, e.data);
        edges.get(e.nodeB.nodeId).add(e2);
    }
    void addEdgePartial(Edge e){
        edgeCount++;
        edges.get(e.nodeA.partialId).add(e);
        edges.get(e.nodeB.partialId).add(e);
    }
    //endregion
    //returns the region that contains the center of the image
    public Node findMidReg() {
        for (Node n : nodes) {
            if (n.points.contains(midPoint)) {
                return n;
            }
        }
        return null;
    }
    public int size() {
        return nodes.size();
    }
    //returns a list of all neighbors of a node
    public LinkedList<Node> getNeighbors(Node n) {
        return edges.get(n.nodeId).stream().map(Edge::getNodeB).collect(Collectors.toCollection(LinkedList::new));
    }
    //returns a formatted string of all of a nodes neighbors
    public String getNeighborsText(Node node) {
        StringBuilder builder = new StringBuilder();
        for (Node n : getNeighbors(node)) {
            builder.append(n);
            builder.append(", ");
        }
        return builder.toString();
    }
    //prints all the edges in the graph
    public void print() {
        for (int i = 0; i < edges.size(); i++) {
            System.out.print("Node " + i + " neighbors : ");
            for (Edge e : edges.get(i)) {
                System.out.print(e);
            }
            System.out.println();
        }
    }
    //returns and array containing all the edge weights (includes pairs)
    public double[] getEdgeWeights() {
        return edges.stream().flatMap(Collection::stream).mapToDouble(Edge::getWeight).toArray();
    }
    //no copy for partial graphs
    public Graph copy(){
        ArrayList<Node> newNodes = this.nodes.stream().map(Node::copy).collect(Collectors.toCollection(ArrayList::new));
        Graph newGraph = new Graph(newNodes, this.midPoint, this.type);
        newGraph.transferEdges(this.edges);
        return newGraph;
    }

    private void transferEdges(ArrayList<HashSet<Edge>> edges) {
        edges.stream().flatMap(Collection::stream).forEach(this::copyEdge);
    }
    void copyEdge(Edge e){
        edgeCount++;
        EdgeData data = e.data.copy();
        addEdge(e.nodeA.nodeId, e.nodeB.nodeId, data);
    }
    private void transferEdgesPartial(List<Edge> edges){
        edges.forEach(this::addEdgePartial);
    }
}
