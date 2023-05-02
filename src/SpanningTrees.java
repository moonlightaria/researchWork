import java.util.*;
//class that stores the various methods of creating spanning trees
public class SpanningTrees {

    //https://www.geeksforgeeks.org/widest-path-problem-practical-application-of-dijkstras-algorithm/
    public static Graph widestPath(Graph graph, Node node){
        PriorityQueue<Node> queue = new PriorityQueue<>(graph.size(), (n1, n2) -> n2.sumWeight - n1.sumWeight);
        double[] widest = new double[graph.size()];
        Edge[] result = new Edge[graph.size()];

        Arrays.fill(widest, Integer.MIN_VALUE);
        queue.add(node);
        widest[node.nodeId] = Integer.MAX_VALUE;

        while (!queue.isEmpty()){
            Node n = queue.poll();
            maxDistance(n, queue, graph, widest, result);
        }
        return createGraph(Graph.GraphTypes.widest, graph, result);
    }
    private static void maxDistance(Node n, PriorityQueue<Node> queue, Graph graph, double[] widest, Edge[] result) {
        double newWidest;
        for (Edge edge : graph.edges.get(n.nodeId)) {
            newWidest = Math.max(widest[edge.nodeB.nodeId], Math.min(widest[edge.nodeA.nodeId], edge.getWeight()));
            if (newWidest > widest[edge.nodeB.nodeId]) {
                widest[edge.nodeB.nodeId] = newWidest;
                result[edge.nodeB.nodeId] = edge; //new Edge(edge.nodeA, edge.nodeB, new EdgeData(newWidest));
                queue.add(edge.nodeB);
            }
        }
    }
    /*https://www.geeksforgeeks.org/java-program-for-dijkstras-shortest-path-algorithm-greedy-algo-7/?ref=lbp
    https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-in-java-using-priorityqueue/
    https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-using-priority_queue-stl/*/
    public static Graph shortestPath(Graph graph, Node node){
        PriorityQueue<Node> queue = new PriorityQueue<>(graph.size(), (n1, n2) -> n2.sumWeight - n1.sumWeight);
        double[] distance = new double[graph.size()];
        boolean[] finished = new boolean[graph.size()];
        Edge[] result = new Edge[graph.size()];

        Arrays.fill(distance, Integer.MAX_VALUE);

        queue.add(node);
        distance[node.nodeId] = 0;
        while (!queue.isEmpty() && containsFalse(finished)){
            Node n = queue.poll();
            minDistance(n, queue, graph, finished, distance, result);
        }
        return createGraph(Graph.GraphTypes.shortest, graph, result);
    }
    private static boolean containsFalse(boolean[] finished) {
        for (boolean b: finished) {
            if (!b){
                return true;
            }
        }
        return false;
    }
    private static void minDistance(Node n, PriorityQueue<Node> queue, Graph graph, boolean[] finished, double[] distance, Edge[] result) {
        double newDistance;
        for (Edge edge : graph.edges.get(n.nodeId)) {
            newDistance = distance[n.nodeId] + edge.getWeight();
            if (newDistance < distance[edge.nodeB.nodeId]) {
                distance[edge.nodeB.nodeId] = newDistance;
                result[edge.nodeB.nodeId] = edge;//new Edge(edge.nodeA, edge.nodeB, new EdgeData(newDistance));
            }
            if (!finished[edge.nodeB.nodeId]) {
                finished[edge.nodeB.nodeId] = true;
                queue.add(edge.nodeB);
            }
        }
    }
    //https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-algorithm-greedy-algo-2/
    public static Graph kruskal(Graph graph){
        Edge[] result = new Edge[graph.size()];
        int e =0;
        ArrayList<Edge> edges = new ArrayList<>();
        for (HashSet<Edge> hashSet: graph.edges) {
            edges.addAll(hashSet);
        }
        edges.sort(Comparator.comparingDouble(Edge::getWeight));
        //x=parent, y=rank
        Point[] subsets = new Point[graph.size()];
        for (int i = 0; i < subsets.length; i++) {
            subsets[i] = new Point(i,0);
        }
        int currentEdge = 0;
        while (e < graph.size() - 1) {
            Edge next_edge = edges.get(currentEdge++);
            int x = find(subsets, next_edge.nodeA.nodeId);
            int y = find(subsets, next_edge.nodeB.nodeId);
            if (x != y) {
                result[e++] = next_edge;
                Union(subsets, x, y);
            }
        }
        return createGraph(Graph.GraphTypes.kruskal, graph, result);
    }
    static int find(Point[] subsets, int i) {
        // find root and make root as parent of i
        // (path compression)
        if (subsets[i].x != i) {
            subsets[i].x= find(subsets, subsets[i].x);
        }
        return subsets[i].x;
    }
    static void Union(Point[] subsets, int x, int y) {
        int xroot = find(subsets, x);
        int yroot = find(subsets, y);

        if (subsets[xroot].y < subsets[yroot].y) {
            subsets[xroot].x= yroot;
        }else if (subsets[xroot].y > subsets[yroot].y) {
            subsets[yroot].x = xroot;
        }else {
            subsets[yroot].x = xroot;
            subsets[xroot].y++;
        }

    }
    private static Graph createGraph(Graph.GraphTypes type, Graph graph, Edge[] result){
        Graph returnGraph = new Graph(graph.nodes, graph.midPoint, type);
        for (Edge edge : result) {
            if (edge != null) {
                returnGraph.addEdge(edge);
            }
        }
        return returnGraph;
    }
    //spanning tree algorithm used by roza in her code
    public static Graph rozaGraph(Graph graph, Node node) {
        PriorityQueue<Node> queue = new PriorityQueue<>(graph.size(), (n1, n2) -> n2.sumWeight - n1.sumWeight);
        double[] widest = new double[graph.size()];
        Edge[] result = new Edge[graph.size()];
        boolean[] finished = new boolean[graph.size()];
        Arrays.fill(widest, Integer.MIN_VALUE);
        queue.add(node);
        widest[node.nodeId] = Integer.MAX_VALUE;

        double avesize = CalcAverage(graph.nodes);
        double midWeight = getMidWeight(graph.getEdgeWeights());
        while (!queue.isEmpty()){
            Node n = queue.poll();
            finished[n.nodeId] = true;
            rozaGraphEdgeCheck(n, queue, graph, widest, finished, result, midWeight, avesize);
        }
        return createGraph(Graph.GraphTypes.roza, graph, result);
    }
    //gets the median edge weight
    static double getMidWeight(double[] weights){
        if (weights.length <=1){
            return 0;
        }else{
            double maxw = Double.MIN_VALUE, minw = Double.MAX_VALUE;
            for (double d: weights){
                if (d >maxw){
                    maxw = d;
                }else if (d<minw){
                    minw = d;
                }
            }
            return (maxw + minw) / 2.0;
        }
    }
    //gets the average segment size
    static double CalcAverage(ArrayList<Node> nodes) {
        double avesize=0.0;
        for (Node n : nodes) {
            avesize += n.points.size();
        }
        return (avesize / nodes.size());
    }
    private static void rozaGraphEdgeCheck(Node n, PriorityQueue<Node> queue, Graph graph, double[] widest, boolean[] finished, Edge[] result, double midWeight, double avesize) {
        double newWidest;
        for (Edge edge : graph.edges.get(n.nodeId)) {
            newWidest = Math.abs(edge.getWeight() - midWeight)  + (edge.nodeB.points.size() / avesize) * Math.abs(edge.getWeight() - midWeight);
            newWidest = Math.min(widest[edge.nodeA.nodeId], newWidest);
            if (newWidest > widest[edge.nodeB.nodeId]) {
                if (!finished[edge.nodeB.nodeId]) {
                    widest[edge.nodeB.nodeId] = newWidest;
                    result[edge.nodeB.nodeId] = edge;//new Edge(edge.nodeA, edge.nodeB, new EdgeData(newWidest));
                    queue.add(edge.nodeB);
                }
            }
        }
    }
    //https://cs.stackexchange.com/questions/106617/multiple-source-shortest-paths-in-a-weighted-graph
    //variation of shortest path to work for multiple starting nodes by adding a temporary node to start from, doesn't work as intended
    //ends up with many separated spanning trees and what appears to be nodes without edges (not properly tested)
    public static Graph partialShortestPath(Graph graph, HashSet<Node> starting, Node added) {
        PriorityQueue<Node> queue = new PriorityQueue<>(graph.size(), (n1, n2) -> n2.sumWeight - n1.sumWeight);
        double[] distance = new double[graph.size()];
        boolean[] finished = new boolean[graph.size()];
        Edge[] result = new Edge[graph.size()];

        Arrays.fill(distance, Integer.MAX_VALUE);

        queue.add(added);
        for (Node n: starting){
            distance[n.partialId] = 0;
            graph.addEdge(added.partialId, n.partialId, new EdgeData(0));
        }
        while (containsFalse(finished) && !queue.isEmpty()){
            Node n = queue.poll();
            finished[n.partialId] = true;
            partialMinDistance(n, queue, graph, finished, distance, result);
        }
        return createPartialGraph(Graph.GraphTypes.shortest, graph, result, added);
    }
    private static void partialMinDistance(Node n, PriorityQueue<Node> queue, Graph graph, boolean[] finished, double[] distance, Edge[] result) {
        double newDistance;
        for (Edge edge : graph.edges.get(n.partialId)) {
            newDistance = distance[n.partialId] + edge.getWeight();
            if (newDistance < distance[edge.nodeB.partialId]) {
                distance[edge.nodeB.partialId] = newDistance;
                result[edge.nodeB.partialId] = edge;//new Edge(edge.nodeA, edge.nodeB, new EdgeData(newDistance));
            }
            if (!finished[edge.nodeB.partialId]) {
                queue.add(edge.nodeB);
            }
        }
    }
    private static Graph createPartialGraph(Graph.GraphTypes type, Graph graph, Edge[] result, Node added){
        graph.nodes.remove(added);
        Graph returnGraph = new Graph(graph.nodes, graph.midPoint, type);
        for (Edge edge : result) {
            if (edge != null) {
                if (edge.nodeA != added && edge.nodeB != added) {
                    returnGraph.addEdgePartial(edge);
                }
            }
        }
        return returnGraph;
    }
}