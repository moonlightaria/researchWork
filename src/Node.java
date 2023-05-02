import processing.core.PImage;

import java.util.HashSet;
import java.util.Set;

public class Node{
    //set of pixels in the segment
    HashSet<Point> points;
    //list of all pixels on the edge scaled to their position on the interface, used to show outline of segments
    HashSet<Point> scaledVertices;
    //sum of all edge weights for edges on this node, used for the priority queue when creating spanning trees
    int sumWeight;
    //id correlating to its position in the nodes and edges fields in the graph
    int nodeId;
    //equivalent to node id for a partial graph
    int partialId;
    //average color of the segment in the original image
    Color avgColor;
    //center most point of the segment, used for connection view along with segment selection
    Point centroid;

    //fields used to store recoloring information
    Color finalColor = null;
    boolean dontColor = false;
    Set<Color> blockColor = new HashSet<>();
    Set<Node> shareColor = new HashSet<>();
    Set<Node> dontShareColor = new HashSet<>();
    Set<Node> biasTowardsSegment = new HashSet<>();
    Set<Node> biasFromSegment = new HashSet<>();
    Set<Color> biasTowardsColor = new HashSet<>();
    Set<Color> biasFromColor = new HashSet<>();
    Set<Node> furthestColorFromSegment = new HashSet<>();

    //fields for storing information during recoloring
    Color color;
    boolean colored;
    Node parent;

    public Node(HashSet<Point> Points, PImage image, int id){
        points = new HashSet<>(Points);
        nodeId = id;
        sumWeight =0;
        int r=0,g=0,b=0;
        Color c;
        for (Point p : points) {
            c = new Color(image.get(p.x,p.y));
            r += c.r;
            g += c.g;
            b += c.b;
        }
        int size = points.size();
        avgColor = new Color (r/size, g/size,b/size, 255);
    }
    public Node() {}
    //finds all the points on the edge of the segment
    private HashSet<Point> outline() {
        HashSet<Point> outline = new HashSet<>();
        for (Point p : points) {
            for (Point check: Algorithms.surroundingPoints(p)){
                if (!points.contains(check)){
                    outline.add(check);
                }
            }
        }
        return outline;
    }
    //scales all the vertices based on the interface scale factor
    public void scaledOutline(double scaleFactor) {
        scaledVertices = outline();
        for (Point p: scaledVertices){
            p.x = (int) Math.floor(p.x*scaleFactor);
            p.y = (int) Math.floor(p.y*scaleFactor);
        }
    }

    @Override
    public String toString(){
        return "Node: " + nodeId;
    }
    //does not copy the recoloring information
    Node copy(){
        Node newNode = new Node();
        newNode.nodeId = this.nodeId;
        newNode.sumWeight = this.sumWeight;
        newNode.centroid = this.centroid.copy();
        newNode.avgColor = this.avgColor.copy();
        newNode.points = new HashSet<>(this.points.size());
        for (Point p: this.points) {
            newNode.points.add(p.copy());
        }
        newNode.scaledVertices = new HashSet<>();
        for (Point p: this.scaledVertices) {
            newNode.scaledVertices.add(p.copy());
        }
        return newNode;
    }
}
