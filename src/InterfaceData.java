import processing.core.PImage;

import java.io.File;
import java.util.*;

import static processing.core.PApplet.dist;

public class InterfaceData {
    //region declarations
    //size of the boxes used for adding extra centers (pixels in original image)
    final int BOXSIZE = 20;
    //width of the display
    final double displayWidth;
    //height of the display
    final double displayHeight;
    //scale factor of the image
    double scaleFactor;
    //width for resizing the image
    int resizeWidth;
    //height for resizing the image
    int resizeHeight;

    //the loaded palette
    Palette palette;
    //the name of the image
    String imageName;
    //the original image
    PImage originalImage;
    //the image alerted through the recoloring process
    PImage alteredImage;
    //the scaled version of the altered image
    PImage scaledAlteredImage;
    //scaled version of the segment image, used to display the segment borders to the user
    PImage scaledSegmentImage;
    //the original graph, used to reset the user inputs
    Graph originalSegmentGraph;
    //altered graph that includes added edges, etc.
    Graph alteredSegmentGraph;
    //the currently generated path
    Graph path;
    //the node at the start of the path
    Node start;
    //whether a path has been generated
    boolean pathSelected;
    //the list of centers used for selecting segments
    HashMap<Point, Integer> centers;
    //whether for not the graph is currently loaded
    boolean graphLoaded;
    //endregion
    //region setup
    public InterfaceData(int displayWidth, int displayHeight, PImage image, String imageName, PImage segmentImage, File palette){
        this.displayWidth = displayWidth;
        this.displayHeight = displayHeight;
        this.originalImage = image;
        this.alteredImage = image.copy();
        this.imageName = imageName;
        this.palette = new Palette(palette);
        setScale(image);
        scaledAlteredImage = scaleImage(originalImage).copy();
        scaledSegmentImage = scaleImage(segmentImage);
        createGraph(segmentImage);
    }
    //creates the graph based off of the images passed in, along with setting the centers and outlines of each node
    private void createGraph(PImage segmentImage) {
        AvgSegmentation originalSegmentation = new AvgSegmentation(segmentImage, originalImage);
        originalSegmentGraph = originalSegmentation.Graph();
        centers(originalSegmentation.nodeGrid);
        SetOutlines();
        alteredSegmentGraph = originalSegmentGraph.copy();
        graphLoaded = true;
    }
    //sets the scale factor along with the resize height and width, if resize height/width is 0 it resizes scaled to the non 0 value
    private void setScale(PImage image){
        if (displayWidth /image.width > displayHeight /image.height){
            scaleFactor = displayWidth / image.width;
            if (scaleFactor*image.height>displayHeight){
                scaleFactor = displayHeight / image.height;
                resizeWidth = 0;
                resizeHeight = (int) displayHeight;
            }else {
                resizeWidth = (int) displayWidth;
                resizeHeight = 0;
            }
        }else{
            scaleFactor = displayHeight /image.height;
            if (scaleFactor*image.width>displayWidth){
                scaleFactor = displayWidth /image.width;
                resizeWidth = (int) displayWidth;
                resizeHeight = 0;
            }else {
                resizeWidth = 0;
                resizeHeight = (int) displayHeight;
            }
        }
    }
    //returns a new scaled image based off of the set scale
    private PImage scaleImage(PImage image){
        PImage scaled = image.copy();
        scaled.resize(resizeWidth, resizeHeight);
        return scaled;
    }
    //sets the outline for all segments
    private void SetOutlines() {
        long start = System.currentTimeMillis();
        for (Node n: originalSegmentGraph.nodes) {
            n.scaledOutline(scaleFactor);
        }
        System.out.println("Graph: set vertices time elapsed: " + (System.currentTimeMillis()-start));
    }
    //endregion
    //region paths
    //creates the spanning tree to be used as the path
    public void generatePath(Collection<Node> selection) {
        if (selection.isEmpty()) {
            start = alteredSegmentGraph.findMidReg();
        }else {
            start = findCenterNode(selection);
        }
        path = SpanningTrees.shortestPath(alteredSegmentGraph, start);
        pathSelected = true;
    }
    //finds the center most segment in a selection
    private Node findCenterNode(Collection<Node> selection) {
        int sumx =0, sumy =0;
        for (Node n: selection){
            Point p = n.centroid;
            sumx += p.x;
            sumy += p.y;
        }
        return closestNode(sumx/selection.size(), sumy/selection.size());
    }
    //generates a partial path based off of the box selection
    public void generatePartialPath(Point p1, Point p2){
        pathSelected = false;
        inOutData data = getNodesInOut(p1,p2);
        data.out.retainAll(data.in);
        Node added = new Node();
        path = SpanningTrees.partialShortestPath(Graph.partialGraph(data.in, alteredSegmentGraph, added), data.out, added);
        pathSelected = true;
    }
    //generates a partial path based off of the shape selection
    public void generatePartialPath(List<Point> points){
        pathSelected = false;
        inOutData data = getNodesInOut(points);
        data.out.retainAll(data.in);
        Node added = new Node();
        path = SpanningTrees.partialShortestPath(Graph.partialGraph(data.in, alteredSegmentGraph, added), data.out, added);
        pathSelected = true;
    }
    //class to store what node are inside and outside a selection, all nodes both inside and outside are considered start nodes to the partial path
    static class inOutData{
        final HashSet<Node> in;
        final HashSet<Node> out;
        public inOutData(){
            in = new HashSet<>();
            out = new HashSet<>();
        }
    }

    private inOutData getNodesInOut(Point p1, Point p2){
        //normalized the points to top right and bottom left to simplify calculations
        Point top, bottom;
        int Xmax, Ymax, Xmin, Ymin;
        Xmax = Math.max(p1.x,p2.x);
        Ymax = Math.max(p1.y,p2.y);
        Xmin = Math.min(p1.x,p2.x);
        Ymin = Math.min(p1.y,p2.y);
        top = new Point(Xmax, Ymax);
        bottom = new Point(Xmin, Ymin);
        inOutData data = new inOutData();
        //if a center is contained in the box it is inside otherwise outside
        for (Map.Entry<Point, Integer> center : centers.entrySet()) {
            Point p = center.getKey();
            if (top.x > p.x && top.y > p.y && bottom.x < p.x && bottom.y < p.y){
                data.in.add(alteredSegmentGraph.nodes.get(center.getValue()));
            }else{
                data.out.add(alteredSegmentGraph.nodes.get(center.getValue()));
            }
        }
        return data;
    }
    private inOutData getNodesInOut(List<Point> points){
        //gets all points inside the shape
        HashSet<Point> interiorPoints = Algorithms.interiorPoints(Algorithms.getOutline(points), Algorithms.findInteriorPoint(points),scaledAlteredImage.width,scaledAlteredImage.height);
        inOutData data = new inOutData();
        //if center is found inside the set of interior points it is inside the shape
        for (Map.Entry<Point, Integer> center : centers.entrySet()) {
            if (interiorPoints.contains(center.getKey())){
                data.in.add(alteredSegmentGraph.nodes.get(center.getValue()));
            }else{
                data.out.add(alteredSegmentGraph.nodes.get(center.getValue()));
            }
        }
        return data;
    }
    //recolors the image based on generated path (not set up for partial path)
    public void recolor() {
        if (path == null){
            generatePath(new LinkedList<>());
        }
        alteredImage = UpdatedRozaCode.recolorImage(this);
        scaledAlteredImage = scaleImage(alteredImage);
    }
    //returns where the file is output along with its name
    public String getOutputLocation() {
        StringBuilder fileName = new StringBuilder();
        fileName.append("C:/Users/edch6/Desktop/code stuff/Research_Work/src/output/");
        fileName.append(imageName);
        fileName.append("_");
        fileName.append(palette.name);
        fileName.append("_");
        switch (path.type) {
            case segment -> fileName.append("segment");
            case shortest -> fileName.append("shortest");
            case widest -> fileName.append("widest");
            case roza -> fileName.append("roza");
            case kruskal -> fileName.append("kruskal");
        }
        fileName.append(".png");
        return fileName.toString();
    }
    //endregion
    //region centers
    //sets the centers for each segment
    void centers(int[][] nodeGrid){
        long start = System.currentTimeMillis();
        centers = new HashMap<>();
        int xOffset=0, yOffset=0;
        HashMap<Integer, LinkedList<Point>> currentBox = new HashMap<>();
        for (int i = 0; i<nodeGrid.length; i+=BOXSIZE){
            //checks if box extends past edge of image if it does set an offset to bring it to the edge of the image
            if (i+BOXSIZE>nodeGrid.length){
                xOffset = i+BOXSIZE-nodeGrid.length;
            }
            for (int j = 0; j<nodeGrid[i].length; j+=BOXSIZE){
                //checks if box extends past edge of image if it does set an offset to bring it to the edge of the image
                if (j+BOXSIZE>nodeGrid[i].length){
                    yOffset = j+BOXSIZE-nodeGrid[i].length;
                }
                //sets an entry for each segment within the box containing all points inside the box
                for (int x = i; x<BOXSIZE+i-xOffset; x++){
                    for (int y = j; y < BOXSIZE+j-yOffset; y++) {
                        LinkedList<Point> temp = currentBox.getOrDefault(nodeGrid[x][y], new LinkedList<>());
                        temp.add(new Point(x,y));
                        currentBox.put(nodeGrid[x][y], temp);
                    }
                }
                //adds the centers and resets the box and offsets
                currentBox.forEach(this::addCenter);
                currentBox.clear();
                yOffset = 0;
            }
            xOffset = 0;
        }
        //adds a center at the centriod of every segment
        for (Node n: originalSegmentGraph.nodes) {
            addCenter(n);
        }
        System.out.println("Graph: centers total time elapsed: " + (System.currentTimeMillis()-start));
    }
    //adds and sets the centroid of a segment
    void addCenter(Node n){
        int xSum =0;
        int ySum =0;
        for (Point p: n.points) {
            xSum += p.x;
            ySum += p.y;
        }
        n.centroid = new Point ((int) (((xSum/n.points.size()))*scaleFactor), (int) ((ySum/n.points.size())*scaleFactor));
        centers.put(n.centroid, n.nodeId);
    }
    //adds a center at the center of all the points contained in a box
    private void addCenter(Integer nodeID, LinkedList<Point> points) {
        int xSum =0;
        int ySum =0;
        for (Point p: points) {
            xSum += p.x;
            ySum += p.y;
        }
        Point centroid = new Point ((int) ((xSum/points.size())*scaleFactor), (int) ((ySum/points.size())*scaleFactor));
        centers.put(centroid, nodeID);
    }
    //method to get selected segments based on the box selection
    public HashSet<Node> getNodes(Point p1, Point p2){
        //normalized the points to top right and bottom left to simplify calculations
        Point top, bottom;
        int Xmax, Ymax, Xmin, Ymin;
        Xmax = Math.max(p1.x,p2.x);
        Ymax = Math.max(p1.y,p2.y);
        Xmin = Math.min(p1.x,p2.x);
        Ymin = Math.min(p1.y,p2.y);
        top = new Point(Xmax, Ymax);
        bottom = new Point(Xmin, Ymin);
        HashSet<Node> val = new HashSet<>();
        //checks every center to see if it is contained within the box
        for (Map.Entry<Point, Integer> center : centers.entrySet()) {
            Point p = center.getKey();
            if (top.x > p.x && top.y > p.y && bottom.x < p.x && bottom.y < p.y){
                val.add(alteredSegmentGraph.nodes.get(center.getValue()));
            }
        }
        return val;
    }
    //method to get selected segments based on the shape selection
    public HashSet<Node> getNodes(List<Point> points){
        //gets the set of all points contained inside the shape
        HashSet<Point> interiorPoints = Algorithms.interiorPoints(Algorithms.getOutline(points), Algorithms.findInteriorPoint(points), scaledAlteredImage.width,scaledAlteredImage.height);
        HashSet<Node> returnVal = new HashSet<>();
        //returns all nodes with a center contained within the set of interior points
        for (Map.Entry<Point, Integer> center : centers.entrySet()) {
            if (interiorPoints.contains(center.getKey())){
                returnVal.add(alteredSegmentGraph.nodes.get(center.getValue()));
            }
        }
        return returnVal;
    }
    //returns the node attached to the center closest to the given point
    public Node closestNode(int x, int y){
        double maxDist=Integer.MAX_VALUE;
        double currentDist;
        int closest=-1;
        for (Map.Entry<Point, Integer> center : centers.entrySet()) {
            currentDist = dist(x, y, center.getKey().x, center.getKey().y);
            if (currentDist < maxDist) {
                maxDist = currentDist;
                closest = center.getValue();
            }
        }
        return alteredSegmentGraph.nodes.get(closest);
    }
    //endregion
    //resets the interface to its state at launch
    public void reset(){
        alteredImage = originalImage.copy();
        scaledAlteredImage = scaleImage(originalImage);
        alteredSegmentGraph = originalSegmentGraph.copy();
    }
}