import processing.core.PImage;

import javax.imageio.ImageIO;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

import static processing.core.PApplet.dist;
//depreciated version of the interface data
public class DisplayData {
    final int BOXSIZE = 20;
    int[][] nodeGrid;
    //region fileNames
    final String imageName;
    String paletteName;
    //endregion
    //region displayInfo
    final double displayWidth;
    final double displayHeight;
    double scaleFactor;
    //endregion
    //region data
    Palette palette;
    PImage originalImage;
    PImage filteredImage;
    PImage segmentImage;
    PImage scaledOriginalImage;
    PImage scaledFilteredImage;
    PImage scaledSegmentImage;
    //endregion
    //region graph
    Graph originalSegmentGraph;
    Graph filteredSegmentGraph;
    boolean filtered = false;
    HashMap<Point, Integer> centers;
    //endregion
    //region loadChecks
    boolean imageLoaded;
    boolean graphLoaded;
    boolean paletteLoaded;
    //endregion

    public DisplayData(int displayWidth, int displayHeight, String imageName, String paletteName){
        this.displayWidth = displayWidth;
        this.displayHeight = displayHeight;
        this.imageName = imageName;
        this.paletteName = paletteName;
        readIn();
        createGraph();
    }
    public DisplayData(int displayWidth, int displayHeight, String imageName) {
        this.displayWidth = displayWidth;
        this.displayHeight = displayHeight;
        this.imageName = imageName;
        readImage();
        createGraph();
    }
    private void createGraph() {
        scaleImage();
        //originalSegmentGraph = new ColorSegmentation(segmentImage, originalImage).Graph();
        //filteredSegmentGraph = new ColorSegmentation(segmentImage, filteredImage).Graph();
        AvgSegmentation originalSegmentation = new AvgSegmentation(segmentImage, originalImage);
        originalSegmentGraph = originalSegmentation.Graph();
        nodeGrid = originalSegmentation.nodeGrid;
        filteredSegmentGraph = new AvgSegmentation(segmentImage, filteredImage).Graph();
        //roughSetVertices();
        centers();
        graphLoaded = true;
    }
    private void scaleImage() {
        scaledOriginalImage = originalImage.copy();
        scaledFilteredImage = filteredImage.copy();
        scaledSegmentImage = segmentImage.copy();
        if (displayWidth /originalImage.width > displayHeight /originalImage.height){
            scaleFactor = displayWidth / originalImage.width;
            if (scaleFactor*originalImage.height>displayHeight){
                scaleFactor = displayHeight / originalImage.height;
                scaledOriginalImage.resize(0, (int) displayHeight);
                scaledFilteredImage.resize(0, (int) displayHeight);
                scaledSegmentImage.resize(0, (int) displayHeight);
            }else {
                scaledOriginalImage.resize((int) displayWidth, 0);
                scaledFilteredImage.resize((int) displayWidth, 0);
                scaledSegmentImage.resize((int) displayWidth, 0);
            }
        }else{
            scaleFactor = displayHeight /originalImage.height;
            if (scaleFactor*originalImage.width>displayWidth){
                scaleFactor = displayWidth /originalImage.width;
                scaledOriginalImage.resize((int) displayWidth,0);
                scaledFilteredImage.resize((int) displayWidth,0);
                scaledSegmentImage.resize((int) displayWidth,0);
            }else {
                scaledOriginalImage.resize(0, (int) displayHeight);
                scaledFilteredImage.resize(0, (int) displayHeight);
                scaledSegmentImage.resize(0, (int) displayHeight);
            }
        }
    }
    /*private void roughSetVertices() {
        for (Node n: originalSegmentGraph.nodes) {
            setVertices(n);
        }
    }
    private void setVertices(Node n) {
        n.scaledVertices = scalePoints(unsortedVertices(n));
    }*/
    private LinkedList<Point> scalePoints(Collection<Point> points) {
        LinkedList<Point> returnList = new LinkedList<>();
        for (Point p: points){
            returnList.add(new Point((int)Math.floor(p.x*scaleFactor),(int)Math.floor(p.y*scaleFactor)));
        }
        return returnList;
    }
    private HashSet<Point> unsortedVertices(Node n) {
        HashSet<Point> returnVal = new HashSet<>();
        for (Point p: n.points) {
            if (!n.points.contains(new Point(p.x-1,p.y-1))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x-1,p.y))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x-1,p.y+1))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x,p.y-1))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x,p.y+1))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x+1,p.y-1))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x+1,p.y))){
                returnVal.add(p);
            }
            if (!n.points.contains(new Point(p.x+1,p.y+1))){
                returnVal.add(p);
            }
        }
        return returnVal;
    }
    void readIn(){
        readImage();
        readPalette();
    }
    void readImage() {
        try{
            originalImage = new PImage(ImageIO.read(new File("./src/images/originals/" + imageName + ".jpg")));
            filteredImage = new PImage(ImageIO.read(new File("./src/images/filtered/" + imageName + ".png")));
            segmentImage = new PImage(ImageIO.read(new File("./src/images/segmentation/" + imageName + ".png")));
            imageLoaded = true;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    void readPalette(){
        try {
            Scanner scanner = new Scanner(new File("./src/palettes/" + paletteName + ".txt"));
            LinkedList<String> lines = new LinkedList<>();
            while (scanner.hasNextLine()){
                String line = scanner.nextLine();
                lines.add(line);
            }
            palette = new Palette(lines.toArray(String[]::new),"unNamed");
            paletteLoaded=true;
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
    public void reset(){
        readImage();
        createGraph();
    }
    public Graph segmentGraph() {
        if (filtered){
            return filteredSegmentGraph;
        }else{
            return originalSegmentGraph;
        }
    }
    public Graph widestPath(Node selectedNode) {
        graphLoaded = false;
        Graph widestPath;
        if (filtered){
            widestPath = SpanningTrees.widestPath(filteredSegmentGraph,selectedNode);
        }else {
            widestPath =SpanningTrees.widestPath(originalSegmentGraph, selectedNode);
        }
        graphLoaded = true;
        return widestPath;
    }
    public Graph shortestPath(Node selectedNode) {
        graphLoaded = false;
        Graph shortestPath;
        if (filtered){
            shortestPath = SpanningTrees.shortestPath(filteredSegmentGraph,selectedNode);
        }else {
            shortestPath =SpanningTrees.shortestPath(originalSegmentGraph, selectedNode);
        }
        graphLoaded = true;
        return shortestPath;
    }
    public Graph kruskal(){
        graphLoaded = false;
        Graph kruskal;
        if (filtered){
            kruskal = SpanningTrees.kruskal(filteredSegmentGraph);
        }else {
            kruskal =SpanningTrees.kruskal(originalSegmentGraph);
        }
        graphLoaded = true;
        return kruskal;
    }
    public Graph rozaGraph(Node selectedNode){
        graphLoaded = false;
        Graph rozaGraph;
        if (filtered){
            rozaGraph = SpanningTrees.rozaGraph(filteredSegmentGraph,selectedNode);
        }else {
            rozaGraph =SpanningTrees.rozaGraph(originalSegmentGraph, selectedNode);
        }
        graphLoaded = true;
        return rozaGraph;
    }
    public void updateImage() {
        scaleImage();
    }
    public void setPalette(String paletteName){
        this.paletteName = paletteName;
        paletteLoaded = false;
        readPalette();
    }
    void centers(){
        long start = System.currentTimeMillis();
        centers = new HashMap<>();
        int xOffset=0, yOffset=0;
        HashMap<Integer, LinkedList<Point>> currentBox = new HashMap<>();
        for (int i = 0; i<nodeGrid.length; i+=BOXSIZE){
            if (i+BOXSIZE>nodeGrid.length){
                xOffset = i+BOXSIZE-nodeGrid.length;
            }
            for (int j = 0; j<nodeGrid[i].length; j+=BOXSIZE){
                if (j+BOXSIZE>nodeGrid[i].length){
                    yOffset = j+BOXSIZE-nodeGrid[i].length;
                }
                for (int x = i; x<BOXSIZE+i-xOffset; x++){
                    for (int y = j; y < BOXSIZE+j-yOffset; y++) {
                        LinkedList<Point> temp = currentBox.getOrDefault(nodeGrid[x][y], new LinkedList<>());
                        temp.add(new Point(x,y));
                        currentBox.put(nodeGrid[x][y], temp);
                    }
                }
                currentBox.forEach(this::addCenter);
                currentBox.clear();
                yOffset = 0;
            }
            xOffset = 0;
        }
        for (Node n: originalSegmentGraph.nodes) {
            addCenter(n);
        }
        System.out.println("Graph: centers total time elapsed: " + (System.currentTimeMillis()-start));
    }
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
    private void addCenter(Integer nodeID, LinkedList<Point> points) {
        int xSum =0;
        int ySum =0;
        for (Point p: points) {
            xSum += p.x;
            ySum += p.y;
        }
        Point centroid = new Point ((int) (((xSum/points.size()))*scaleFactor), (int) ((ySum/points.size())*scaleFactor));
        centers.put(centroid, nodeID);
    }
    public Node closestNode(int mouseX, int mouseY){
        double maxDist=Integer.MAX_VALUE;
        double currentDist;
        int closest=-1;
        for (Map.Entry<Point, Integer> center: centers.entrySet()) {
            currentDist = dist(mouseX, mouseY, center.getKey().x, center.getKey().y); /// Math.log10(Math.sqrt(originalSegmentGraph.nodes.get(center.getValue()).points.size()));
            if (currentDist<maxDist) {
                maxDist = currentDist;
                closest = center.getValue();
            }
        }
        return originalSegmentGraph.nodes.get(closest);
    }
}
