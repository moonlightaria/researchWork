import processing.core.PApplet;
import java.util.LinkedList;
//class to test if segmentation model works properly, set up for the old version the display data
public class SegmentationModelTester extends PApplet {
    static DisplayData data;
    static LinkedList<Point> missing;
    static LinkedList<Point> duplicate;
    static LinkedList<Node> edgeless;
    public static void main(String[] args) {
        PApplet.main("SegmentationModelTester");
    }
    public void settings(){
        size(1800,900);
    }
    public void setup(){
        long start = System.currentTimeMillis();
        data = new DisplayData(1400,900,"bike1");
        System.out.println("ModelCreation time elapsed: " + (System.currentTimeMillis() - start));
        verify(data.originalSegmentGraph);
    }
    //draws where the errors are on the image
    public void draw(){
        //image(data.originalImage,0,0);
        stroke(0);
        for (Point p: missing){
            point(p.x,p.y);
        }
        stroke(255);
        for (Point p: duplicate){
            point(p.x,p.y);
        }
        stroke(color(255,255,255));
        for (Node n: edgeless){
            for (Point p: n.points){
                point(p.x,p.y);
            }
        }
    }
    //traversed the graph comparing it to the image to find any errors
    public static void verify(Graph g) {
        System.out.println("num of nodes: " + g.size());
        System.out.println("num of edge: " + g.edgeCount);
        boolean[][] traversed = new boolean[data.originalImage.height][data.originalImage.width];
        for (int i =0; i<data.originalImage.height; i++) {
            for (int j = 0; j < data.originalImage.width; j++) {
                traversed[i][j] = false;
            }
        }
        duplicate = new LinkedList<>();
        edgeless = new LinkedList<>();
        for (Node n: g.nodes) {
            //check if segments are merged if merging is being used
            /*if (n.points.size() < 5) {
                System.out.println(n.nodeId + " not merged");
            }*/
            //checks if the pixel exists on the image and if it has already been added to a segment
            for (Point p : n.points) {
                if (traversed[p.y][p.x]) {
                    System.out.println(p.x + " ," + p.y + " |duplicate data");
                    duplicate.add(p);
                } else {
                    traversed[p.y][p.x] = true;
                }
            }
            //checks if a segment has no edges
            if (g.edges.get(n.nodeId).size() <1){
                System.out.println(n + " | no edges");
                for (Point p: n.points){
                    System.out.println(p.x + " ," + p.y + " |edgeless point");
                }
                edgeless.add(n);
            }
        }
        //checks if any pixels hasn't been added to a segment
        missing = new LinkedList<>();
        for(int i =0; i<data.originalImage.height; i++) {
            for (int j = 0; j < data.originalImage.width; j++) {
                if (!traversed[i][j]) {
                    missing.add(new Point(i,j));
                    System.out.println(i + " ," + j + " |missing data");
                }
            }
        }
    }
}
