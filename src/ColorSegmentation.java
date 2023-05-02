import processing.core.PImage;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
//creates edges based off of the average color distance of every connection between pair of two pixels
public class ColorSegmentation {
    private final PImage originalImage;
    private final PImage segmentImage;
    //representation of the image as map of what pixel has been assigned to what node
    final int[][] nodeGrid;
    //hashmap storing a point holding the id of the two nodes and the edge data for that edge
    final HashMap<Point, EdgeData> neighbors;
    //list storing hashset of points for each segment
    final LinkedList<HashSet<Point>> points;

    public ColorSegmentation(PImage segmentImage, PImage originalImage) {
        this.segmentImage = segmentImage;
        this.originalImage = originalImage;
        nodeGrid = new int[segmentImage.width][segmentImage.height];
        neighbors = new HashMap<>();
        points = new LinkedList<>();
        for (int i = 0; i < segmentImage.width; i++) {
            for (int j = 0; j < segmentImage.height; j++) {
                nodeGrid[i][j] = -1;
            }
        }
    }
    //returns the graph based on the segmentation
    public Graph graph() {
        long start = System.currentTimeMillis();
        setPoints();
        System.out.println("Segmentation time elapsed: " + (System.currentTimeMillis() - start));
        return new Graph(points, neighbors, originalImage);
    }
    //sets the points based on the segmentation
    void setPoints(){
        int currentID = 0;
        for (int i = 0; i < segmentImage.width; i++) {
            for (int j = 0; j < segmentImage.height; j++) {
                //searches for the next unassigned pixel, when found floodfills around it to find the rest of the pixels that are a part of that segment
                if (nodeGrid[i][j] == -1) {
                    points.add(new HashSet<>());
                    traverse(new Point(i,j), segmentImage.pixels[j * segmentImage.width + i], currentID);
                    currentID++;
                }
            }
        }
    }
    //iterative version of flood fill
    void traverse(Point p, int color, int id){
        LinkedList<Point> q = new LinkedList<>();
        q.add(p);
        nodeGrid[p.x][p.y] = id;
        while (!q.isEmpty()) {
            p = q.poll();
            points.get(id).add(p);
            for (Point point: Algorithms.surroundingPoints(p)){
                if (checkValidity(p, point, color, id)){
                    q.add(point);
                    nodeGrid[point.x][point.y] = id;
                }
            }
        }
    }
    boolean checkValidity(Point p, Point point, int c, int id) {
        //bounds check
        if ((point.y < 0) || (point.y > segmentImage.height - 1) || (point.x < 0) || (point.x > segmentImage.width - 1)) {
            return false;
        }
        //if traversed
        if (nodeGrid[point.x][point.y] != -1) {
            // if not part of this segment add the connection
            if (nodeGrid[point.x][point.y] != id) {
                Point edge = new Point(id, nodeGrid[point.x][point.y]);
                EdgeData d = new EdgeData(1, Color.rozaColorDiff(originalImage.pixels[p.y*originalImage.width+p.x],originalImage.pixels[point.y*originalImage.width+point.x]));
                neighbors.merge(edge, d, EdgeData::merge);
            }
            return false;
        }
        //if share color, in the segment image all pixels of a segment share a color
        return Color.rozaColorDiff(c, segmentImage.pixels[point.y * segmentImage.width + point.x]) == 0;
    }
}
