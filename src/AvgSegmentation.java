import processing.core.PImage;

import java.util.HashSet;
import java.util.LinkedList;
//creates edges based off of the color distance of the average color of the two segments
public class AvgSegmentation {
    private final PImage originalImage;
    private final PImage segmentImage;
    //representation of the image as map of what pixel has been assigned to what node
    final int[][] nodeGrid;
    //hashset storing pairs of node ids representing edges
    final HashSet<Point> neighbors;
    //list storing hashset of points for each segment
    final LinkedList<HashSet<Point>> points;

    public AvgSegmentation(PImage segmentImage, PImage originalImage) {
        this.segmentImage = segmentImage;
        this.originalImage = originalImage;
        nodeGrid = new int[segmentImage.width][segmentImage.height];
        neighbors = new HashSet<>();
        points = new LinkedList<>();
        for (int i = 0; i < segmentImage.width; i++) {
            for (int j = 0; j < segmentImage.height; j++) {
                nodeGrid[i][j] = -1;
            }
        }
    }
    //returns the graph based on the segmentation
    public Graph Graph() {
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
                if (nodeGrid[i][j] == -1) {
                    points.add(new HashSet<>());
                    traverse(new Point(i, j), segmentImage.pixels[j * segmentImage.width + i], currentID);
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
            for (Point point : Algorithms.surroundingPoints(p)) {
                if (checkValidity(point, color, id)) {
                    q.add(point);
                    nodeGrid[point.x][point.y] = id;
                }
            }
        }
    }
    boolean checkValidity(Point p, int c, int id) {
        //bounds check
        if ((p.y < 0) || (p.y > segmentImage.height - 1) || (p.x < 0) || (p.x > segmentImage.width - 1)) {
            return false;
        }
        //if traversed
        if (nodeGrid[p.x][p.y] != -1) {
            //if not part of this segment add the connection
            if (nodeGrid[p.x][p.y] != id) {
                neighbors.add(new Point(id, nodeGrid[p.x][p.y]));
            }
            return false;
        }
        //if share color, in the segment image all pixels of a segment share a color
        return Color.rozaColorDiff(c, segmentImage.pixels[p.y * segmentImage.width + p.x]) == 0;
    }
}