import processing.core.PImage;

import java.util.*;
// depreciated method of segmenting images, now used another image that contain all the segments
public class SegmentationModel5 {
    private static final int MIN_NODE_SIZE = 50;
    final Color[][] pixels;
    final int[][] traversed;
    final LinkedList<LinkedList<Point>> points;
    final LinkedList<HashMap<Integer, EdgeData>> neighbors;
    int currentID = 0;
    final PImage image;

    public SegmentationModel5(PImage image) {
        this.image = image;
        pixels = new Color[image.height][image.width];
        traversed = new int[image.height][image.width];
        neighbors = new LinkedList<>();
        points = new LinkedList<>();
        for (int i = 0; i < image.height; i++) {
            for (int j = 0; j < image.width; j++) {
                pixels[i][j] = new Color(image.pixels[i * image.width + j]);
                traversed[i][j] = -1;
            }
        }
    }
    /*public Graph Graph() {
        long start = System.currentTimeMillis();
        setPoints();
        System.out.println("Segmentation time elapsed: " + (System.currentTimeMillis() - start));
        return new Graph(points, neighbors, image);
    }*/
    public void setPoints() {
        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels[i].length; j++) {
                if (traversed[i][j] == -1) {
                    points.add(new LinkedList<>());
                    neighbors.add(new HashMap<>());
                    traverse(i, j, pixels[i][j], currentID);
                    int neighbor=0;
                    int currentMax=0;
                    if (points.get(currentID).size() < MIN_NODE_SIZE) {
                        for (Map.Entry<Integer,EdgeData> entry: neighbors.get(currentID).entrySet()){
                            if (entry.getValue().connections>currentMax){
                                currentMax = entry.getValue().connections;
                                neighbor = entry.getKey();
                            }
                        }
                        points.get(neighbor).addAll(points.get(currentID));
                        int finalNeighbor = neighbor;
                        neighbors.get(currentID).forEach((k, v) -> neighbors.get(finalNeighbor).put(k, EdgeData.merge(neighbors.get(currentID).get(k), v)));
                        for (Point point : points.get(currentID)) {
                            traversed[point.y][point.x] = neighbor;
                        }
                        points.remove(currentID);
                        neighbors.remove(currentID);
                    } else {
                        currentID += 1;
                    }
                }
            }
        }
    }
    public void traverse(int y, int x, Color c, int id) {
        Queue<Point> q = new LinkedList<>();
        Point p;
        q.add(new Point(x, y));
        while (!q.isEmpty()) {
            p = q.poll();
            points.get(id).add(p);

            if (checkValidity(p.y-1, p.x, c, id)){
                q.add(new Point(p.x,p.y-1));
                traversed[p.y-1][p.x] = id;
            }
            if (checkValidity(p.y+1, p.x, c, id)){
                q.add(new Point(p.x,p.y+1));
                traversed[p.y+1][p.x] = id;
            }
            if (checkValidity(p.y, p.x-1, c, id)){
                q.add(new Point(p.x-1,p.y));
                traversed[p.y][p.x-1] = id;
            }
            if (checkValidity(p.y, p.x+1, c, id)){
                q.add(new Point(p.x+1,p.y));
                traversed[p.y][p.x+1] = id;
            }
        }
    }
    private boolean checkValidity(int y, int x, Color c, int id) {
        if ((y < 0) || (y > pixels.length - 1) || (x < 0) || (x > pixels[y].length - 1)) {
            return false;
        }
        if (traversed[y][x] != -1) {
            if (traversed[y][x] != id) {
                EdgeData d = neighbors.get(id).getOrDefault(traversed[y][x], new EdgeData(1, Color.rozaColorDiff(pixels[y][x], c)));
                d.connections += 1;
                d.sumWeight += Color.colorDistance(pixels[y][x],c);
                neighbors.get(id).put(traversed[y][x], d);
            }
            return false;
        }
        return pixels[y][x].equals(c);
    }
}
