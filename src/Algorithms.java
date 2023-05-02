import processing.core.PShape;
import processing.core.PVector;

import java.util.*;

public class Algorithms {
    //http://www.jeffreythompson.org/collision-detection/poly-point.php
    //algorithm for checking collision
    public static boolean collision(int px, int py, PShape shape) {
        boolean collision = false;
        int next;

        for (int current = 0; current< shape.getVertexCount(); current++) {
            next = current+1;
            if (next == shape.getVertexCount()) next = 0;

            PVector vc = shape.getVertex(current);
            PVector vn = shape.getVertex(next);

            if (((vc.y >= py && vn.y < py) || (vc.y < py && vn.y >= py)) &&
                    (px < (vn.x-vc.x)*(py-vc.y) / (vn.y-vc.y)+vc.x)) {
                collision = !collision;
            }
        }
        return collision;
    }
    public static boolean collision(int px, int py, List<Point> shape) {
        boolean collision = false;
        int next;

        for (int current = 0; current< shape.size(); current++) {
            next = current+1;
            if (next == shape.size()) next = 0;

            Point vc = shape.get(current);
            Point vn = shape.get(next);

            if (((vc.y >= py && vn.y < py) || (vc.y < py && vn.y >= py)) &&
                    (px < (vn.x-vc.x)*(py-vc.y)*1.0 / (vn.y-vc.y)+vc.x)) {
                collision = !collision;
            }
        }
        return collision;
    }
    //returns all points surrounding the passed point
    static Point[] surroundingPoints(Point p){
        return new Point[]{
                new Point(p.x-1,p.y-1),
                new Point(p.x-1,p.y),
                new Point(p.x-1,p.y+1),
                new Point(p.x, p.y-1),
                new Point(p.x,p.y+1),
                new Point(p.x+1,p.y-1),
                new Point(p.x+1,p.y),
                new Point(p.x+1,p.y)};
    }
    //returns a hashset containing all point on the exterior of a shape
    static HashSet<Point> getOutline(List<Point> points){
        HashSet<Point> val = new HashSet<>();
        for (int current = 0; current< points.size(); current++) {
            int next = current + 1;
            if (next == points.size()) next = 0;

            Point p1 = points.get(current);
            Point p2 = points.get(next);
            val.addAll(plotLine(p1.x,p1.y,p2.x,p2.y));
        }
        return val;
    }
    //https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm#All_cases
    //algorithm for getting all points on a line
    static List<Point> plotLine(int x0, int y0, int x1, int y1) {
        if (Math.abs(y1 - y0) < Math.abs(x1 - x0)) {
            if (x0 > x1) {
                return plotLineLow(x1, y1, x0, y0);
            } else {
                return plotLineLow(x0, y0, x1, y1);
            }
        } else{
            if (y0 > y1) {
                return plotLineHigh(x1, y1, x0, y0);
            }else {
                return plotLineHigh(x0, y0, x1, y1);
            }
        }
    }
    static List<Point> plotLineLow(int x0, int y0, int x1, int y1) {
        List<Point> points = new LinkedList<>();
        int dx = x1 - x0;
        int dy = y1 - y0;
        int yi = 1;
        if (dy < 0) {
            yi = -1;
            dy = -dy;
        }
        int D = (2 * dy) - dx;
        int y = y0;

        for (int x = x0; x <= x1; x++) {
            points.add(new Point(x, y));
            if (D > 0) {
                y = y + yi;
                D = D + (2 * (dy - dx));
            } else{
                D = D + 2 * dy;
            }
        }
        return points;
    }
    static List<Point> plotLineHigh(int x0, int y0, int x1, int y1) {
        List<Point> points = new LinkedList<>();
        int dx = x1 - x0;
        int dy = y1 - y0;
        int xi = 1;
        if (dx < 0){
            xi = -1;
            dx = -dx;
        }
        int D = (2 * dx) - dy;
        int x = x0;
        for (int y = y0; y <= y1; y++) {
            points.add(new Point(x, y));
            if (D > 0) {
                x = x + xi;
                D = D + 2 * (dx - dy);
            } else {
                D = D + 2 * dx;
            }
        }
        return points;
    }
    //returns a random point inside a shape
    public static Point findInteriorPoint(List<Point> points) {
        int xmax = points.stream().mapToInt(p -> p.x).max().orElse(-1);
        int ymax = points.stream().mapToInt(p -> p.y).max().orElse(-1);
        int xmin = points.stream().mapToInt(p -> p.x).min().orElse(-1);
        int ymin = points.stream().mapToInt(p -> p.y).min().orElse(-1);
        //error check
        if (xmax <= xmin || ymax <= ymin || xmin == -1 || ymin == -1){
            return new Point (-1,-1);
        }
        Point p;
        Random rand = new Random();
        int count =0;
        do{
            p = new Point (xmin + (int) Math.floor((xmax-xmin) * rand.nextDouble()), ymin + (int) Math.floor((ymax-ymin) * rand.nextDouble()));
            count++;
        }while (count < 100 && !collision(p.x,p.y,points));
        //error condition if point can't be found
        if (count >= 100){
            return new Point(-1,-1);
        }
        return p;
    }
    //finds all points inside a shape using floodfill
    public static HashSet<Point> interiorPoints(HashSet<Point> outline, Point startPoint, int width, int height) {
        HashSet<Point> interiorPoints = new HashSet<>();
        Queue<Point> queue = new LinkedList<>();
        queue.add(startPoint);
        while (!queue.isEmpty()) {
            Point point = queue.poll();
            Point[] checks = {new Point(point.x+1,point.y), new Point(point.x-1,point.y), new Point(point.x,point.y+1), new Point(point.x,point.y-1)};
            for (Point p: checks){
                //error condition
                if (p.x < 0 || p.y < 0 || p.x > width || p.y > height){
                    return new HashSet<>();
                }
                if (!outline.contains(p) && !interiorPoints.contains(p)){
                    queue.add(p);
                    interiorPoints.add(p);
                }
            }
        }
        return interiorPoints;
    }

    //original port of roza's code, kept to make sure nothing changed too much from the original, no longer works due to refactorings
    //changed to only handle one palette
    /*public static void makeGraphonRegions(PImage img, Graph graph, Palette palette){
        PImage draw = img.copy();
        RecoloringByBottleNeck(graph, draw, palette);
        draw.save("C:/Users/edch6/Desktop/code stuff/Research_Work/src/test.png");
    }

     static void RecoloringByBottleNeck(Graph graph, PImage screen, Palette palette) {
        ArrayList<Integer> regsizes = new ArrayList<>(graph.size());
        for (int i =0; i< graph.size(); i++){
            regsizes.add(graph.nodes[i].points.size());
        }
        double avesize = CalcAverage(regsizes);

        double[][] table = new double[palette.size()][palette.size()];
        for (int i = 0; i < palette.size(); i++){
            for (int j = 0; j < palette.size(); j++){
                table[i][j] = palette.get(i).rozaColorDiff(palette.get(j));
            }
        }
        ArrayList<Integer> lookUp_table = new ArrayList<>();
        HistogramMatching(graph.getEdgeWeights(), table);
        double midweight = getMidWeight(graph.getEdgeWeights());
        int mid_reg = graph.findMidReg();

        Graph widestGraph =DijkstraAlgorithms.widestPathColoring(graph, graph.nodes[mid_reg], midweight, avesize);
        for (Node node: widestGraph.nodes){
            node.colored = false;
        }
        AssignColors(widestGraph, widestGraph.nodes[mid_reg], palette, table, lookUp_table, graph.nodes[mid_reg]);
        int ind=0;
        ArrayList<Integer> bin = new ArrayList<>();
        for (int i = 0; i < palette.size(); i++) {
            bin.add(0);
        }
        int index;
        int value;
        for (int i = 0; i < graph.size(); i++) {
            Color pcol = graph.nodes[i].color;
            for (int p = 0; p < palette.size(); p++) {
                if (palette.get(p) == pcol)
                    ind = p;
            }
            for (Point p: graph.nodes[i].points)
            {
                screen.set(p.x,p.y,pcol.processingForm());
                index = ind;
                value = bin.get(ind);
                bin.set(index, ++value);
            }
        }
    }
    static double CalcAverage(int[] scores) {
        double avesize=0.0;
        for (Integer score : scores) {
            avesize += score;
        }
        return (avesize / scores.size());
    }
    static void HistogramMatching(Double[] edgeWeights, double[][] table){
        Arrays.sort(edgeWeights);
        ArrayList<Double> rArr = new ArrayList<>(edgeWeights.length);
        rArr.addAll(Arrays.asList(edgeWeights));

        int numbins = table.length * table[0].length;
        ArrayList<Double> tArr = new ArrayList<>(numbins);
        for (double[] doubles : table) {
            for (int j = 0; j < table[0].length; j++) {
                tArr.add(doubles[j]);
            }
        }

        int bsize = 500;
        ArrayList<Integer> rfreq = new ArrayList<>(bsize);
        ArrayList<Integer> tfreq = new ArrayList<>(bsize);
        for (int i =0; i<bsize; i++){
            rfreq.add(0);
            tfreq.add(0);
        }

        int index;
        int value;
        for (Double aDouble : rArr) {
            index = (int) (aDouble.doubleValue());
            value = rfreq.get(index);
            rfreq.set(index, ++value);
        }
        for (Double aDouble : tArr) {
            index = (int) (aDouble.doubleValue());
            value = tfreq.get(index);
            tfreq.set(index, ++value);
        }
        ArrayList<Double> normalized_rcdf = new ArrayList<>(bsize);
        ArrayList<Double> normalized_tcdf = new ArrayList<>(bsize);
        for (int i = 0; i < bsize; i++) {
            normalized_rcdf.add(0.0);
            normalized_tcdf.add(0.0);
        }
        calculate_cdf(rfreq, normalized_rcdf);
        calculate_cdf(tfreq, normalized_tcdf);
        ArrayList<Integer> lookup_table = lookupTable(normalized_rcdf, normalized_tcdf);

        ArrayList<Integer> Mappedfreq = new ArrayList<>(bsize);
        for (int i = 0; i < bsize; i++) {
            Mappedfreq.add(0);
        }
        int id;
        for (Double edgeWeight : edgeWeights) {
            id = (int) Math.round(edgeWeight);
            index = lookup_table.get(id);
            value = Mappedfreq.get(index);
            Mappedfreq.set(index, ++value);
        }

        //calculate cdf after mapping
        ArrayList<Double> normalized_Mappedtcdf = new ArrayList<>();
        for (int i =0; i< Mappedfreq.size(); i++){
            normalized_Mappedtcdf.add(0.0);
        }
        calculate_cdf(Mappedfreq, normalized_Mappedtcdf);
    }
    static double getMidWeight(Double[] weights){
        ArrayList<Double> nweights = new ArrayList<>(Arrays.asList(weights));
        double midweight;
        if (nweights.size() <=1){
            midweight =0;
        }else{
            double maxw = Collections.max(nweights);
            double minw = Collections.min(nweights);
            midweight = (maxw + minw) / 2.0;
        }
        return midweight;
    }

    static void calculate_cdf(ArrayList<Integer> hist, ArrayList<Double> normalized_cdf){
        double sum = 0;
        double sH = 0;
        for (Integer i: hist){
            sH += i;
        }
        ArrayList<Double> cdf = new ArrayList<>(hist.size());
        for (Integer i: hist) {
            sum += i/sH;
            cdf.add(sum);
        }
        for (int i = 0; i < hist.size(); i++) {
            normalized_cdf.set(i, cdf.get(i));
        }
    }
    static ArrayList<Integer> lookupTable(ArrayList<Double> normalized_rcdf, ArrayList<Double> normalized_tcdf){
        int size = normalized_rcdf.size();
        ArrayList<Integer> lookup_table = new ArrayList<>(size);
        for (int i = 0; i < size;i++) {
            lookup_table.add(0);
        }

        int j = 0;
        for (int i = 0; i <size;i++) {
            while (j+1 < normalized_tcdf.size() && normalized_tcdf.get(j) < normalized_rcdf.get(i)) {
                j++;
            }
            lookup_table.set(i, j);
        }
        return lookup_table;
    }
    static void AssignColors(Graph graph, Node root, Palette palette, double[][] table, ArrayList<Integer> lookUp_table, Node node){
        int index = ReturnClosestPaletteIndex1(graph, root, palette, table, lookUp_table, node);

        node.color = palette.get(index);
        node.assignedColors.add(palette.get(index));
        node.colored = true;
        for (int i = 0; i < node.children.size(); i++) {
            if (!node.children.get(i).colored){
                AssignColors(graph, root, palette, table, lookUp_table, node.children.get(i));
            }
        }
    }

    private static int ReturnClosestPaletteIndex1(Graph graph, Node root, Palette palette, double[][] table, ArrayList<Integer> lookUp_table, Node node) {
        double mindif = 1000;
        double difp;
        int index = -1;
        boolean signr;
        boolean signp;
        double g1, g2, p1, p2;
        int pidex = -1;
        double colorDif;
        Color col1, col2;
        double minPaleDif = Double.MAX_VALUE, MaxPaleDif = Double.MIN_VALUE, maxEdgeDif = Double.MIN_VALUE;
        for (double[] arr : table) {
            for (double d : arr) {
                if (minPaleDif > d) {
                    minPaleDif = d;
                } else if (MaxPaleDif < d) {
                    MaxPaleDif = d;
                }
            }
        }
        for (double d : graph.getEdgeWeights()) {
            if (d > maxEdgeDif) {
                maxEdgeDif = d;
            }
        }
        col2 = graph.avgColor(node);
        double diff;
        if (node.parent == null || node.nodeId == root.nodeId) {
            for (Edge e : graph.edges.get(node.nodeId)) {
                if (e.nodeB.colored) {
                    LinkedList<Color> colors = palette.colors;
                    for (int j = 0; j < colors.size(); j++) {
                        Color c = colors.get(j);
                        if (e.nodeB.color == c) {
                            pidex = j;
                        }
                    }
                }
            }
            if (pidex == -1) {
                LinkedList<Color> colors = palette.colors;
                for (int i = 0; i < colors.size(); i++) {
                    Color c = colors.get(i);
                    diff = col2.rozaColorDiff(c);
                    if (diff < mindif) {
                        mindif = diff;
                        index = i;
                    }
                }
            }
        }
        double regiondif;
        int id;
        double dif;
        Color palette_j;
        if (node.parent != null){
            col1 = graph.avgColor(node.parent);
            for (int i =0; i < palette.size(); i++){
                if (node.parent.color == palette.get(i)){
                    pidex = i;
                    break;
                }
            }
            p1 = (0.2126 *palette.get(pidex).r) + (0.7152 *palette.get(pidex).g) + (0.0722 * palette.get(pidex).b);
            g1 = (0.2126 *col1.b) + (0.7152 *col1.b) + (0.0722 * col1.b);	// return luminance
            g2 = (0.2126 *col2.b) + (0.7152 *col2.b) + (0.0722 * col2.b);	// return luminance
            regiondif = col1.rozaColorDiff(col2);
            if (!lookUp_table.isEmpty()){
                id = (int) Math.round(regiondif);
                colorDif = lookUp_table.get(id);
            }else {
                colorDif = regiondif;
            }
            signr = g1 >= g2;
            mindif = 1000;
            for (int j =0; j< palette.size(); j++) {
                palette_j = palette.get(j);
                difp = palette.get(pidex).rozaColorDiff(palette_j);
                p2 = (((0.2126 *palette_j.r) + (0.7152 *palette_j.g) + (0.0722 * palette_j.b)));
                signp = p1>=p2;
                dif = Math.abs(colorDif - difp);
                if (signr == signp) {
                    if (dif <mindif) {
                        mindif = dif;
                        index = j;
                    }
                }
            }
            if (index == -1) {
                mindif = 1000;
                for (int j = 0; j < palette.size(); j++) {
                    difp = palette.get(j).rozaColorDiff(palette.get(pidex));
                    dif = Math.abs (colorDif - difp);
                    if (dif < mindif) {
                        mindif = dif;
                        index = j;
                    }
                }
            }
        }
        return index;
    }*/
}
