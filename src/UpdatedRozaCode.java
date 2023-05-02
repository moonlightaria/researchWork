import processing.core.PImage;

import java.util.HashMap;
import java.util.Map;

//altered version of roza's code to handle various user inputs
public class UpdatedRozaCode {
    static double minPaleDif = Double.MAX_VALUE;
    static double MaxPaleDif = Double.MIN_VALUE;
    static double[][] table;
    static Palette palette;
    public static PImage recolorImage(InterfaceData data){
        long start = System.currentTimeMillis();
        PImage draw = data.alteredImage.copy();
        palette = data.palette;
        data.start.colored = true;
        traverseGraph(data.path, data.start);
        //unmarks all nodes as traversed and assigns the final color to any node that shares a color with on assigned a final color
        for (Node n : data.path.nodes) {
            n.colored = false;
            if (n.finalColor != null) {
                for (Node node: n.shareColor){
                    node.finalColor = n.finalColor;
                }
            }
        }
        RecoloringByBottleNeck(data.path, data.start, draw);
        //resets all the fields used in recoloring
        for (Node n: data.path.nodes){
            n.color = null;
            n.colored = false;
            n.parent = null;
        }
        System.out.println("Recoloring: total time elapsed: " + (System.currentTimeMillis()-start));
        return draw;
    }
    //traverses the graph assigning parent to each node (work around for old implementation, maybe could be moved back into spanning tree creation)
    private static void traverseGraph(Graph graph, Node current) {
        for (Edge e: graph.edges.get(current.nodeId)){
            if (!e.nodeB.colored) {
                e.nodeB.parent = e.nodeA;
                e.nodeB.colored = true;
                traverseGraph(graph, e.nodeB);
            }
        }
    }
    //main function of the recoloring processes
    static void RecoloringByBottleNeck(Graph spanningTree, Node source, PImage image) {
        //calculates the table (difference between each color in the palette, and gets the min and max difference
        table = new double[palette.size()][palette.size()];
        for (int i = 0; i < palette.size(); i++){
            for (int j = 0; j < palette.size(); j++){
                table[i][j] = Color.rozaColorDiff(palette.get(i), palette.get(j));
                if (minPaleDif > table[i][j]) {
                    minPaleDif = table[i][j];
                } else if (MaxPaleDif < table[i][j]) {
                    MaxPaleDif = table[i][j];
                }
            }
        }
        int[] lookUp_table = HistogramMatching(spanningTree.getEdgeWeights(), table);
        AssignColors(spanningTree, source, palette, lookUp_table, source);

        //recolors the image based off of the assigned colors
        for (Node node: spanningTree.nodes) {
            if (!node.dontColor) {
                for (Point p : node.points) {
                    image.set(p.x, p.y, node.color.processingForm());
                }
            }
        }
    }
    //get the frequency of each value and calculates cdf to pass to create the lookup table
    static int[] HistogramMatching(double[] edgeWeights, double[][] table){
        int bsize = 500;
        int[] rfreq = new int[bsize];
        int[] tfreq = new int[bsize];

        //gets frequency of each value for both edge weights and the table
        for (double d : edgeWeights) {
            rfreq[(int) d]++;
        }
        for (double[] arr: table) {
            for (double d: arr){
                tfreq[(int) d]++;
            }
        }
        double[] normalized_rcdf = new double[bsize];
        double[] normalized_tcdf = new double[bsize];

        calculate_cdf(rfreq, normalized_rcdf);
        calculate_cdf(tfreq, normalized_tcdf);
        return lookupTable(normalized_rcdf, normalized_tcdf);
    }
    //assertion normalized_cdf is bigger or equal in size to hist
    static void calculate_cdf(int[] hist, double[] normalized_cdf){
        double sum = 0;
        double sumHist = 0;
        for (Integer i: hist){
            sumHist += i;
        }
        for (int i = 0; i < hist.length; i++) {
            sum += hist[i] / sumHist;
            normalized_cdf[i] = sum;
        }
    }
    //creates the lookup table based on the cdf of the edge weights and palette differences
    static int[] lookupTable(double[] normalized_rcdf, double[] normalized_tcdf){
        int size = normalized_rcdf.length;
        int[] lookup_table = new int[size];
        int j = 0;

        for (int i = 0; i <size;i++) {
            while (j+1 < normalized_tcdf.length && normalized_tcdf[j] < normalized_rcdf[i]) {
                j++;
            }
            lookup_table[i] =j;
        }
        return lookup_table;
    }
    //assigns colors to segments as it traversed the spanning tree
    static void AssignColors(Graph graph, Node root, Palette palette, int[] lookUp_table, Node node){
        //checks if a node has been assigned a color, otherwise assigns a color
        if (node.finalColor != null) {
            node.color = node.finalColor;
        }else{
            node.color = returnColor(graph, root, palette, lookUp_table, node);
        }
        //assigns the color it was recolored to all nodes that share a color with it
        for (Node n: node.shareColor){
            n.finalColor = node.color;
        }
        //recursively navigate to the next node unrecolored node in the graph
        node.colored = true;
        for (Edge e: graph.edges.get(node.nodeId)) {
            if (!e.nodeB.colored){
                AssignColors(graph, root, palette, lookUp_table, e.nodeB);
            }
        }
    }
    //assigns a color to a given node
    private static Color returnColor(Graph graph, Node root, Palette palette, int[] lookUp_table, Node node) {
        HashMap<Color, Double> distanceMap = new HashMap<>(palette.size());
        double g1, g2, p1, p2;
        double colorDif, diff, difp;
        Color col1, col2, pColor;
        col2 = node.avgColor;
        //if no parent exists or is root node, if a neighboring node has been color recolor it to the same color as that neighbor
        if (node.parent == null || node.nodeId == root.nodeId) {
            for (Edge e : graph.edges.get(node.nodeId)) {
                if (e.nodeB.colored) {
                    return e.nodeB.color;
                }
            }
            //otherwise get the difference between the nodes average color and the palette colors
            for (Color c : palette.colors) {
                diff = Color.rozaColorDiff(col2, c) * bias(node,c);
                distanceMap.put(c,diff);
            }
            return validColor(node, distanceMap);
        }
        //not exactly sure the significance of most of these calculations, does some luminance stuff and if one is both are brighter or dimer recolor the node to that
        col1 = node.parent.avgColor;
        pColor = node.parent.color;
        p1 = (0.2126 * pColor.r) + (0.7152 * pColor.g) + (0.0722 * pColor.b);
        g1 = (0.2126 *col1.r) + (0.7152 *col1.g) + (0.0722 * col1.b);	// return luminance
        g2 = (0.2126 *col2.r) + (0.7152 *col2.g) + (0.0722 * col2.b);	// return luminance
        colorDif = lookUp_table[(int) Math.round(Color.rozaColorDiff(col1, col2))];
        boolean signr = g1 >= g2;
        for (Color c: palette.colors) {
            difp = Color.rozaColorDiff(pColor, c);
            p2 = (0.2126 *c.r) + (0.7152 *c.g) + (0.0722 * c.b);
            boolean signp = p1>=p2;
            diff = Math.abs(colorDif - difp) * bias(node,c);
            if (signr == signp) {
                distanceMap.put(c,diff);
            }
        }
        if (!distanceMap.isEmpty()) {
            return validColor(node, distanceMap);
        }
        //if the above check failed to add any colors, run the above distance calculation on all colors
        for (Color c: palette.colors) {
            difp = Color.rozaColorDiff(pColor, c);
            diff = Math.abs(colorDif - difp) * bias(node,c);
            distanceMap.put(c,diff);
        }
        return validColor(node, distanceMap);
    }
    //a function to return the bias towards or away from a node
    private static double bias(Node node, Color c) {
        double biasA = 2;
        double biasB = 10;
        if (node.biasTowardsColor.contains(c) || node.biasTowardsSegment.stream().anyMatch(n -> n.finalColor == c || n.color == c)){
            return biasA/biasB;
        }
        if (node.biasFromColor.contains(c) || node.biasFromSegment.stream().anyMatch(n -> n.finalColor == c || n.color == c)){
            return biasB/biasA;
        }
        return 1;
    }

    //function to return the color based on the given distances
    private static Color validColor(Node node, HashMap<Color, Double> distanceMap) {
        Color val = null;
        //check if a color has been blocked and if it has been moves to the next color, if no valid color defaults to the furthest color
        while (!distanceMap.isEmpty()){
            val = distanceMap.entrySet().stream().min(Map.Entry.comparingByValue()).orElseThrow().getKey();
            Color finalVal = val;
            if (node.dontShareColor.stream().anyMatch(n -> n.color == finalVal || n.finalColor == finalVal) || node.blockColor.contains(finalVal)){
                distanceMap.remove(val);
                continue;
            }
            break;
        }
        //if furthest color gets the furthest color based on color distances stored in the table
        if (!node.furthestColorFromSegment.isEmpty()) {
            double maxDist = 0;
            int index = -1;
            double[] distances = table[palette.colors.indexOf(val)];
            for (int i = 0; i < distances.length; i++) {
                double d = distances[i];
                if (d > maxDist) {
                    maxDist = 0;
                    index = i;
                }
            }
            //assigns that furthest color to all nodes that are furthest from and the color it was recolored to all nodes furthest from those nodes (assumed the color it was assign is the furthest)
            Color furthest = palette.get(index);
            for (Node n : node.furthestColorFromSegment) {
                n.finalColor = furthest;
                for (Node n2: n.furthestColorFromSegment){
                    n2.finalColor = val;
                }
            }
        }
        return val;
    }
}
