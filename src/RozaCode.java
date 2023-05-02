import processing.core.PImage;
//optimized and adapted version of roza's original code, currently unused
public class RozaCode {
    static double minPaleDif = Double.MAX_VALUE;
    static double MaxPaleDif = Double.MIN_VALUE;
    public static void makeGraphonRegions(InterfaceData data){
        long start = System.currentTimeMillis();
        PImage draw = data.alteredImage.copy();
        StringBuilder fileName = new StringBuilder();
        fileName.append(data.imageName);
        fileName.append("_");
        fileName.append(data.palette.name);
        fileName.append("_");
        switch (data.path.type) {
            case segment -> fileName.append("segment");
            case shortest -> fileName.append("shortest");
            case widest -> fileName.append("widest");
            case roza -> fileName.append("roza");
            case kruskal -> fileName.append("kruskal");
        }
        data.start.colored = true;
        traverseGraph(data.path, data.start);
        for (Node n : data.path.nodes) {
            n.colored = false;
        }
        RecoloringByBottleNeck(data.path, data.start, draw, data.palette);
        for (Node n: data.path.nodes){
            n.color = null;
            n.colored = false;
            n.parent = null;
        }
        draw.save("C:/Users/edch6/Desktop/code stuff/Research_Work/src/output/" + fileName + ".png");
        System.out.println("Recoloring: total time elapsed: " + (System.currentTimeMillis()-start));
    }

    public static void makeGraphonRegions(DisplayData data, Graph selectedGraph, Node selectedNode){
        long start = System.currentTimeMillis();
        PImage draw;
        StringBuilder fileName = new StringBuilder();
        fileName.append(data.imageName);
        fileName.append("_");
        fileName.append(data.paletteName);
        fileName.append("_");
        if (data.filtered){
            fileName.append("filtered_");
            draw = data.filteredImage.copy();
        }else{
            fileName.append("original_");
            draw = data.originalImage.copy();
        }
        switch (selectedGraph.type) {
            case segment -> fileName.append("segment");
            case shortest -> fileName.append("shortest");
            case widest -> fileName.append("widest");
            case roza -> fileName.append("roza");
            case kruskal -> fileName.append("kruskal");
        }
        selectedNode.colored = true;
        traverseGraph(selectedGraph, selectedNode);
        for (Node n : selectedGraph.nodes) {
            n.colored = false;
        }
        RecoloringByBottleNeck(selectedGraph, selectedNode, draw, data.palette);
        for (Node n: selectedGraph.nodes){
            n.color = null;
            n.colored = false;
            n.parent = null;
        }
        draw.save("C:/Users/edch6/Desktop/code stuff/Research_Work/src/output/" + fileName + ".png");
        System.out.println("Recoloring: total time elapsed: " + (System.currentTimeMillis()-start));
    }

    private static void traverseGraph(Graph graph, Node current) {
        for (Edge e: graph.edges.get(current.nodeId)){
            if (!e.nodeB.colored) {
                e.nodeB.parent = e.nodeA;
                e.nodeB.colored = true;
                traverseGraph(graph, e.nodeB);
            }
        }
    }

    static void RecoloringByBottleNeck(Graph spanningTree, Node source, PImage image, Palette palette) {
        double[][] table = new double[palette.size()][palette.size()];
        for (int i = 0; i < palette.size(); i++){
            for (int j = 0; j < palette.size(); j++){
                table[i][j] = Color.rozaColorDiff(palette.get(i),palette.get(j));
                if (minPaleDif > table[i][j]) {
                    minPaleDif = table[i][j];
                } else if (MaxPaleDif < table[i][j]) {
                    MaxPaleDif = table[i][j];
                }
            }
        }
        int[] lookUp_table = HistogramMatching(spanningTree.getEdgeWeights(), table);
        AssignColors(spanningTree, source, palette, lookUp_table, source);

        for (Node node: spanningTree.nodes) {
            if (!node.dontColor) {
                for (Point p : node.points) {
                    image.set(p.x, p.y, node.color.processingForm());
                }
            }
        }
    }
    static int[] HistogramMatching(double[] edgeWeights, double[][] table){
        int bsize = 500;
        int[] rfreq = new int[bsize];
        int[] tfreq = new int[bsize];

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
    static void AssignColors(Graph graph, Node root, Palette palette, int[] lookUp_table, Node node){
        if (node.finalColor != null){
            node.color = node.finalColor;
        }else {
            node.color = ReturnClosestColor(graph, root, palette, lookUp_table, node);
        }
        node.colored = true;
        for (Edge e: graph.edges.get(node.nodeId)) {
            if (!e.nodeB.colored){
                AssignColors(graph, root, palette, lookUp_table, e.nodeB);
            }
        }
    }
    private static Color ReturnClosestColor(Graph graph, Node root, Palette palette, int[] lookUp_table, Node node) {
        Color returnVal = null;
        double minDif = 1000;
        boolean signr,signp;
        double g1, g2, p1, p2;
        Color pidex = null;
        double colorDif, diff, difp;
        Color col1, col2;
        col2 = node.avgColor;
        if (node.parent == null || node.nodeId == root.nodeId) {
            for (Edge e : graph.edges.get(node.nodeId)) {
                if (e.nodeB.colored) {
                     returnVal = e.nodeB.color;
                }
            }
            if (returnVal ==null) {
                for (Color c : palette.colors) {
                    diff = Color.rozaColorDiff(col2,c);
                    if (diff < minDif) {
                        minDif = diff;
                        returnVal = c;
                    }
                }
            }
        }
        if (node.parent != null){
            col1 = node.parent.avgColor;
            pidex = node.parent.color;
            p1 = (0.2126 * pidex.r) + (0.7152 * pidex.g) + (0.0722 * pidex.b);
            g1 = (0.2126 *col1.r) + (0.7152 *col1.g) + (0.0722 * col1.b);	// return luminance
            g2 = (0.2126 *col2.r) + (0.7152 *col2.g) + (0.0722 * col2.b);	// return luminance
            colorDif = lookUp_table[(int) Math.round(Color.rozaColorDiff(col1,col2))];
            signr = g1 >= g2;
            minDif = 1000;
            for (Color c: palette.colors) {
                difp = Color.rozaColorDiff(pidex,c);
                p2 = (0.2126 *c.r) + (0.7152 *c.g) + (0.0722 * c.b);
                signp = p1>=p2;
                diff = Math.abs(colorDif - difp);
                if (signr == signp) {
                    if (diff <minDif) {
                        minDif = diff;
                        returnVal = c;
                    }
                }
            }
            if (returnVal == null) {
                minDif = 1000;
                for (Color c: palette.colors) {
                    difp = Color.rozaColorDiff(pidex,c);
                    diff = Math.abs(colorDif - difp);
                    if (diff < minDif) {
                        minDif = diff;
                        returnVal = c;
                    }
                }
            }
        }
        return returnVal;
    }
}