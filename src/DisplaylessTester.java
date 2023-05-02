import java.util.Arrays;
import java.util.List;
//depreciated tester class
public class DisplaylessTester {
    static final String[] imageNames = {"12003", "65019", "175083", "bike1", "bird", "dog1", "fall3", "indian1024", "metal2"};
    static final String[] paletteNames = {"basic", "purple", "brown"};
    static void recolorTest(){
        long startTime = System.currentTimeMillis();
        DisplayData data = new DisplayData(1400,900, "bike1", "brown");
        RozaCode.makeGraphonRegions(data, data.widestPath(data.originalSegmentGraph.findMidReg()), data.originalSegmentGraph.findMidReg());
        System.out.println("Recolor Test: Total time elapsed: " +  ((System.currentTimeMillis()-startTime)));
    }
    static void recolorAll(){
        long startTime = System.currentTimeMillis();
        for (String imageName: imageNames) {
            DisplayData data = new DisplayData(1400,900,imageName);
            for (String paletteName: paletteNames){
                data.filtered = false;
                Node selectedNode = data.originalSegmentGraph.findMidReg();
                data.setPalette(paletteName);

                RozaCode.makeGraphonRegions(data, data.shortestPath(selectedNode),  selectedNode);
                RozaCode.makeGraphonRegions(data, data.widestPath(selectedNode),  selectedNode);
                RozaCode.makeGraphonRegions(data, data.kruskal(), selectedNode);
                RozaCode.makeGraphonRegions(data, data.rozaGraph(selectedNode),  selectedNode);

                data.filtered = true;
                selectedNode = data.filteredSegmentGraph.findMidReg();
                RozaCode.makeGraphonRegions(data, data.shortestPath(selectedNode),  selectedNode);
                RozaCode.makeGraphonRegions(data, data.kruskal(), selectedNode);
                RozaCode.makeGraphonRegions(data, data.widestPath(selectedNode),  selectedNode);
                RozaCode.makeGraphonRegions(data, data.rozaGraph(selectedNode),  selectedNode);
            }
        }
        System.out.println("time for recoloring: " + (System.currentTimeMillis()-startTime));
    }
    static void debug(){
        List<Point> points = Arrays.asList(new Point(5, 0), new Point(5,50), new Point(5,100), new Point(10,50));
        System.out.println(Algorithms.findInteriorPoint(points));
    }
    public static void main(String[] args){
        debug();
    }
}
