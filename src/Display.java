import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PShape;

import java.util.*;
//depreciated version of the interface
public class Display extends PApplet {
    //region declarations
    //region initialize fields
    final int WIDTH = 1800;
    final int HEIGHT = 900;
    final int CONTROL_SPACE = 400;
    final String imageName = "12003";
    final String paletteName = "purple";
    //endregion
    //region draw fields
    final int INDENT = WIDTH - 375;
    final int TEXT_OFFSET = 35;
    final int HORIZONTAL_SPACING = 30;
    final int VERTICAL_SPACING = 30;
    final int CIRCLE_RADIUS = 5;
    PShape buttonShape;
    PShape paletteShape;
    PFont textFont;
    //endregion
    //region state data
    DisplayData data;
    Graph selectedGraph;
    Color selectedColor;
    Node selectedNode;
    String inputText = "";
    boolean selectingEdges;
    //endregion
    //region viewStates
    boolean segmentView;
    boolean graphView;
    //endregion
    //region actions
    int action;
    String[] buttonText;
    final byte NUM_ACTIONS = 21;
    final byte ADD_COLOR_TO_PALETTE = 0;
    final byte REMOVE_COLOR_FROM_PALETTE = 1;
    final byte CHANGE_COLORS_RED = 2;
    final byte CHANGE_COLORS_GREEN = 3;
    final byte CHANGE_COLORS_BLUE = 4;
    final byte CHANGE_COLORS_ALPHA = 5;
    final byte PRINT_NEIGHBORS = 6;
    final byte ADD_EDGES_TO_SELECTED = 7;
    final byte ASSIGN_SELECTED_COLOR_TO_SELECTED = 8;
    final byte ASSIGN_RANDOM_COLOR_TO_SELECTED = 9;
    final byte DONT_COLOR_SELECTED = 10;
    final byte FILTERED_IMAGE = 11;
    final byte TOGGLE_SEGMENT_VIEW = 12;
    final byte TOGGLE_GRAPH_VIEW = 13;
    final byte TOGGLE_WIDEST_PATH = 14;
    final byte TOGGLE_SHORTEST_PATH = 15;
    final byte TOGGLE_KRUSKAL = 16;
    final byte TOGGLE_ROZA_GRAPH = 17;
    final byte TOGGLE_SEGMENT_GRAPH = 18;
    final byte RESET_IMAGE = 19;
    final byte RECOLOR_IMAGE = 20;
    final ArrayList<Byte> toggable = new ArrayList<>();
    boolean[] toggled;
    //endregion
    //endregion
    //region setup
    public void settings() {
        size(WIDTH,HEIGHT);
    }
    public void setup() {
        Collections.addAll(toggable, ADD_EDGES_TO_SELECTED, FILTERED_IMAGE,TOGGLE_SEGMENT_VIEW, TOGGLE_GRAPH_VIEW,
                TOGGLE_WIDEST_PATH, TOGGLE_SHORTEST_PATH, TOGGLE_KRUSKAL, TOGGLE_ROZA_GRAPH, TOGGLE_SEGMENT_GRAPH);
        toggled = new boolean[toggable.size()];
        data = new DisplayData(WIDTH-CONTROL_SPACE,HEIGHT, imageName, paletteName);
        //region paletteShape
        paletteShape = createShape();
        paletteShape.beginShape();
        paletteShape.vertex(0,0);
        paletteShape.vertex(20,0);
        paletteShape.vertex(20,20);
        paletteShape.vertex(0,20);
        paletteShape.endShape();
        //endregion
        //region buttonShape
        buttonShape = createShape();
        buttonShape.beginShape();
        buttonShape.fill(50);
        buttonShape.vertex(0,0);
        buttonShape.vertex(20,0);
        buttonShape.vertex(20,20);
        buttonShape.vertex(0,20);
        buttonShape.endShape();
        //endregion
        textFont = createFont("Arial",14);
        selectedGraph = data.originalSegmentGraph;
        selectedNode = selectedGraph.nodes.get(0);
        selectedColor = data.palette.get(0);
    }
    //endregion
    //region draw
    public void draw(){
        background(255);
        if (data.graphLoaded) {
            drawGraph();
        }
        drawControls();
    }
    private void drawControls() {
        fill(255);
        rect(WIDTH-CONTROL_SPACE,0,CONTROL_SPACE,HEIGHT);
        fill(0);
        textFont(textFont);
        buttonText = setText();
        text("Input: " + inputText, INDENT,40);
        int textPosY = 75;
        for (int i = 0; i < buttonText.length; i++) {
            buttonShape.setFill(color(0));
            int check = toggable.indexOf((byte) i);
            if (check != -1){
                if (toggled[check]){
                    buttonShape.setFill(color(255,0,0));
                }
            }
            shape(buttonShape, INDENT, 60 + i * VERTICAL_SPACING);
            text(buttonText[i], INDENT + TEXT_OFFSET, textPosY + 30*i);
        }
        if(data.paletteLoaded) {
            for (int i = 0; i < data.palette.size(); i++) {
                paletteShape.setFill(data.palette.get(i).processingForm());
                shape(paletteShape, INDENT + HORIZONTAL_SPACING * i, HEIGHT-50);
            }
        }
    }
    private String[] setText() {
        String[] buttonText = new String[NUM_ACTIONS];
        buttonText[ADD_COLOR_TO_PALETTE] = "add color to palette";
        buttonText[REMOVE_COLOR_FROM_PALETTE] = "remove current color from palette";
        buttonText[CHANGE_COLORS_RED] = "red: ";
        buttonText[CHANGE_COLORS_GREEN] = "green: ";
        buttonText[CHANGE_COLORS_BLUE] = "blue: ";
        buttonText[CHANGE_COLORS_ALPHA] = "alpha: ";
        if (selectedColor != null) {
            buttonText[CHANGE_COLORS_RED] = buttonText[CHANGE_COLORS_RED] + str(selectedColor.getRed());
            buttonText[CHANGE_COLORS_GREEN] = buttonText[CHANGE_COLORS_GREEN] + str(selectedColor.getGreen());
            buttonText[CHANGE_COLORS_BLUE] = buttonText[CHANGE_COLORS_BLUE] + str(selectedColor.getBlue());
            buttonText[CHANGE_COLORS_ALPHA] = buttonText[CHANGE_COLORS_ALPHA] +str(selectedColor.getAlpha());
        }
        buttonText[PRINT_NEIGHBORS] = "neighbor: ";
        buttonText[ADD_EDGES_TO_SELECTED] = "add edges: ";
        buttonText[ASSIGN_SELECTED_COLOR_TO_SELECTED] = "assign current color to selected";
        buttonText[ASSIGN_RANDOM_COLOR_TO_SELECTED] = "assign random color to selected";
        buttonText[DONT_COLOR_SELECTED] = "don't recolor selected: ";
        if (selectedNode != null){
            buttonText[PRINT_NEIGHBORS] = buttonText[PRINT_NEIGHBORS] + selectedNode;
            buttonText[DONT_COLOR_SELECTED] = buttonText[DONT_COLOR_SELECTED] + str(selectedNode.dontColor);
            if (selectedGraph !=null) {
                buttonText[ADD_EDGES_TO_SELECTED] = buttonText[ADD_EDGES_TO_SELECTED] + selectedGraph.getNeighborsText(selectedNode);
            }
            if (selectedNode.color != null) {
                buttonText[ASSIGN_SELECTED_COLOR_TO_SELECTED] = buttonText[ASSIGN_SELECTED_COLOR_TO_SELECTED] + selectedNode.color;
                buttonText[ASSIGN_SELECTED_COLOR_TO_SELECTED] = buttonText[ASSIGN_RANDOM_COLOR_TO_SELECTED] + selectedNode.color;
            }
        }
        buttonText[FILTERED_IMAGE] = "filtered: " + str(data.filtered);
        buttonText[TOGGLE_SEGMENT_VIEW] = "segment view";
        buttonText[TOGGLE_GRAPH_VIEW] = "graph view";
        buttonText[TOGGLE_WIDEST_PATH] = "widest path";
        buttonText[TOGGLE_SHORTEST_PATH] = "shortest path";
        buttonText[TOGGLE_KRUSKAL] = "kruskal";
        buttonText[TOGGLE_ROZA_GRAPH] = "roza graph";
        buttonText[TOGGLE_SEGMENT_GRAPH] = "segment graph";
        buttonText[RESET_IMAGE] = "reset image";
        buttonText[RECOLOR_IMAGE] = "recolorImage";
        return buttonText;
    }
    private void drawGraph() {
        if (segmentView){
            drawSegmentView();
        }else {
            drawNormalView();
        }
        if (graphView) {
            drawGraphView();
        }
    }
    private void drawGraphView() {
        for (Node n: data.segmentGraph().nodes) {
            circle(n.centroid.x,n.centroid.y,CIRCLE_RADIUS);
        }
        if (segmentView){
            stroke(0);
        }
        for (HashSet<Edge> hashSet: selectedGraph.edges) {
            for (Edge e : hashSet) {
                line(e.nodeA.centroid.x, e.nodeA.centroid.y, e.nodeB.centroid.x, e.nodeB.centroid.y);
            }
        }
    }
    private void drawSegmentView() {
        if (data.scaledSegmentImage != null) {
            image(data.scaledSegmentImage, 0, 0);
        }
        drawSegmentBorder(data.closestNode(mouseX,mouseY));
    }
    private void drawNormalView() {
        if (data.scaledOriginalImage != null) {
            image(data.scaledOriginalImage, 0, 0);
        }
        drawSegmentBorder(data.closestNode(mouseX,mouseY));
    }
    private void drawSegmentBorder(Node n) {
        stroke(color(255,0,0,150));
        fill(color(255,0,0,150));
        for (Point p: n.scaledVertices) {
            rect(p.x, p.y, 2, 2);
        }
    }
    //endregion
    //region input
    public void mousePressed(){
        if (!inControlSpace(mouseX)) {
            if (selectingEdges) {
                selectedGraph.addEdge(selectedNode, data.closestNode(mouseX, mouseY));
            } else {
                selectedNode = data.closestNode(mouseX, mouseY);
            }
            return;
        }
        for (int i =0; i<data.palette.size(); i++){
            if (Algorithms.collision(mouseX - (INDENT + i*HORIZONTAL_SPACING),mouseY-(HEIGHT-50), paletteShape)){
                selectedColor = data.palette.get(i);
                return;
            }
        }
        for (int i = 0; i < buttonText.length; i++) {
            if (Algorithms.collision(mouseX - INDENT, mouseY - 60 - VERTICAL_SPACING * i, buttonShape)) {
                action = i;
                action();
            }
        }
    }
    private boolean inControlSpace(int mouseX) {
        return mouseX>WIDTH-CONTROL_SPACE;
    }
    public void keyPressed() {
        if (key == BACKSPACE){
            inputText = inputText.substring(0,inputText.length()-1);
        }else if (Character.isDigit(key)){
            inputText = inputText + key;
        }else{
            shortcuts(key);
        }
    }
    private void shortcuts(char key) {
        switch (Character.toLowerCase(key)){
            case ('+') -> action = ADD_COLOR_TO_PALETTE;
            case ('-') -> action = REMOVE_COLOR_FROM_PALETTE;
            case ('r') -> action = CHANGE_COLORS_RED;
            case ('g') -> action = CHANGE_COLORS_GREEN;
            case ('b') -> action = CHANGE_COLORS_BLUE;
            case ('a') -> action = CHANGE_COLORS_ALPHA;
            case ('p') -> action = PRINT_NEIGHBORS;
            case ('e') -> action = ADD_EDGES_TO_SELECTED;
            case ('c') -> action = ASSIGN_SELECTED_COLOR_TO_SELECTED;
            case (',') -> action = ASSIGN_RANDOM_COLOR_TO_SELECTED;
            case ('s') -> action = TOGGLE_SEGMENT_VIEW;
            case ('d') -> action = TOGGLE_GRAPH_VIEW;
            case ('/') -> action = RESET_IMAGE;
        }
        action();
    }
    //endregion
    //region action
    private void action() {
        switch (action){
            case (ADD_COLOR_TO_PALETTE) -> data.palette.add(new Color());
            case (REMOVE_COLOR_FROM_PALETTE) -> data.palette.remove(selectedColor);
            case (CHANGE_COLORS_RED) -> changeColorsRed();
            case (CHANGE_COLORS_GREEN) -> changeColorsGreen();
            case (CHANGE_COLORS_BLUE) -> changeColorsBlue();
            case (CHANGE_COLORS_ALPHA) -> changeColorsAlpha();
            case (PRINT_NEIGHBORS) -> System.out.println(selectedGraph.getNeighborsText(selectedNode));
            case (ADD_EDGES_TO_SELECTED) -> addEdgesToSelected();
            case (ASSIGN_SELECTED_COLOR_TO_SELECTED) -> assignColor(selectedColor,selectedNode);
            case (ASSIGN_RANDOM_COLOR_TO_SELECTED) -> assignColor(data.palette.random(), selectedNode);
            case (DONT_COLOR_SELECTED) -> selectedNode.dontColor = !selectedNode.dontColor;
            case (FILTERED_IMAGE) -> filtered();
            case (TOGGLE_SEGMENT_VIEW) -> toggleSegmentView();
            case (TOGGLE_GRAPH_VIEW) -> toggleGraphView();
            case (TOGGLE_WIDEST_PATH) -> toggleWidestPath();
            case (TOGGLE_SHORTEST_PATH) -> toggleShortestPath();
            case (TOGGLE_KRUSKAL) -> toggleKruskal();
            case (TOGGLE_ROZA_GRAPH) -> toggleRozaGraph();
            case (TOGGLE_SEGMENT_GRAPH) -> toggleSegmentGraph();
            case (RESET_IMAGE) -> reset();
            case (RECOLOR_IMAGE) -> RozaCode.makeGraphonRegions(data, selectedGraph, selectedNode);
            default -> throw new Error("no valid action");
        }
        System.out.println("action: " + action);
    }

    private void filtered() {
        data.filtered = !data.filtered;
        toggled[toggable.indexOf(FILTERED_IMAGE)] = data.filtered;
        selectedNode = data.segmentGraph().nodes.get(selectedNode.nodeId);
        switch(selectedGraph.type){
            case segment -> selectedGraph = data.segmentGraph();
            case shortest -> selectedGraph = data.shortestPath(selectedNode);
            case widest -> selectedGraph = data.widestPath(selectedNode);
            case kruskal -> selectedGraph = data.kruskal();
            case roza -> selectedGraph = data.rozaGraph(selectedNode);
        }
    }
    private void changeColorsRed() {
        selectedColor.setRed(inputText);
        inputText = "";
    }
    private void changeColorsGreen() {
        selectedColor.setGreen(inputText);
        inputText = "";
    }
    private void changeColorsBlue(){
        selectedColor.setBlue(inputText);
        inputText = "";
    }
    private void changeColorsAlpha(){
        selectedColor.setAlpha(inputText);
        inputText = "";
    }
    private void addEdgesToSelected() {
        selectingEdges = !selectingEdges;
        toggled[toggable.indexOf(ADD_EDGES_TO_SELECTED)] = !toggled[toggable.indexOf(ADD_EDGES_TO_SELECTED)];
    }
    private void toggleSegmentView() {
        segmentView = !segmentView;
        toggled[toggable.indexOf(TOGGLE_SEGMENT_VIEW)] = !toggled[toggable.indexOf(TOGGLE_SEGMENT_VIEW)];
    }
    private void toggleGraphView(){
        graphView = !graphView;
        toggled[toggable.indexOf(TOGGLE_GRAPH_VIEW)] = !toggled[toggable.indexOf(TOGGLE_GRAPH_VIEW)];
    }
    private void toggleWidestPath() {
        selectedGraph = data.widestPath(selectedNode);
        setAllGraphsFalse();
        toggled[toggable.indexOf(TOGGLE_WIDEST_PATH)] = true;
    }
    private void toggleShortestPath() {
        selectedGraph = data.shortestPath(selectedNode);
        setAllGraphsFalse();
        toggled[toggable.indexOf(TOGGLE_SHORTEST_PATH)] = true;
    }
    private void toggleKruskal() {
        selectedGraph = data.kruskal();
        setAllGraphsFalse();
        toggled[toggable.indexOf(TOGGLE_KRUSKAL)] = true;
    }
    private void toggleRozaGraph(){
        selectedGraph = data.rozaGraph(selectedNode);
        setAllGraphsFalse();
        toggled[toggable.indexOf(TOGGLE_ROZA_GRAPH)] = true;
    }
    private void toggleSegmentGraph() {
        selectedGraph = data.segmentGraph();
        setAllGraphsFalse();
        toggled[toggable.indexOf(TOGGLE_SEGMENT_GRAPH)] = true;
    }
    private void setAllGraphsFalse(){
        toggled[toggable.indexOf(TOGGLE_WIDEST_PATH)] = false;
        toggled[toggable.indexOf(TOGGLE_SHORTEST_PATH)] = false;
        toggled[toggable.indexOf(TOGGLE_KRUSKAL)] = false;
        toggled[toggable.indexOf(TOGGLE_ROZA_GRAPH)] = false;
        toggled[toggable.indexOf(TOGGLE_SEGMENT_GRAPH)] = false;
    }
    private void rbRecoloring(){
        LinkedList<Node> current;
        LinkedList<Node> next = new LinkedList<>();
        HashSet<Node> queued = new HashSet<>();
        boolean red = true;
        next.add(selectedNode);
        queued.add(selectedNode);
        while (!next.isEmpty()) {
            current = next;
            next = new LinkedList<>();
            for (Node n : current) {
                if (red){
                    assignColor(new Color(255,0,0,255),n);
                }else{
                    assignColor(new Color(0,0,0,255),n);
                }
                for (Node node: selectedGraph.getNeighbors(n)){
                    if (!queued.contains(node)){
                        queued.add(node);
                        next.add(node);
                    }
                }
            }
            red = !red;
        }
        data.originalImage.updatePixels();
        data.updateImage();
    }
    private void assignColor(Color c, Node n){
        n.finalColor = c;
    }
    private void randomizeColors(){
        for (Node n: selectedGraph.nodes){
            assignColor(data.palette.random(),n);
        }
    }
    private void reset(){
        data.reset();
        selectedGraph =data.segmentGraph();
        selectedNode = selectedGraph.nodes.get(0);
    }
    //endregion
    public static void main(String[] args){
        PApplet.main("Display");
    }
}