import processing.core.*;

import javax.imageio.ImageIO;
import java.io.File;
import java.io.IOException;
import java.util.*;

//the user interface for the user to select actions to change the recoloring
public class UserInterface extends PApplet {
    //region initialize fields
    //width and height of UI
    final int WIDTH = 1800;
    final int HEIGHT = 900;
    //size of the space for the control
    final int CONTROL_SPACE = 400;
    //spot where button are drawn
    final int INDENT = WIDTH - 375;
    //offset for where the text is drawn for each button
    final int TEXT_OFFSET = 35;
    //y offset from top for where button text should be drawn
    final int TEXT_POSY = 45;
    //space between palette colors on display
    final int HORIZONTAL_SPACING = 30;
    //offset from top of screen to where buttons are drawn
    final int BUTTON_POSY = 30;
    //space between action button on display
    final int VERTICAL_SPACING = 30;
    //size of the circles used in connection view
    final int CIRCLE_RADIUS = 5;

    //file name of the image without file extension
    String imageName;
    //field to pass the image to the data class (work around to processing input methods returning void)
    PImage selectedImage;
    //field to pass the segment image to the data class
    PImage segmentImage;
    //field to pass the palette file to the data class
    File selectedPalette;
    //the shape used to draw the action buttons
    PShape buttonShape;
    //the shape use to draw the palette
    PShape paletteShape;
    //the font used for all text in the UI
    PFont textFont;
    //Color of the first selection of nodes
    Color selection1Color;
    //Color of the second selection of nodes
    Color selection2Color;
    //Color of the mouse drag box
    Color targetColor;
    //endregion
    //region state data
    //stores all of the data (graph, palette, etc.)
    InterfaceData data;
    //the currently selected color
    Color selectedColor;
    //the first selection of nodes
    HashSet<Node> selection;
    //the second selection of nodes
    HashSet<Node> selection2;
    //the path draw when using the abstract selection
    LinkedList<Point> selectionPath;
    //the start and end point when using the box selection
    Point boxStart;
    Point boxEnd;
    //whether the user is using the second selection
    boolean secondSelection;
    //whether the user should be able to use the controls (disabled as lazy way to prevent odd behavior)
    boolean controlActive;
    //if the user is using the box method to select nodes
    boolean boxSelection;
    //if the user is using the abstract drawing to select nodes
    boolean abstractSelection;
    //if the user is adding the nodes to the existing selection
    boolean addToExisting;
    //if the user is removing nodes from the existing selection
    boolean removeFromExisting;
    //booleans for what view the user is currently using
    boolean segmentView;
    boolean connectionView;
    boolean pathView;
    //enum to store all the possible actions a user could take
    enum Action{
        CONNECT_SEGMENTS(0,"connect all selected segments to target segments: ", true),
        CONNECT_ALL(1, "connect all selected segments to each other", false),
        COLOR_SELECTED(2, "assign selected color to all selected segments: ", false),
        BLOCK_COLOR(3, "prevents selected color from being assigned to selected segments", false),
        DONT_COLOR_SELECTED(4, "don't recolor selected region, selected regions will still be assigned a color for purpose of recolor other segments", false),
        SHARE_COLOR(5, "all selected region will be recolored to the same color", false),
        DONT_SHARE_COLOR(6, "prevents all selected segments from sharing a color with target segments", false),
        BIAS_TOWARDS_COLOR(7, "bias the recoloring of all selected segments towards a selected color", false),
        BIAS_FROM_COLOR(8, "bias the recoloring of all selected segments away from a selected color", false),
        BIAS_TOWARDS_SEGMENT(9, "bias the recoloring of all selected segments towards the color of a target segment, does nothing if a selected segment is assigned a color before the target", false),
        BIAS_FROM_SEGMENT(10, "bias the recoloring of all selected segments away from the color of a target segment, does nothing if a selected segment is assigned a color before the target", false),
        FURTHEST_COLOR_FROM_SEGMENT(11, "recolors the two segment to the colors furthest away from each other in the palette", false),
        TOGGLE_SEGMENT_VIEW(12, "toggle the view to show the segments", true),
        TOGGLE_CONNECTION_VIEW(13, "toggle the view to show the connections between regions", true),
        TOGGLE_PATH(14,"toggle the view to show the generated path", true),
        GENERATE_PATH(15, "generates the path for recoloring based on the current connections starting from selected segments", false),
        PARTIAL_PATH(16, "generates a path solely including the nodes in the selection using the nodes on the border as start points", false),
        RECOLOR_IMAGE(17, "recolor the image based on an existing path, if no path has been generate a path will be generated based on default settings", false),
        RESET_IMAGE(18, "reset all user selections", false),
        SAVE_IMAGE(19,"save the current image", false);

        //corresponds to the index of the array, used to look up actions and for index in corresponding arrays
        final int num;
        //the text explaining the action to the user, appear next to the button
        final String text;
        //if the action is toggable (has two states i.e. on or off for a view or if the user is selection a second selection)
        final boolean toggable;

        Action(int num, String text, boolean toggable){
            this.num = num;
            this.text = text;
            this.toggable = toggable;
        }
        //returns the corresponding action based on its position in the array
        static Action getAction(int n){
            return Arrays.stream(Action.values()).filter(a -> a.num == n).findFirst().orElse(null);
        }
    }
    //the last action the user took
    Action lastAction;
    //array to store whether an action has been toggled on
    boolean[] toggled;
    //endregion
    //region setup
    //set processing settings
    public void settings() {
        size(WIDTH,HEIGHT);
    }
    //setup before display runs
    public void setup() {
        //first argument is text prompt second argument is the name of the method being called as a string
        /*selectInput("select a palette for recoloring", "selectPalette");
        selectInput("select a corresponding set of segments", "selectSegments");
        selectInput("select a image to recolor", "selectImage");*/
        defaultInit();
        data = new InterfaceData(WIDTH-CONTROL_SPACE,HEIGHT, selectedImage, imageName, segmentImage, selectedPalette);
        //region paletteShape
        //defines the shape used for displaying the palette
        paletteShape = createShape();
        paletteShape.beginShape();
        paletteShape.vertex(0,0);
        paletteShape.vertex(20,0);
        paletteShape.vertex(20,20);
        paletteShape.vertex(0,20);
        paletteShape.endShape();
        //endregion
        //region buttonShape
        //defines the shape for displaying the buttons
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
        selection = new HashSet<>();
        selection2 = new HashSet<>();
        selection1Color = new Color (255,0,0,255);
        selection2Color = new Color (0,255,0,255);
        targetColor = new Color (0,0,255,255);
        toggled = new boolean[Action.values().length];
        controlActive = true;
        //sets display to update once per 15 seconds, changing can help deal with lag
        frameRate(15);
    }
    //region IO
    //default option to remove having to select image files when testing
    void defaultInit(){
        String paletteName = "brown";
        imageName = "12003";
        try {
            selectedImage = new PImage(ImageIO.read(new File("./src/images/originals/" + imageName + ".jpg")));
            segmentImage = new PImage(ImageIO.read(new File("./src/images/segmentation/" + imageName + ".png")));
            selectedPalette = new File("./src/palettes/" + paletteName + ".txt");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    //methods to select files
    public void selectImage(File file){
        if (file == null){
            System.out.println("Error: no input selected");
            System.exit(0);
        }
        imageName = file.getName();
        selectedImage = loadImage(file.getAbsolutePath());
        if (selectedImage == null){
            System.out.println("ERROR: selected file not in supported format");
            System.exit(0);
        }
    }
    public void selectSegments(File file){
        if (file == null){
            System.out.println("Error: no input selected");
            System.exit(0);
        }
        segmentImage = loadImage(file.getAbsolutePath());
        if (selectedImage == null){
            System.out.println("ERROR: selected file not in supported format");
            System.exit(0);
        }
    }
    public void selectPalette(File file){
        if (file == null){
            System.out.println("Error: no input selected");
            System.exit(0);
        }
        selectedPalette = file;
    }
    //endregion
    //endregion
    //region draw
    //overrides processing draw, called once per frame
    public void draw(){
        if (data != null) {
            //clear the screen
            background(255);
            if (data.graphLoaded) {
                drawGraph();
            }
            drawControls();
        }
    }
    //draws the controls
    private void drawControls() {
        //clears the space for the display
        fill(255);
        rect(WIDTH-CONTROL_SPACE,0,CONTROL_SPACE,HEIGHT);
        //sets color for text
        fill(0);
        textFont(textFont);
        Action[] actions = Action.values();
        //draws all buttons and text
        for (int i = 0; i < actions.length; i++) {
            //sets button color,
            if (toggled[i]){
                buttonShape.setFill(color(255,0,0));
            }else{
                buttonShape.setFill(color(0));
            }
            //draws button shape and text
            shape(buttonShape, INDENT, BUTTON_POSY + i * VERTICAL_SPACING);
            text(actions[i].text, INDENT + TEXT_OFFSET, TEXT_POSY + VERTICAL_SPACING * i);
        }
        //draws palette
        if(data.palette != null) {
            for (int i = 0; i < data.palette.size(); i++) {
                //sets color to that of palette and draws it to screen
                paletteShape.setFill(data.palette.get(i).processingForm());
                shape(paletteShape, INDENT + HORIZONTAL_SPACING * i, HEIGHT-(VERTICAL_SPACING+paletteShape.height));
            }
        }
    }
    //draws the images and the graph overlay
    private void drawGraph() {
        //draw the image (segment image or current image)
        if (segmentView){
            drawSegmentView();
        }else {
            drawNormalView();
        }
        //checks whether to draw the path or segment graph based if connection view is active
        if (connectionView) {
            if (pathView && data.pathSelected) {
                drawConnectionView(data.path);
            }else{
                drawConnectionView(data.alteredSegmentGraph);
            }
        }
        //draws the outline for the selections
        for (Node n: selection){
            drawSegmentBorder(n, selection1Color);
        }
        for (Node n: selection2){
            drawSegmentBorder(n, selection2Color);
        }
        //draws the outline of the current outline based on the current mode
        stroke(targetColor.processingForm());
        fill(targetColor.processingForm());
        //draws 4 lines forming a box around the selected area
        if (boxSelection){
            line(boxStart.x, boxStart.y, boxStart.x, boxEnd.y);
            line(boxStart.x, boxStart.y, boxEnd.x, boxStart.y);

            line(boxEnd.x, boxEnd.y, boxEnd.x, boxStart.y);
            line(boxEnd.x, boxEnd.y, boxStart.x, boxEnd.y);
        }
        //draws the shape around what has been selected
        if (abstractSelection){
            for (int current = 0; current< selectionPath.size(); current++) {
                int next = current + 1;
                if (next == selectionPath.size()) next = 0;

                line(selectionPath.get(current).x, selectionPath.get(current).y, selectionPath.get(next).x, selectionPath.get(next).y);
            }
        }
    }
    //draws a representation of the graph
    private void drawConnectionView(Graph graph) {
        //draws a circle at the center of each segment (can be off of the segment in case of holes in the segment)
        for (Node n: graph.nodes) {
            circle(n.centroid.x,n.centroid.y,CIRCLE_RADIUS);
        }
        //changes color if lines to black if segment view is enabled
        if (segmentView){
            stroke(0);
        }
        //draws a line to show every edge on the graph (can cause lag)
        for (HashSet<Edge> hashSet: graph.edges) {
            for (Edge e : hashSet) {
                line(e.nodeA.centroid.x, e.nodeA.centroid.y, e.nodeB.centroid.x, e.nodeB.centroid.y);
            }
        }
    }
    //draws the segment image to the screen
    private void drawSegmentView() {
        if (data.scaledSegmentImage != null) {
            image(data.scaledSegmentImage, 0, 0);
        }
    }
    //draws the current image to the screen
    private void drawNormalView() {
        if (data.scaledAlteredImage != null) {
            image(data.scaledAlteredImage, 0, 0);
        }
    }
    //draws a rough outline around a segment of the given color
    private void drawSegmentBorder(Node n, Color c) {
        stroke(c.processingForm());
        fill(c.processingForm());
        for (Point p: n.scaledVertices) {
            rect(p.x, p.y, 1, 1);
        }
    }
    //endregion
    //region input
    //processing method that runs every time the mouse is pressed
    public void mousePressed(){
        //check if the mouse isn't in the control space (area designated for the controls)
        if (!inControlSpace(mouseX)) {
            //if right mouse is pressed start the box selection, otherwise start the abstract shape selection
            if (mouseButton == RIGHT) {
                boxSelection = true;
                abstractSelection = false;
                boxStart = new Point(mouseX, mouseY);
                boxEnd = boxStart;
            } else {
                boxSelection = false;
                abstractSelection = true;
                selectionPath = new LinkedList<>();
            }
            //check if current selection needs to be cleared
            if (!(addToExisting || removeFromExisting)) {
                if (secondSelection) {
                    selection2 = new HashSet<>();
                }else {
                    selection = new HashSet<>();
                }
            }
            return;
        }
        if (controlActive) {
            //check collision for each color to see if it is selected, if it is changes the current color to that color
            for (int i =0; i<data.palette.size(); i++){
                if (Algorithms.collision(mouseX - (INDENT + i*HORIZONTAL_SPACING),mouseY-(HEIGHT-50), paletteShape)){
                    selectedColor = data.palette.get(i);
                    System.out.println(selectedColor);
                    return;
                }
            }
            //checks collision for each action button, if is selected runs that action
            for (int i = 0; i < Action.values().length; i++) {
                if (Algorithms.collision(mouseX - INDENT, mouseY - (30 + i * VERTICAL_SPACING), buttonShape)) {
                    action(Action.getAction(i));
                    return;
                }
            }
        }
    }
    //updates the user selection when the mouse is moved
    public void mouseDragged(){
        if (boxSelection){
            boxEnd = new Point (mouseX, mouseY);
        }else if (abstractSelection){
            selectionPath.add(new Point (mouseX, mouseY));
        }
    }
    //updates the selection with the nodes contained in the user selection
    public void mouseReleased(){
        HashSet<Node> foundCenters = new HashSet<>();
        if (boxSelection){
            foundCenters = data.getNodes(boxStart, boxEnd);
        }else if (abstractSelection){
            foundCenters = data.getNodes(selectionPath);
        }
        HashSet<Node> sel;
        if (secondSelection){
            sel = selection2;
        }else{
            sel = selection;
        }
        if (removeFromExisting) {
            sel.removeAll(foundCenters);
        }else{
            sel.addAll(foundCenters);
        }
    }
    //checks if mouse is inside the control space
    private boolean inControlSpace(int mouseX) {
        return mouseX>WIDTH-CONTROL_SPACE;
    }
    //processing method that runs whenever a key is pressed
    public void keyPressed() {
        switch (keyCode){
            case CONTROL -> addToExisting = true;
            case SHIFT -> removeFromExisting = true;
            case RETURN, ENTER -> lastAction();
        }
    }
    //processing method that runs whenever a key is released
    public void keyReleased(){
        switch (keyCode){
            case CONTROL -> addToExisting = false;
            case SHIFT -> removeFromExisting = false;
        }
    }
    //endregion
    //region action
    //run the action corresponding to the action passed
    private void action(Action action) {
        lastAction = action;
        switch (action) {
           case CONNECT_SEGMENTS -> connectSegments();
           case CONNECT_ALL -> connectAll();
           case COLOR_SELECTED -> colorSelected();
           case BLOCK_COLOR -> blockColor();
           case DONT_COLOR_SELECTED -> dontColorSelected();
           case SHARE_COLOR -> shareColor();
           case DONT_SHARE_COLOR -> dontShareColor();
           case BIAS_TOWARDS_COLOR -> biasTowardsColor();
           case BIAS_FROM_COLOR -> biasFromColor();
           case BIAS_TOWARDS_SEGMENT -> biasTowardsSegment();
           case BIAS_FROM_SEGMENT -> biasFromSegment();
           case FURTHEST_COLOR_FROM_SEGMENT -> furthestColorFromSegment();
           case TOGGLE_SEGMENT_VIEW -> toggleSegmentView();
           case TOGGLE_CONNECTION_VIEW -> toggleConnectionView();
           case TOGGLE_PATH -> togglePath();
           case GENERATE_PATH -> generatePath();
           case PARTIAL_PATH -> partialPath();
           case RECOLOR_IMAGE -> recolorImage();
           case RESET_IMAGE -> data.reset();
           case SAVE_IMAGE -> data.alteredImage.save(data.getOutputLocation());
        }
        System.out.println("action: " + action);
    }
    //runs the last action
    private void lastAction() {
        action(lastAction);
    }
    //first run switches to second selection, second run adds edges to the segment graph between nodes in selection 1 and selection 2,
    // run again using the last action key, disables control as precaution against unwanted behavior
    private void connectSegments() {
        if (selection.isEmpty()){
            return;
        }
        if (secondSelection){
            secondSelection = false;
            controlActive = true;
            for (Node n1: selection){
                for (Node n2: selection2){
                    data.alteredSegmentGraph.addEdge(n1,n2);
                }
            }
            selection2.clear();
        }else {
            secondSelection = true;
            controlActive = false;
        }
    }
    //adds edges between all segments in a selection, can cause lag
    private void connectAll() {
        for (Node n1 : selection) {
            for (Node n2 : selection) {
                data.alteredSegmentGraph.addEdge(n1, n2);
            }
        }
    }
    //assign the selection color to all segments in the selection
    private void colorSelected() {
        for (Node n: selection){
            n.finalColor = selectedColor;
        }
    }
    //prevents all segments in the selection from being recolored to the selected color
    private void blockColor() {
        for (Node n: selection){
            n.blockColor.add(selectedColor);
        }
    }
    //unused action to assign a random color to every segment in the selection
    private void randomlyColorSelected(){
        for (Node n: selection){
            n.finalColor = data.palette.random();
        }
    }
    //stops all segments in the selection from being recolored
    private void dontColorSelected(){
        for (Node n: selection){
            n.dontColor = !n.dontColor;
        }
    }
    //sets all segments in the selection to share a color upon recoloring
    private void shareColor(){
        for (Node n1: selection){
            n1.shareColor.addAll(selection);
        }
    }
    //prevents all segments in the first selection from sharing a color with a segment in the second selection,
    //causes issues if all colors end up present in one selection before recoloring any node in the other selection,
    //run again by using the last action key
    private void dontShareColor(){
        if (selection.isEmpty()){
            return;
        }
        if (secondSelection){
            secondSelection = false;
            controlActive = true;
            for (Node n1: selection){
                n1.dontShareColor.addAll(selection2);
            }
            for (Node n2: selection2){
                n2.dontShareColor.addAll(selection);
            }
            selection2.clear();
        }else {
            secondSelection = true;
            controlActive = false;
        }
    }
    //biases all segments in the selections recoloring towards the selected color
    private void biasTowardsColor(){
        for (Node n: selection){
            n.biasTowardsColor.add(selectedColor);
        }
    }
    //biases all segments in the selections recoloring away the selected color
    private void biasFromColor(){
        for (Node n: selection){
            n.biasFromColor.add(selectedColor);
        }
    }
    //bias all segments in both selections towards the colors present in the other selection
    private void biasTowardsSegment(){
        if (selection.isEmpty()){
            return;
        }
        if (secondSelection){
            secondSelection = false;
            controlActive = true;
            for (Node n: selection){
                n.biasTowardsSegment.addAll(selection2);
            }
            for (Node n: selection2){
                n.biasTowardsSegment.addAll(selection);
            }
            selection2.clear();
        }else {
            secondSelection = true;
            controlActive = false;
        }
    }
    //bias all segments in both selections away from the colors present in the other selection
    private void biasFromSegment(){
        if (selection.isEmpty()){
            return;
        }
        if (secondSelection){
            secondSelection = false;
            controlActive = true;
            for (Node n: selection){
                n.biasFromSegment.addAll(selection2);
            }
            for (Node n: selection2){
                n.biasFromSegment.addAll(selection);
            }
            selection2.clear();
        }else {
            secondSelection = true;
            controlActive = false;
        }
    }
    //set the segments to recolor to the furthest color from the other segment
    private void furthestColorFromSegment(){
        if (selection.isEmpty()){
            return;
        }
        if (secondSelection){
            secondSelection = false;
            controlActive = true;
            for (Node n: selection){
                n.furthestColorFromSegment.addAll(selection2);
            }
            for (Node n: selection2){
                n.furthestColorFromSegment.addAll(selection);
            }
            selection2.clear();
        }else {
            secondSelection = true;
            controlActive = false;
        }
    }
    //toggles segment view
    private void toggleSegmentView() {
        segmentView = !segmentView;
        toggled[Action.TOGGLE_SEGMENT_VIEW.num] = segmentView;
    }
    //toggles connection view
    private void toggleConnectionView(){
        connectionView = !connectionView;
        toggled[Action.TOGGLE_CONNECTION_VIEW.num] = connectionView;
    }
    //toggles path view
    private void togglePath() {
        pathView = !pathView;
        toggled[Action.TOGGLE_PATH.num] = pathView;
    }
    //generates a path starting from the center most segment of the selection
    private void generatePath(){
         data.generatePath(selection);
         pathView = true;
         toggled[Action.TOGGLE_PATH.num] = true;
    }
    //generates a path starting from all segments that have selection points inside and outside the user selection
    private void partialPath(){
        if (boxSelection){
            data.generatePartialPath(boxStart,boxEnd);
        }else if (abstractSelection){
            data.generatePartialPath(selectionPath);
        }
        pathView = true;
        toggled[Action.TOGGLE_PATH.num] = true;
    }
    //recolors the image based on the current path, if no path exists generates a default path
    private void recolorImage(){
        data.recolor();
    }

    //endregion
    public static void main(String[] args){
        PApplet.main("UserInterface");
    }
}
