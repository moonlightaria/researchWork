import static processing.core.PApplet.loadStrings;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
//class to store a palette of colors
public class Palette {
    //name of the palette
    final String name;
    //list to store all colors in the palette
    List<Color> colors;
    //color lover info currently unused
    int id;
    int votes;
    double hearts;
    public Palette(String[] data, String name){
        this.name = name;
        id = Integer.parseInt(data[0].substring(12));
        votes = Integer.parseInt(data[1].substring(14));
        hearts = Double.parseDouble(data[2].substring(15));
        colors = new ArrayList<>();
        for (int i = 4; i < data.length; i++) {
            colors.add(new Color(data[i]));
        }
    }

    public Palette(File file) {
        this(loadStrings(file),file.getName());
    }

    public int size() {
        return colors.size();
    }
    public Color get(int i){
        return colors.get(i);
    }
    //methods for adding and removing from the palette, not used in current version
    public void remove(int i) {
        colors.remove(i);
    }
    public void remove(Color color) {
        colors.remove(color);
    }
    public void add(Color color) {
        colors.add(color);
    }
    //returns a random color in the palette
    public Color random(){
        Random rand = new Random();
        return colors.get(rand.nextInt(colors.size()));
    }
    private void debug(){
        System.out.println(id);
        System.out.println(votes);
        System.out.println(hearts);
        System.out.println(colors.size());
    }
}
