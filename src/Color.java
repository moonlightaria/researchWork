import java.util.Random;
//class for storing colors, processing stores colors as an integer, first byte represents alpha, second red, third, green, last blue
public class Color {
    int r;
    int b;
    int g;
    int a;
    //region constructors
    public Color(int red, int green, int blue, int alpha) {
        r = red;
        g = green;
        b = blue;
        a = alpha;
    }
    public Color(int color){
        a = ((byte) (color >>> 24))+128;
        r = ((byte) (color >>> 16))+128;
        g = ((byte) (color >>> 8))+128;
        b = ((byte) color)+128;
    }
    public Color(String hexColor) {
        a = 255;
        r = Integer.parseInt(hexColor.substring(0,2), 16);
        g = Integer.parseInt(hexColor.substring(2,4), 16);
        b = Integer.parseInt(hexColor.substring(4,6), 16);
    }
    public Color() {
        r = 255;
        g = 255;
        b = 255;
        a = 255;
    }
    //endregion
    @Override
    public String toString(){
        return "Color a: " + a + " r: " + r + " g: " + g + " b: " + b;
    }
    @Override
    public boolean equals(Object o){
        if (this == o) return true;
        if (o == null) return false;
        if (getClass() != o.getClass()) return false;
        Color c = (Color) o;
        return this.r == c.r && this.g == c.g && this.b == c.b && this.a == c.a;
    }
    //color distance in rbga space
    public static double colorDistance(Color c1, Color c2){
        int deltaR = c1.r - c2.r;
        int deltaG = c1.g - c2.g;
        int deltaB = c1.b - c2.b;
        int deltaAlpha = c1.a - c2.a;
        double rgbDistanceSquared = (deltaR * deltaR + deltaG * deltaG + deltaB * deltaB) / 3.0;
        return deltaAlpha * deltaAlpha / 2.0 + rgbDistanceSquared * c1.a * c2.a / (255 * 255);
    }
    public static double colorDistance(int c1 ,int c2){
        int deltaAlpha = (((byte) (c1 >>> 24))+128) - (((byte) (c2 >>> 24))+128);
        int deltaR = (((byte) (c1 >>> 16))+128) - (((byte) (c2 >>> 16))+128);
        int deltaG = (((byte) (c1 >>> 8))+128) - (((byte) (c2 >>> 8))+128);
        int deltaB = (((byte) (c1))+128) - (((byte) (c2))+128);
        double rgbDistanceSquared = (deltaR * deltaR + deltaG * deltaG + deltaB * deltaB) / 3.0;
        return deltaAlpha * deltaAlpha / 2.0 + rgbDistanceSquared * (c1 << 24) * (c2 << 24) / (255 * 255);
    }
    //color distance found in roza's code, works in rgb space
    static double rozaColorDiff(Color c1, Color c2){
        return Math.sqrt((c1.r - c2.r)*(c1.r - c2.r) + (c1.g - c2.g)*(c1.g - c2.g) + (c1.b - c2.b)*(c1.b - c2.b));
    }
    static double rozaColorDiff(int c1, int c2){
        int deltaR = (((byte) (c1 >>> 16))+128) - (((byte) (c2 >>> 16))+128);
        int deltaG = (((byte) (c1 >>> 8))+128) - (((byte) (c2 >>> 8))+128);
        int deltaB = (((byte) (c1))+128) - (((byte) (c2))+128);
        return Math.sqrt((deltaR*deltaR) + (deltaG*deltaG) + (deltaB*deltaB));
    }

    public int getRed() {
        return r;
    }
    public void setRed(int r) {
        this.r = r;
    }
    public int getBlue() {
        return b;
    }
    public void setBlue(int b) {
        this.b = b;
    }
    public int getGreen() {
        return g;
    }
    public void setGreen(int g) {
        this.g = g;
    }
    public int getAlpha(){
        return a;
    }
    public void setAlpha(int a){
        this.a = a;
    }
    //returns a color in processing form
    public int processingForm() {
        int value = a & 0xFF;
        value = (value << 8) + (r & 0xFF);
        value = (value << 8) + (g & 0xFF);
        value = (value << 8) + (b & 0xFF);
        return value;
    }
    //returns a random color
    public static int random(){
        Random rand = new Random();
        Color c = new Color(rand.nextInt(256),rand.nextInt(256),rand.nextInt(256),255);
        return c.processingForm();
    }
    //returns a random color with a random alpha
    public static int randomA(){
        Random rand = new Random();
        Color c = new Color(rand.nextInt(256),rand.nextInt(256),rand.nextInt(256),rand.nextInt(256));
        return c.processingForm();
    }

    public void setRed(String s){
        setRed(validateColorPart(s));
    }
    public void setGreen(String s){
        setGreen(validateColorPart(s));
    }
    public void setBlue(String s){
        setBlue(validateColorPart(s));
    }
    public void setAlpha(String s){
        setAlpha(validateColorPart(s));
    }
    //validates if the string input is a valid integer (between 0 and 255) and corrects it if it isn't
    public int validateColorPart(String s){
        try{
            if (Integer.parseInt(s) <255){
                return Math.max(Integer.parseInt(s), 0);
            }else{
                return 255;
            }
        }catch (NumberFormatException e){
            return 0;
        }
    }

    public Color copy() {
        return new Color(r,g,b,a);
    }
}
