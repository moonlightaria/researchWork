//class to store a point on an image or on the interface, also used to store a pairs of integers for other purposes
public class Point {
    int x;
    int y;

    public Point(int x, int y) {
        this.x = x;
        this.y = y;
    }
    public Point copy() {
        return new Point(x,y);
    }
    public String toString(){
        return "Point x: " + x + " y: " + y;
    }
    public int hashCode(){
        return x*10000+y;
    }
    public boolean equals(Object o){
        if (this == o) return true;
        if (o == null) return false;
        if (getClass() != o.getClass()) return false;
        Point p = (Point) o;
        if (this.x == p.x) {
            return this.y == p.y;
        }
        return false;
    }
}
