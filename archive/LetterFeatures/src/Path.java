
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.awt.Point;
import java.awt.geom.Point2D;

public class Path {
	private static final double ADJ_LEN = Math.sqrt(2.0) / 2.0 - 1;
	private final Direction[] directions;
	private final List<Direction> directionList;
	private double length;
	private final int originX;
	private final int originY;
	private final int terminalX;
	private final int terminalY;
	
	//constructor
	public Path(Path that, int deltaX, int deltaY) {
		super();
		this.directions = that.directions;
		this.directionList = that.directionList;
		this.length = that.length;
		this.originX = that.originX;
		this.originY = that.originY;
		this.terminalX = that.terminalX;
		this.terminalY = that.terminalY;
	}
	
	public Path(int startX, int startY, Direction[] directions) {
		this.originX = startX;
		this.originY = startY;
		this.directions = directions.clone();
		this.directionList = Collections.unmodifiableList(Arrays.asList(directions));
		int endX = startX;
		int endY = startY;
		int diagonals = 0;
		for (Direction direction : directions) {
			endX += direction.screenX;
			endY += direction.screenY;
			if (direction.screenX != 0 && direction.screenY != 0) {
				diagonals++;
			}
		}
		
		this.terminalX = endX;
		this.terminalY = endY;
		this.length = directions.length + diagonals * ADJ_LEN;
	}
	
	public Path(int startX, int startY, List<Direction> directions) {
		this(startX, startY, directions.toArray(new Direction[directions.size()]));
	}
	
	// accessors
	public double getLength() {
		return length;
	}

	public Direction[] getDirections() {
		return directions;
	}
	
	public ArrayList<Point> getPath(){
		ArrayList<Point> path = new ArrayList<Point>();
		int x = originX;
		int y = originY;
		
		for(Direction d : directions){
			x += d.planeX;
			y += d.planeY;
			path.add(new Point(Math.abs(x),Math.abs(y)));
		}
		
		return path;
	}

	public int getOriginX() {
		return originX;
	}

	public int getOriginY() {
		return originY;
	}

	public int getTerminalX() {
		return terminalX;
	}

	public int getTerminalY() {
		return terminalY;
	}
	
	public boolean isClosed() {
		return originX == terminalX && originY == terminalY;
	}
	
	public Path translate(int deltaX, int deltaY) {
		return new Path(this, deltaX, deltaY);
	}
	
	/*
	 * Two paths are equal if they have the same origin and the same directions.
	 */
	
	@Override
	public boolean equals(Object obj) {
		if (obj == this) {
			return true;
		}
		if (!(obj instanceof Path)) {
			return false;
		}
		Path that = (Path) obj;
		
		if (this.originX != that.originX) {
			return false;
		}
		if (this.originY != that.originY) {
			return false;
		}
		if (this.terminalX != that.terminalX) {
			return false;
		}
		if (this.terminalY != that.terminalY) {
			return false;
		}
		if (!Arrays.equals(this.directions, that.directions)) {
			return false;
		}
		return true;
	}
	
	@Override
	public int hashCode() {
		return originX ^ 7 * originY ^ directions.hashCode();
	}

	@Override
	public String toString() {
		return "X: " + originX + ", Y: " + originY + " " + 
				Arrays.toString(directions);
	}
	

}
