// From https://github.com/QuanCat/marching_squares_algorithms

import java.util.ArrayList;

public class MarchingSquares {
	/*
	 * @param width
	 *            the width of the data matrix
	 * @param height
	 *            the width of the data matrix
	 * @param data
	 *            the data elements
	 */
	private final byte[] data;
	private final int height;
	private final int width;
	
	public MarchingSquares(byte[] data, int height, int width) {
		this.data = data;
		this.height = height;
		this.width = width;
	}
	// accessors

	public byte[] getData() {
		return data;
	}

	public int getHeight() {
		return height;
	}

	public int getWidth() {
		return width;
	}
	
	//method
	public Path identifyPerimeter(int initialX, int initialY) {
		if (initialX < 0) {initialX = 0;}
		if (initialX > width) {initialX = width;}
		if (initialY < 0) {initialY = 0;}
		if (initialY > height) {initialY = height;}
		
		int initialValue = value(initialX, initialY);
		if (initialValue == 0 || initialValue == 15)
			throw new IllegalArgumentException(String.format("Supplied initial coordinates (%d, %d) do not lie on a perimeter.", initialX, initialY));

		ArrayList<Direction> directions = new ArrayList<Direction>();
		int x = initialX;
		int y = initialY;
		Direction previous = null;
		
		do {
			final Direction direction;
			switch (value(x, y)) {
				case  1: direction = Direction.N; break;
				case  2: direction = Direction.E; break;
				case  3: direction = Direction.E; break;
				case  4: direction = Direction.W; break;
				case  5: direction = Direction.N; break;
				case  6: direction = previous == Direction.N ? Direction.W : Direction.E; break;
				case  7: direction = Direction.E; break;
				case  8: direction = Direction.S; break;
				case  9: direction = previous == Direction.E ? Direction.N : Direction.S; break;
				case 10: direction = Direction.S; break;
				case 11: direction = Direction.S; break;
				case 12: direction = Direction.W; break;
				case 13: direction = Direction.N; break;
				case 14: direction = Direction.W; break;
				default: throw new IllegalStateException();
			}
			directions.add(direction);
			x += direction.screenX;
			y += direction.screenY; // accomodate change of basis
			previous = direction;
		} while (x != initialX || y != initialY);

		return new Path(initialX, -initialY, directions);

	}
	
	public Path identifyPerimeter() {
		int size = width * height;
		for (int i = 0; i < size; i++) {
			if (data[i] != 0) {
				return identifyPerimeter(i % width, i / width);
			}
		}
		return null;
	}

	// private utility methods
	
	private int value(int x, int y) {
		int sum = 0;
		if (isSet(x, y)) sum |= 1;
		if (isSet(x + 1, y)) sum |= 2;
		if (isSet(x, y + 1)) sum |= 4;
		if (isSet(x + 1, y + 1)) sum |= 8;
		return sum;
	}

	private boolean isSet(int x, int y) {
		return x <= 0 || x > width || y <= 0 || y > height ?
			false :
			data[(y - 1) * width + (x - 1)] != 0;
	}
	
	
}
