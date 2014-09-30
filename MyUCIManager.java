import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import acm.graphics.GOval;

public class MyUCIManager extends ImageManager {
	private static final int PEN_SIZE = 2;
	private final List<Character> chars = new ArrayList<Character>();

	public MyUCIManager() {
		System.out.println("new instance");
	}

	public void paint(Graphics gra) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(
					"ujipenchars2.txt"));
			String line;
			int linecount = 0;
			int filename = 0;
			while ((line = br.readLine()) != null) {
				if (!line.startsWith("WORD"))
					continue;
				char character = line.charAt("WORD ".length());
				chars.add(character);

				int numStrokes = Integer.parseInt(br.readLine().substring(
						"  NUMSTROKES ".length()));
				GOval[][] pts = new GOval[numStrokes][];
				
				int maxX = Integer.MIN_VALUE;
				int maxY = Integer.MIN_VALUE;
				int minX = Integer.MAX_VALUE;
				int minY = Integer.MAX_VALUE;
				
				for (int i = 0; i < pts.length; i++) {
					String ptsLine = br.readLine();
					int searchStart = "  POINTS ".length();
					int searchEnd = ptsLine.indexOf(" ", searchStart);
					int numPts = Integer.parseInt(ptsLine.substring(
							searchStart, searchEnd));
					pts[i] = new GOval[numPts];
					String[] split = ptsLine
							.substring(ptsLine.indexOf("#") + 2).split(" ");
					if (split.length != pts[i].length * 2)
						throw new IllegalArgumentException(
								"Reported number of points does not match actual number of points.");

					for (int j = 0; j < pts[i].length; j++) {
						int x = Integer.parseInt(split[j * 2]);
						if (x > maxX)
							maxX = x;
						if (x < minX)
							minX = x;
						int y = Integer.parseInt(split[j * 2 + 1]);
						if (y > maxY)
							maxY = y;
						if (y < minY)
							minY = y;
						pts[i][j] = new GOval(x, y, PEN_SIZE, PEN_SIZE);
					}
				}
					
					minX = minX - 10;
					minY = minY - 10;

					BufferedImage bi = new BufferedImage(220, 220,
							BufferedImage.TYPE_INT_RGB);
					Graphics2D g = bi.createGraphics();
					g.setColor(Color.WHITE);
					g.fillRect(0, 0, 220, 220);
					g.setColor(Color.BLACK);
					
					int xRange = maxX - minX;
					int yRange = maxY - minY;
					double xRate = xRange / 200.0;
					double yRate = yRange / 200.0;
					
					int startX, startY, endX, endY;
					
					for (int i = 0; i < pts.length; i++) {
						for (int j = 1; j < pts[i].length; j++) {
							startX = (int) Math.round(((int)pts[i][j - 1].getX() - minX) / xRate);
							startY = (int) Math.round(((int)pts[i][j - 1].getY() - minY) / yRate);
							endX = (int) Math.round(((int)pts[i][j].getX() - minX) / xRate);
							endY = (int) Math.round(((int)pts[i][j].getY() - minY) / yRate);
							
							g.setStroke(new BasicStroke(15));
							g.drawLine(startX, startY, endX, endY);
						}
					}
					
					gra.drawImage(bi, 0, 0, this);
					System.out.println(Resources.filenames[filename % 194] + Integer.toString(filename + 1) + ".jpeg");
					File file = new File(Resources.filenames[filename % 194] + Integer.toString(filename + 1) + ".jpeg");
					ImageIO.write(bi, "jpg", file);

					linecount++;
					filename++;

				}
			System.out.println("finished");
			System.out.println(linecount);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public int[] getLabels(int limit) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ImageData getImages(int limit) {
		// TODO Auto-generated method stub
		return null;
	}

}