import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontFormatException;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import javax.swing.JComponent;

import math.geom2d.*;

public class Viewer extends JComponent {
	private static final long serialVersionUID = 1L;

	private Font myfont;

	public Viewer() {
		try {
			InputStream is = new FileInputStream("VAGRoundedBold.ttf");
			this.myfont = Font.createFont(Font.TRUETYPE_FONT, is).deriveFont(
					60f);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (FontFormatException e) {
			e.printStackTrace();
		}
	}

	@Override
	public Dimension getMinimumSize() {
		return new Dimension(600, 600);
	}

	@Override
	public Dimension getPreferredSize() {
		return new Dimension(600, 600);
	}

	@Override
	public void paintComponent(Graphics g) {
		Dimension dim = getSize();
		super.paintComponent(g);
		
		String letter = "a";

		// Draw on a buffered image
		BufferedImage img = new BufferedImage(dim.width, dim.height,
				BufferedImage.TYPE_INT_ARGB);
		{
			Graphics2D g2 = (Graphics2D) img.getGraphics();

			boolean smoothDrawing = false;
			if (smoothDrawing) {
				RenderingHints rh = new RenderingHints(
						RenderingHints.KEY_TEXT_ANTIALIASING,
						RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
				g2.setRenderingHints(rh);
			}

			g2.setColor(Color.white);
			g2.fillRect(0, 0, dim.width, dim.height);
			g2.setFont(myfont.deriveFont((float) dim.width));
			g2.setColor(Color.black);
			g2.drawString(letter, 40, dim.width - 50);
		}

		// Process
		{
			int iw = img.getWidth();
			int ih = img.getHeight();
			byte[][] imgData = new byte[ih][iw];

			// note that image is processed row by row top to bottom
			for (int y = 0; y < ih; y++) {
				for (int x = 0; x < iw; x++) {
					Color c = new Color(img.getRGB(x, y));
					imgData[y][x] = (byte) (255 - c.getRed());
				}
			}

			byte[] oneDArray = new byte[imgData.length * imgData.length];
			int s = 0;
			for (int i = 0; i < ih; i++) {
				for (int j = 0; j < iw; j++) {
					oneDArray[s] = imgData[i][j];
					s++;
				}
			}
			
			MarchingSquares ms = new MarchingSquares(oneDArray, ih, iw);
			List<java.awt.geom.Point2D.Double> path = ms.identifyPerimeter().getPath();

			Graphics2D g2 = (Graphics2D) img.getGraphics();
			for(java.awt.geom.Point2D.Double point : path){
				
				img.setRGB((int)point.x, (int)point.y, Color.red.getRGB());
			}
			
			System.out.println("Processed.");
		}

		g.drawImage(img, 0, 0, img.getWidth(), img.getHeight(), null);
	}
}
