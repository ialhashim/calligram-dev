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

import javax.swing.JComponent;

import math.geom2d.*;

public class Viewer extends JComponent {
	private static final long serialVersionUID = 1L;

	private Font myfont;

	public Viewer() {
		try {
			InputStream is = new FileInputStream("VAGRoundedBold.ttf");
			this.myfont = Font.createFont(Font.TRUETYPE_FONT, is).deriveFont(60f);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (FontFormatException e) {
			e.printStackTrace();
		}
	}

	@Override
	public Dimension getMinimumSize() {
		return new Dimension(600,600);
	}

	@Override
	public Dimension getPreferredSize() {
		return new Dimension(600,600);
	}

	@Override
	public void paintComponent(Graphics g) {
		Dimension dim = getSize();
		super.paintComponent(g);

		// Draw on a buffered image
		BufferedImage img = new BufferedImage(dim.width, dim.height, BufferedImage.TYPE_INT_ARGB);
		{
			Graphics2D g2 = (Graphics2D) img.getGraphics();
			
			boolean smoothDrawing = false;
			if( smoothDrawing )
			{
				RenderingHints rh = new RenderingHints(
						RenderingHints.KEY_TEXT_ANTIALIASING,
						RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
				g2.setRenderingHints(rh);
			}
	
			g2.setColor(Color.white);
			g2.fillRect(0, 0, dim.width, dim.height);
			g2.setFont(myfont.deriveFont((float)dim.width));
			g2.setColor(Color.black);
			g2.drawString("b", 40, dim.width - 50);
		}
		
		// Process		
		{
			int iw = img.getWidth();
	        int ih = img.getHeight();
	        int[][] imgData = new int[ih][iw];
	        
	        // note that image is processed row by row top to bottom
	        for(int y = 0; y < ih; y++) {
	            for(int x = 0; x < iw; x++) {
	            	Color c = new Color(img.getRGB(x,y));
	            	imgData[y][x] = c.getRed();
	            }
	        }
	        
	        System.out.println("Processed.");
		}
        
		g.drawImage(img, 0, 0, img.getWidth(), img.getHeight(), null);
	}
}
