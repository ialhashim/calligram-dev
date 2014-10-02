import math.geom2d.*; // documentation http://geom-java.sourceforge.net/api/index.html

import java.awt.*;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
 
public class LetterFeatures implements Runnable {
 
    @Override
    public void run() {
        // Create the window
        JFrame f = new JFrame("Hello, World!");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setLayout(new FlowLayout());
        
        try{
        	InputStream is = new FileInputStream("VAGRoundedBold.ttf");
        	Font myfont = Font.createFont(Font.TRUETYPE_FONT, is).deriveFont(42f);
        	f.setFont(myfont);
        	
            JTextField userTextField=new JTextField(20);
            userTextField.setFont(myfont);
            userTextField.setText("Testing!");
            f.add(userTextField);
        }catch(IOException e)
        {
        	e.printStackTrace();
        } catch (FontFormatException e) {
			e.printStackTrace();
		}       
        
        f.pack();
        f.setVisible(true);
    }
 
    public static void main(String[] args) {
    	LetterFeatures se = new LetterFeatures();
        // Schedules the application to be run at the correct time in the event queue.
        SwingUtilities.invokeLater(se);
    }
 
}