import java.awt.*;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;
 
public class LetterFeatures implements Runnable {
 
    @Override
    public void run() {
        // Create the window
        JFrame f = new JFrame("Hello, World!");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setLayout(new FlowLayout());
        
        f.add(new Viewer());
        
        f.pack();
        f.setVisible(true);
    }
 
    public static void main(String[] args) {
    	LetterFeatures se = new LetterFeatures();
        // Schedules the application to be run at the correct time in the event queue.
        SwingUtilities.invokeLater(se);
    }
 
}