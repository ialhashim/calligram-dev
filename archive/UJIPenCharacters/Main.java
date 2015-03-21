import javax.swing.JFrame;


public class Main {

	public static void main(String[] args) {
		
		MyUCIManager myuci = new MyUCIManager();
		
		JFrame frame = new JFrame();
		frame.getContentPane().add(myuci);
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    frame.setSize(400,400);
	    frame.setVisible(true);
	}
}
