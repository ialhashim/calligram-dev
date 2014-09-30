import javax.swing.JPanel;

public abstract class ImageManager extends JPanel{
	public final int[] getLabels() {
		return getLabels(-1);
	}
	
	public abstract int[] getLabels(int limit);
	
	public final ImageData getImages() {
		return getImages(-1);
	}
	
	public abstract ImageData getImages(int limit);
}