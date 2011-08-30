/**
 * 
 */
package rb.java.affinefcn;

/**
 * Interface for affinely decomposed initial value conditions.
 * 
 * @author dwirtz
 *
 */
public interface IAffineInitials {
	
	/**
	 * Gets the number of input functionals in the affine expansion.
	 * 
	 * @return The number of input functionals in the affine expansion.
	 */
	public int getQu0();
	
	/**
	 * Evaluates \theta^{u0}_i for the current parameter p
	 * 
	 * @param i
	 *            The coefficient function i
	 * @param p
	 *            The current parameter
	 * @return The value of \theta^{u0}_i for parameter p
	 */
	public double thetaQu0(int i, double[] p);

}
