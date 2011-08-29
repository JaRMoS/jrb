/**
 * Created on Aug 28, 2011 in Project JRB
 * Location: rb.java.affinefcn.ITransient.java
 */
package rb.java.affinefcn;

/**
 * Interface for AffineFunctions in unsteady rb systems.
 * 
 * Provides access to the coefficient functions for the mass matrices.
 * 
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 * 
 */
public interface ITransient {

	/**
	 * 
	 * @param i
	 *            The mass coefficient function i
	 * @param p
	 *            The current paramter p
	 * @param t
	 *            The time t
	 * @return
	 */
	public double thetaQm(int i, double[] p, double t);

	/**
	 * 
	 * @return
	 */
	public int getQm();
}
