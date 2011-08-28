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
	 * @param p
	 * @return
	 */
	public double evaluateM(int i, double[] p);
	
	/**
	 * 
	 * @return
	 */
	public int getQm();
}
