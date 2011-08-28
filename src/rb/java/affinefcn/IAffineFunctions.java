/**
 * Created on Aug 28, 2011 in Project JRB
 * Location: rb.java.affinefcn.IAffineFunctions.java
 */
package rb.java.affinefcn;

/**
 * Base interface for any affine functions used as an externally loaded class.
 * 
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 *
 */
public interface IAffineFunctions {
	
	/**
	 * Returns the size of the affine linear combination for each output functional  
	 *
	 * @return Ql for output functional "output_number" 
	 */
	public int[] getQl();
	
	/**
	 * Evaluates \theta^{q_l}_i for the output q_l and the current parameter p.
	 * 
	 * @param i
	 * @param q_l
	 * @param p
	 * @return
	 */
	public double evaluateL(int i, int q_l, double[] p);
	
	/**
	 * Gets the number of input functionals in the affine expansion.
	 * @return The number of input functionals in the affine expansion.
	 */
	public int getQf();
	
	/**
	 * Evaluates \theta^f_i for the current parameter p
	 * @param i
	 * @param p
	 * @return
	 */
	public double evaluateF(int i, double[] p);
	
	/**
	 * Gets the number of bilinear forms in the affine expansion.
	 * @return The number of bilinear forms in the affine expansion.
	 */
	public int getQa();
	
	/**
	 * Evaluates \theta^a_i for the current parameter p
	 * @param i
	 * @param p
	 * @return
	 */
	public double evaluateA(int i, double[] p);
}
