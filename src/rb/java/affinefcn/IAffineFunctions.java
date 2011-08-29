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
	 * Flag that indicates whether the output functional affine theta
	 * coefficient functions are time dependent or not. The time parameter has
	 * to be passed every time anyways, but computations can be made mode
	 * efficient if the matrices stay identical for all time-steps
	 * 
	 * @return Whether \theta^l are time-dependent or not.
	 */
	public boolean isTimeDependentL();

	/**
	 * Flag that tells JRB if the bilinear form/input functional affine theta
	 * coefficient functions are time-dependent or not. The time parameter has
	 * to be passed every time anyways, but computations can be made mode
	 * efficient if the matrices stay identical for all time-steps
	 * 
	 * @return Whether \theta^{a,f} are time-dependent or not.
	 */
	public boolean isTimeDependentAF();

	/**
	 * Returns the size of the affine linear combination for each output
	 * functional
	 * 
	 * @return Ql for output functional "output_number"
	 */
	public int[] getQl();
	
	/**
	 * Returns the number of output functionals.
	 * @return The number of output functionals
	 */
	public int getNumOutputs();

	/**
	 * Evaluates \theta^{q_l}_i for the output q_l and the current parameter p.
	 * 
	 * @param k
	 *            The output index
	 * @param i
	 *            The number of the \theta function
	 * @param p
	 *            The current parameter
	 * @param t
	 *            The current time
	 * @return The value of \theta^l_i at (p,t) for output k
	 */
	public double thetaQl(int k, int i, double[] p, double t);

	/**
	 * Gets the number of input functionals in the affine expansion.
	 * 
	 * @return The number of input functionals in the affine expansion.
	 */
	public int getQf();

	/**
	 * Evaluates \theta^f_i for the current parameter p
	 * 
	 * @param i
	 *            The coefficient function i
	 * @param p
	 *            The current parameter
	 * @param t
	 *            The current time
	 * @return The value of \theta^f_i at (p,t)
	 */
	public double thetaQf(int i, double[] p, double t);

	/**
	 * Gets the number of bilinear forms in the affine expansion.
	 * 
	 * @return The number of bilinear forms in the affine expansion.
	 */
	public int getQa();

	/**
	 * Evaluates \theta^a_i for the current parameter p
	 * 
	 * @param i
	 *            The bilinar form i
	 * @param p
	 *            The current parameter
	 * @param t
	 *            The current time t
	 * @return The value of \theta^a_i at (p,t)
	 */
	public double thetaQa(int i, double[] p, double t);
}
