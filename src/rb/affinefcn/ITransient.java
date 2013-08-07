package rb.affinefcn;

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
	 * Flag that determines if the theta^m coefficient functions are time dependent.
	 * 
	 * @return True if any theta^m_i is time dependent, false otherwise
	 */
	public boolean isTimeDependentM();

	/**
	 * 
	 * @param i
	 * The mass coefficient function i
	 * @param p
	 * The current paramter p
	 * @param t
	 * The time t
	 * @return
	 */
	public double thetaQm(int i, double[] p, double t);

	/**
	 * 
	 * @return
	 */
	public int getQm();
}
