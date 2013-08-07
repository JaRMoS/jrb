package models.rbm_advec_tc;

import rb.affinefcn.IAffineFunctions;
import rb.affinefcn.IAffineInitials;
import rb.affinefcn.ITransient;

/**
 * Affine coefficient functions for time-independent advection/diffusion problem
 * 
 * @author Daniel Wirtz
 * @date Aug 29, 2011
 * 
 */
public class AffineFunctions implements IAffineFunctions, ITransient, IAffineInitials {

	/**
	 * @param p
	 * @return Some value
	 */
	public double get_SCM_LB(double[] p) {
		double p_min = p[0];

		for (int i = 1; i < p.length; i++) {
			if (p[i] < p_min) {
				p_min = p[i];
			}
		}

		return 1.;// p_min;
	}

	// ///////////// Implementation of new IAffineFunctions interface
	/**
	 * @see rb.affinefcn.IAffineFunctions#getQl()
	 */
	public int[] getQl() {
		return new int[] { 1 };
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#getQf()
	 */
	public int getQf() {
		return 6;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#getQa()
	 */
	public int getQa() {
		return 2;
	}

	// ///////// Implementation of new ITransient interface

	/**
	 * @see rb.affinefcn.ITransient#getQm()
	 */
	public int getQm() {
		return 1;
	}

	/**
	 * @see rb.affinefcn.ITransient#thetaQm(int, double[], double)
	 */
	@Override
	public double thetaQm(int i, double[] p, double t) {
		return 1.;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#isTimeDependentL()
	 */
	@Override
	public boolean isTimeDependentL() {
		return false;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#isTimeDependentAF()
	 */
	@Override
	public boolean isTimeDependentA() {
		return false;
	}

	/**
	 * @see rb.affinefcn.ITransient#isTimeDependentM()
	 */
	@Override
	public boolean isTimeDependentM() {
		return false;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#getNumOutputs()
	 */
	@Override
	public int getNumOutputs() {
		return 1;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#thetaQl(int, int, double[], double)
	 */
	@Override
	public double thetaQl(int k, int i, double[] p, double t) {
		if ((k < 0) || (k > getQl().length - 1)) {
			throw new RuntimeException("Input parameter is invalid in evaluateL()," + " k = " + k
					+ " but get_n_outputs() = " + getQl());
		}

		if ((i < 0) || (i > getQl()[i] - 1)) {
			throw new RuntimeException("Input parameter is invalid in evaluateL()," + " q_l = " + i
					+ " but get_Q_l(i) = " + getQl()[i]);
		}

		return 1;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#thetaQf(int, double[], double)
	 */
	@Override
	public double thetaQf(int i, double[] p, double t) {
		if ((i < 0) || (i > getQf() - 1)) {
			throw new RuntimeException("Input parameter is invalid in evaluateF()," + " i = " + i
					+ " but get_n_F_functions() = " + getQf());
		}

		switch (i) {
		case 0:
			return 0;
		case 1:
			return 0;
		case 2:
			return p[0];
		case 3:
			return p[1];
		case 4:
			return 0;
		case 5:
			return 0;
		default:
			throw new Error("Should not reach here");
		}
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#thetaQa(int, double[], double)
	 */
	@Override
	public double thetaQa(int i, double[] p, double t) {
		if ((i < 0) || (i > getQa() - 1)) {
			throw new RuntimeException("Input parameter is invalid in evaluateA()," + " i = " + i
					+ " but get_n_A_functions() = " + getQa());
		}

		switch (i) {
		case 0:
			return p[0];
		case 1:
			return p[1];
		default:
			throw new Error("Should not reach here");
		}
	}

	/**
	 * @see rb.affinefcn.IAffiniInitials#getQu0()
	 */
	@Override
	public int getQu0() {
		return 3;
	}

	/**
	 * @see rb.affinefcn.IAffiniInitials#thetaQu0(int i, double[] p)
	 */
	@Override
	public double thetaQu0(int i, double[] p) {
		switch (i) {
		case 0:
			return 0;
		case 1:
			return 1;
		case 2:
			return 0;
		default:
			throw new Error("Should not reach here");
		}
	}
}
