package models.rbm_advec_cheatbase;

import rb.affinefcn.IAffineFunctions;
import rb.affinefcn.IAffineInitials;
import rb.affinefcn.ITransient;

/**
 * Affine coefficient functions for time-dependent advection/diffusion problem
 * 
 * @author Daniel Wirtz
 * @date Aug 29, 2011
 * 
 */
public class AffineFunctions implements IAffineFunctions, ITransient, IAffineInitials {

	/**
	 * TODO: precompute spline coefficients for each i=0..2
	 * 
	 * @param i
	 * @param p
	 * @param t
	 * @return
	 */
	private double dirichletBoundCoeffs(int i, double[] p, double t) {
		// double Q = 3.0;
		// double h = 1. / (Q - 1);
		double h = 0.5;
		double x0 = h * i;

		// Spline coefficients
		double h3 = h * h * h;
		double a1 = 2 / -h3;
		double b1 = (6 * x0 - 3 * h) / h3;
		double c1 = (6 * x0 * x0 - 6 * x0 * h) / -h3;
		double d1 = 1 + (-2 * x0 * x0 * x0 + 3 * x0 * x0 * h) / -h3;
		double a2 = 2 / h3;
		double b2 = (6 * x0 + 3 * h) / -h3;
		double c2 = (6 * x0 * x0 + 6 * x0 * h) / h3;
		double d2 = (2 * x0 * x0 * x0 + 3 * x0 * x0 * h - h3) / -h3;

		double x = p[0];
		double S1 = a1 * x * x * x + b1 * x * x + c1 * x + d1;
		double S2 = a2 * x * x * x + b2 * x * x + c2 * x + d2;
		double res = 0;
		if ((x0 - h <= x) && (x < x0))
			res += S1 * (1 - t);
		if ((x0 <= x) && (x <= x0 + h))
			res += S2 * (1 - t);
		return res;
	}

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

	/**
	 * @see rb.affinefcn.IAffineFunctions#getNumOutputs()
	 */
	@Override
	public int getNumOutputs() {
		return 1;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#getQa()
	 */
	public int getQa() {
		return 2;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#getQf()
	 */
	public int getQf() {
		return 6;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#getQl()
	 */
	public int[] getQl() {
		return new int[] { 1 };
	}

	/**
	 * @see rb.affinefcn.ITransient#getQm()
	 */
	public int getQm() {
		return 1;
	}

	/**
	 * @see rb.affinefcn.IAffiniInitials#getQu0()
	 */
	@Override
	public int getQu0() {
		return 3;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#isTimeDependentAF()
	 */
	@Override
	public boolean isTimeDependentA() {
		return true;
	}

	/**
	 * @see rb.affinefcn.IAffineFunctions#isTimeDependentL()
	 */
	@Override
	public boolean isTimeDependentL() {
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
	 * @see rb.affinefcn.IAffineFunctions#thetaQa(int, double[], double)
	 */
	@Override
	public double thetaQa(int i, double[] p, double t) {
		if ((i < 0) || (i > getQa() - 1)) {
			throw new RuntimeException("Input parameter is invalid in thetaQa()," + " i = " + i + " but getQa() = "
					+ getQa());
		}
		return p[i + 1] * (1 - t);
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

		// flux_mat entspricht Theta_a
		// Udir sind identisch mit den Theta_u0 (zu finden in
		// my_dirichlet_values_coefficients)
		//
		// -------------------------------------------------
		// Q_v = length(flux_mat);
		// Udir = model.dirichlet_values_ptr([],model);
		// Q_Udir = length(Udir);
		// bdir_E_conv = zeros(Q_Udir * Q_v,1);
		// for q1 = 1:Q_Udir
		// for q2 = 1:Q_v
		// bdir_E_conv((q1-1)*Q_v+ q2) = Udir(q1)*flux_mat(q2);
		// end;
		// end;
		// ---------------------------------------------
		//
		// bdir_E_conv liefert dann die endg�ltigen Theta_f

		int didx = (int) Math.floor(i / getQa());
		int Qaidx = i - getQa() * didx;
		return thetaQa(Qaidx, p, t) * dirichletBoundCoeffs(didx, p, t);
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
	 * @see rb.affinefcn.ITransient#thetaQm(int, double[], double)
	 */
	@Override
	public double thetaQm(int i, double[] p, double t) {
		return 1.;
	}

	/**
	 * @see rb.affinefcn.IAffiniInitials#thetaQu0(int i, double[] p)
	 */
	@Override
	public double thetaQu0(int i, double[] p) {
		return dirichletBoundCoeffs(i, p, 0);
	}
}
