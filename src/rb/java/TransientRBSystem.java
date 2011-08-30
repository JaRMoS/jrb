package rb.java;

//    rbAPPmit: An Android front-end for the Certified Reduced Basis Method
//    Copyright (C) 2010 David J. Knezevic and Phuong Huynh
//
//    This file is part of rbAPPmit
//
//    rbAPPmit is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rbAPPmit is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rbAPPmit.  If not, see <http://www.gnu.org/licenses/>. 

import java.io.BufferedReader;
import java.io.IOException;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import rb.java.affinefcn.ITransient;
import rmcommon.Log;
import rmcommon.geometry.DiscretizationType;
import rmcommon.io.AModelManager;
import rmcommon.io.MathObjectReader;

// This class provides the Online reduced
// basis functionality for linear parabolic
// problems.
// This class is modeled on the TransientRBSystem
// class in rbOOmit

/**
 * Changes made by
 * 
 * @author Daniel Wirtz
 * @date Aug 26, 2011
 * 
 */
public class TransientRBSystem extends RBSystem {

	// Logging tag
	private static final String DEBUG_TAG = "TransientRBSystem";

	/**
	 * boolean flag that determines whether or not we impose a temporal filter
	 * on each output
	 */
	private boolean apply_temporal_filter_flag;

	private double[][][][] Aq_Mq_representor_norms;

	private double[][] cached_Aq_Aq_matrix;

	private double[][] cached_Aq_Mq_matrix;

	private double[] cached_Fq_Aq_vector;
	private double[] cached_Fq_Mq_vector;
	private double cached_Fq_term;

	private double[][] cached_Mq_Mq_matrix;
	/**
	 * Time step size.
	 */
	private double dt;
	/**
	 * The error bound for the field variable at each time level.
	 */
	protected double[] error_bound_all_k;
	/**
	 * Parameter that defines the temporal discretization: euler_theta = 0 --->
	 * Forward Euler euler_theta = 0.5 ---> Crank-Nicolson euler_theta = 1 --->
	 * Backward Euler
	 */
	private double euler_theta;
	/**
	 * The standard deviation of the filter (in terms of time steps).
	 */
	private double filter_width;
	/**
	 * Vectors storing the residual representor inner products to be used in
	 * computing the residuals online.
	 */
	private double[][][] Fq_Mq_representor_norms;

	/**
	 * The number of terms in the affine expansion of the mass matrix
	 */
	private int fQm;

	/**
	 * Total number of time-steps.
	 */
	protected int fTotalTimesteps;

	private double[][][] Mq_Mq_representor_norms;

	/**
	 * The number of time steps we actually plot in the output plotter. This is
	 * sometimes less than K so that we don't plot the effect of the filter near
	 * the final time.
	 */
	public int n_plotting_steps;

	/**
	 * RBSystem has a member RB_solution, we also need old_RB_solution here.
	 */
	protected RealVector old_RB_solution;

	// @SuppressWarnings({ "unused" })
	// /**
	// * A secondary SCM object since we might need a lower bound for the mass
	// * matrix and the stiffness matrix.
	// */
	// private RBSCMSystem mSecondRbScmSystem;

	/**
	 * Dense RB mass matrix.
	 */
	protected RealMatrix RB_L2_matrix;
	/**
	 * Dense matrices for the RB computations.
	 */
	protected RealMatrix[] RB_M_q_matrix;
	/**
	 * The current time-level, 0 <= _k <= _K.
	 */
	protected int timestep;
	/**
	 * The solution coefficients at each time level from the most recent
	 * RB_solve.
	 */
	protected RealVector[] timestepRBSolutions;

	/**
	 * Apply the temporal filter to the outputs
	 */
	public void apply_temporal_filter() {
		double[][] RB_convolved_outputs_all_k = new double[getNumOutputs()][getTotalTimesteps() + 1];
		double[][] RB_convolved_error_bounds_all_k = new double[getNumOutputs()][getTotalTimesteps() + 1];

		boolean tdL = affineFunctionsInstance.isTimeDependentL();
		double output_dual_norm = 0;
		for (int n = 0; n < getNumOutputs(); n++) {
			if (!tdL)
				output_dual_norm = compute_output_dual_norm(n, 0);// Zero
																	// is
																	// current
																	// time

			for (int time_level = 0; time_level <= getTotalTimesteps(); time_level++) {
				if (tdL)
					output_dual_norm = compute_output_dual_norm(n, time_level
							* getdt());

				double conv_weight_integral_sq = 0.;
				RB_convolved_outputs_all_k[n][time_level] = 0.;

				for (int k_prime = 0; k_prime <= getTotalTimesteps(); k_prime++) {
					double time_diff = getdt() * (time_level - k_prime);
					RB_convolved_outputs_all_k[n][time_level] += getdt()
							* conv_weight(time_diff)
							* RB_outputs_all_k[n][k_prime];

					conv_weight_integral_sq += getdt()
							* Math.pow(conv_weight(time_diff), 2.);
				}

				RB_convolved_error_bounds_all_k[n][time_level] = error_bound_all_k[getTotalTimesteps()]
						* output_dual_norm * Math.sqrt(conv_weight_integral_sq);
			}
		}

		RB_outputs_all_k = RB_convolved_outputs_all_k;
		RB_output_error_bounds_all_k = RB_convolved_error_bounds_all_k;
	}

	/**
	 * Helper function that caches the time-independent residual quantities.
	 */
	protected void cache_online_residual_terms(int N) {

		cached_Fq_term = 0.;
		int q = 0;
		for (int q_f1 = 0; q_f1 < getQf(); q_f1++) {

			double cached_theta_q_f1 = thetaQf(q_f1);
			for (int q_f2 = q_f1; q_f2 < getQf(); q_f2++) {
				double delta = (q_f1 == q_f2) ? 1. : 2.;
				cached_Fq_term += delta * cached_theta_q_f1 * thetaQf(q_f2)
						* Fq_representor_norms[q];

				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			double cached_theta_q_f = thetaQf(q_f);
			for (int q_a = 0; q_a < getQa(); q_a++) {
				double cached_theta_q_a = thetaQa(q_a);
				for (int i = 0; i < N; i++) {
					// Clear the entries on the first pass
					if ((q_f == 0) && (q_a == 0))
						cached_Fq_Aq_vector[i] = 0.;

					cached_Fq_Aq_vector[i] += 2. * cached_theta_q_f
							* cached_theta_q_a
							* Fq_Aq_representor_norms[q_f][q_a][i];
				}
			}
		}

		q = 0;
		for (int q_a1 = 0; q_a1 < getQa(); q_a1++) {
			double cached_theta_q_a1 = thetaQa(q_a1);
			for (int q_a2 = q_a1; q_a2 < getQa(); q_a2++) {
				double cached_theta_q_a2 = thetaQa(q_a2);
				double delta = (q_a1 == q_a2) ? 1. : 2.;

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						// Clear the entries on the first pass
						if (q == 0)
							cached_Aq_Aq_matrix[i][j] = 0.;

						cached_Aq_Aq_matrix[i][j] += delta * cached_theta_q_a1
								* cached_theta_q_a2
								* Aq_Aq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			double cached_theta_q_f = thetaQf(q_f);

			for (int q_m = 0; q_m < getQm(); q_m++) {
				double cached_theta_q_m = thetaQm(q_m);

				for (int i = 0; i < N; i++) {
					// Clear the entries on the first pass
					if ((q_f == 0) && (q_m == 0))
						cached_Fq_Mq_vector[i] = 0.;

					cached_Fq_Mq_vector[i] += 2. * cached_theta_q_f
							* cached_theta_q_m
							* Fq_Mq_representor_norms[q_f][q_m][i];
				}
			}
		}

		for (int q_a = 0; q_a < getQa(); q_a++) {
			double cached_theta_q_a = thetaQa(q_a);

			for (int q_m = 0; q_m < getQm(); q_m++) {
				double cached_theta_q_m = thetaQm(q_m);

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						// Clear the entries on the first pass
						if ((q_a == 0) && (q_m == 0))
							cached_Aq_Mq_matrix[i][j] = 0.;

						cached_Aq_Mq_matrix[i][j] += 2. * cached_theta_q_a
								* cached_theta_q_m
								* Aq_Mq_representor_norms[q_a][q_m][i][j];
					}
				}
			}
		}

		q = 0;
		for (int q_m1 = 0; q_m1 < getQm(); q_m1++) {
			double cached_theta_q_m1 = thetaQm(q_m1);
			for (int q_m2 = q_m1; q_m2 < getQm(); q_m2++) {
				double cached_theta_q_m2 = thetaQm(q_m2);
				double delta = (q_m1 == q_m2) ? 1. : 2.;

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						if (q == 0)
							cached_Mq_Mq_matrix[i][j] = 0.;

						cached_Mq_Mq_matrix[i][j] += delta * cached_theta_q_m1
								* cached_theta_q_m2
								* Mq_Mq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

	}

	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector.
	 * 
	 * This assumes that the time-independent quantities were cached in
	 * RB_solve.
	 * 
	 * TODO: fit to time-dependent residual norms for error estimation!
	 */
	@Override
	protected double compute_residual_dual_norm(int N) {
		// This assembly assumes we have already called
		// cache_online_residual_terms
		// and that the RB_solve parameter is constant in time

		RealVector RB_u_euler_theta = RB_solution.mapMultiply(getEulerTheta())
				.add(old_RB_solution.mapMultiply(1. - getEulerTheta()));
		RealVector mass_coeffs = RB_solution.subtract(old_RB_solution)
				.mapMultiply(-1. / getdt());

		double residual_norm_sq = cached_Fq_term;

		for (int i = 0; i < N; i++) {
			residual_norm_sq += RB_u_euler_theta.getEntry(i)
					* cached_Fq_Aq_vector[i];
			residual_norm_sq += mass_coeffs.getEntry(i)
					* cached_Fq_Mq_vector[i];
		}

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++) {
				residual_norm_sq += RB_u_euler_theta.getEntry(i)
						* RB_u_euler_theta.getEntry(j)
						* cached_Aq_Aq_matrix[i][j];
				residual_norm_sq += mass_coeffs.getEntry(i)
						* mass_coeffs.getEntry(j) * cached_Mq_Mq_matrix[i][j];
				residual_norm_sq += RB_u_euler_theta.getEntry(i)
						* mass_coeffs.getEntry(j) * cached_Aq_Mq_matrix[i][j];
			}

		if (residual_norm_sq < 0) {
			Log.d(DEBUG_TAG, "Warning: Square of residual norm is negative "
					+ "in TransientRBSystem::compute_residual_dual_norm()");

			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect result
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}

	protected double conv_weight(double x) {
		// Specify a Gaussian with standard deviation sigma

		double sigma = filter_width * getdt();

		return 1. / Math.sqrt(2. * Math.PI * sigma * sigma)
				* Math.exp(-x * x / (2. * sigma * sigma));
	}

	/**
	 * Gets dt, the time-step size.
	 */
	public double getdt() {
		return dt;
	}

	/**
	 * Get/set euler_theta, parameter that determines the temporal
	 * discretization. euler_theta = 0 ---> Forward Euler euler_theta = 0.5 --->
	 * Crank-Nicolson euler_theta = 1 ---> Backward Euler
	 */
	public double getEulerTheta() {
		return euler_theta;
	}

	/**
	 * return truth solution
	 * 
	 * @see rb.java.RBSystem#getFullSolution()
	 */
	@Override
	public float[][][] getFullSolution() {
		// RB solution size
		int N = timestepRBSolutions[1].getDimension();
		// Number of time-steps to show
		int nt = getVisualNumTimesteps();
		// Dimension of the full solution
		int fullDim = fullBasisVectors[0][0].length;

		float[][][] truth_sol = new float[getNumFields()][1][fullDim * nt];
		for (int fieldNr = 0; fieldNr < getNumFields(); fieldNr++) {
			double tmpval;
			for (timestep = 1; timestep <= nt; timestep++) {
				// Choose equally spaced indices
				int solidx = (int) Math.round(Math.floor(timestep
						* fTotalTimesteps / nt));
				for (int dim = 0; dim < fullDim; dim++) {
					tmpval = 0;
					for (int j = 0; j < N; j++) {
						tmpval += fullBasisVectors[fieldNr][j][dim]
								* timestepRBSolutions[solidx].getEntry(j);
					}
					truth_sol[fieldNr][0][(timestep - 1) * fullDim + dim] = (float) tmpval;
				}
			}
		}
		return truth_sol;
	}

	/**
	 * @return Q_m, the number of terms in the affine expansion of the mass
	 *         matrix
	 */
	public int getQm() {
		return fQm;
	}

	/**
	 * Get K, the total number of time-steps.
	 */
	public int getTotalTimesteps() {
		return fTotalTimesteps;
	}

	/*
	 * public double[] get_RBsolution(int nt){ int N = get_N(); double[] tmpsol
	 * = new double[N*nt]; for (_k = 1; _k <= nt; _k++) for (int j = 0; j < N;
	 * j++) tmpsol[(_k-1)*N+j] =
	 * RB_temporal_solution_data[(int)Math.round(Math.floor
	 * (_k*_K/nt))].getEntry(j); return tmpsol; }
	 */
	public int getVisualNumTimesteps() {
		/*
		 * int nt = Math.round(50000/get_calN()); nt = nt>_K?_K:nt; return nt;
		 */
		int nt = (int) Math.round(75000 / getGeometry().nodes
				/ (1 + 0.4 * (getNumFields() - 1))); // can go up to 150000
//		nt = nt > 150 ? 150 : nt; // cap nt at 25
		nt = nt > 25 ? 25 : nt; // cap nt at 25
		return nt > fTotalTimesteps ? fTotalTimesteps : nt;
	}

	/**
	 * Resize the output vectors according to n_outputs.
	 */
	@Override
	protected void initialize_data_vectors() {
		super.initialize_data_vectors();

		RB_outputs_all_k = new double[getNumOutputs()][getTotalTimesteps() + 1];
		RB_output_error_bounds_all_k = new double[getNumOutputs()][getTotalTimesteps() + 1];

		// Resize the error bound vector
		error_bound_all_k = new double[getTotalTimesteps() + 1];

		// Resize the array that stores the solution data at all time levels
		timestepRBSolutions = new RealVector[getTotalTimesteps() + 1];
	}

	/**
	 * Override read_offline_data_from_files in order to read in the mass matrix
	 * and initial condition data as well.
	 */
	@Override
	public void loadOfflineData_rbappmit(AModelManager m) throws IOException {

		super.loadOfflineData_rbappmit(m);

		// Initialize the residual caching data storage
		cached_Fq_Aq_vector = new double[getNBF()];
		cached_Aq_Aq_matrix = new double[getNBF()][getNBF()];
		cached_Fq_Mq_vector = new double[getNBF()];
		cached_Aq_Mq_matrix = new double[getNBF()][getNBF()];
		cached_Mq_Mq_matrix = new double[getNBF()][getNBF()];

		{
			BufferedReader reader = m.getBufReader("RB_L2_matrix.dat");

			String[] tokens = reader.readLine().split(" ");

			// Set the size of the inner product matrix
			RB_L2_matrix = new Array2DRowRealMatrix(getNBF(), getNBF());

			// Fill the matrix
			int count = 0;
			for (int i = 0; i < getNBF(); i++)
				for (int j = 0; j < getNBF(); j++) {
					RB_L2_matrix.setEntry(i, j,
							Double.parseDouble(tokens[count]));
					count++;
				}
			reader.close();
			reader = null;

			Log.d(DEBUG_TAG, "Finished reading RB_L2_matrix.dat");
		}

		// Read in the M_q matrices
		{
			RB_M_q_matrix = new RealMatrix[getQm()];
			for (int q_m = 0; q_m < getQm(); q_m++) {

				BufferedReader reader = m.getBufReader("RB_M_"
						+ String.format("%03d", q_m) + ".dat");

				String line = reader.readLine();
				reader.close();
				reader = null;
				String[] tokens = line.split(" ");

				// Set the size of the inner product matrix
				RB_M_q_matrix[q_m] = new Array2DRowRealMatrix(getNBF(),
						getNBF());

				// Fill the vector
				int count = 0;
				for (int i = 0; i < getNBF(); i++)
					for (int j = 0; j < getNBF(); j++) {
						RB_M_q_matrix[q_m].setEntry(i, j,
								Double.parseDouble(tokens[count]));
						count++;
					}
			}
			Log.d(DEBUG_TAG, "Finished reading RB_M_q data");
		}

		// Read in Fq_Mq representor norm data
		{
			BufferedReader reader = m.getBufReader("Fq_Mq_norms.dat");

			String line = reader.readLine();
			reader.close();
			reader = null;
			String[] tokens = line.split(" ");

			// Declare the array
			Fq_Mq_representor_norms = new double[getQf()][getQm()][getNBF()];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < getQf(); q_f++)
				for (int q_m = 0; q_m < getQm(); q_m++)
					for (int i = 0; i < getNBF(); i++) {
						Fq_Mq_representor_norms[q_f][q_m][i] = Double
								.parseDouble(tokens[count]);
						count++;
					}

			Log.d(DEBUG_TAG, "Finished reading Fq_Mq_norms.dat");
		}

		// Read in M_M representor norm data
		{
			BufferedReader reader = m.getBufReader("Mq_Mq_norms.dat");

			String line = reader.readLine();
			reader.close();
			reader = null;
			String[] tokens = line.split(" ");

			// Declare the array
			int Q_m_hat = getQm() * (getQm() + 1) / 2;
			Mq_Mq_representor_norms = new double[Q_m_hat][getNBF()][getNBF()];

			// Fill it
			int count = 0;
			for (int q = 0; q < Q_m_hat; q++)
				for (int i = 0; i < getNBF(); i++)
					for (int j = 0; j < getNBF(); j++) {
						Mq_Mq_representor_norms[q][i][j] = Double
								.parseDouble(tokens[count]);
						count++;
					}
			Log.d(DEBUG_TAG, "Finished reading Mq_Mq_norms.dat");
		}

		// Read in Aq_M representor norm data
		{
			try {
				BufferedReader reader = m.getBufReader("Aq_Mq_norms.dat");

				String line = reader.readLine();
				reader.close();
				reader = null;
				String[] tokens = line.split(" ");

				// Declare the array
				Aq_Mq_representor_norms = new double[getQa()][getQm()][getNBF()][getNBF()];

				// Fill it
				int count = 0;
				for (int q_a = 0; q_a < getQa(); q_a++)
					for (int q_m = 0; q_m < getQm(); q_m++)
						for (int i = 0; i < getNBF(); i++)
							for (int j = 0; j < getNBF(); j++) {
								Aq_Mq_representor_norms[q_a][q_m][i][j] = Double
										.parseDouble(tokens[count]);
								count++;
							}

			} catch (IOException iae) {
				// Declare the array
				Aq_Mq_representor_norms = new double[getQa()][getQm()][getNBF()][getNBF()];

				MathObjectReader mr = m.getMathObjReader();
				int count = 0;
				for (int i = 0; i < getQa(); i++)
					for (int j = 0; j < getQm(); j++) {
						String file = "Aq_Mq_" + String.format("%03d", i) + "_"
								+ String.format("%03d", j) + "_norms.bin";
						Aq_Mq_representor_norms[i][j] = mr.readRawDoubleMatrix(
								m.getInStream(file), getNBF(), getNBF());

						// for (int k = 0; k < get_n_basis_functions(); k++)
						// for (int l = 0; l < get_n_basis_functions(); l++)
						// Aq_Mq_representor_norms[i][j][k][l] = dis
						// .ReadDouble();
						count++;
						// dis.close();
					}
			}

			Log.d(DEBUG_TAG, "Finished reading Aq_Mq_norms.dat");
		}
	}

	/**
	 * Override read_offline_data_from_files in order to read in the mass matrix
	 * and initial condition data as well.
	 */
	@Override
	public void loadOfflineDataJRB(AModelManager m) throws IOException {
		super.loadOfflineDataJRB(m);

		// Initialize the residual caching data storage
		cached_Fq_Aq_vector = new double[getNBF()];
		cached_Aq_Aq_matrix = new double[getNBF()][getNBF()];
		cached_Fq_Mq_vector = new double[getNBF()];
		cached_Aq_Mq_matrix = new double[getNBF()][getNBF()];
		cached_Mq_Mq_matrix = new double[getNBF()][getNBF()];

		MathObjectReader mr = m.getMathObjReader();
		/*
		 * Read L2 matrix
		 */
		RB_L2_matrix = mr.readMatrix(m.getInStream("RB_L2_matrix.bin"));

		/*
		 * Read in the M_q matrices
		 */
		String filename;
		RB_M_q_matrix = new RealMatrix[fQm];
		for (int q_m = 0; q_m < fQm; q_m++) {
			filename = "RB_M_" + String.format("%03d", q_m) + ".bin";
			RB_M_q_matrix[q_m] = mr.readMatrix(m.getInStream(filename));
		}
		Log.d(DEBUG_TAG, "Finished reading RB_M_q data");

		/*
		 * Read in Fq_Mq representor norm data
		 */
		// Declare the array
		Fq_Mq_representor_norms = new double[getQf()][fQm][];
		for (int q_f = 0; q_f < getQf(); q_f++) {
			for (int q_m = 0; q_m < fQm; q_m++) {
				filename = "Fq_Mq_" + String.format("%03d", q_f) + "_"
						+ String.format("%03d", q_m) + ".bin";
				Fq_Mq_representor_norms[q_f][q_m] = mr.readRawDoubleVector(m
						.getInStream(filename));
			}
		}
		Log.d(DEBUG_TAG, "Finished reading Fq_Mq_norms.dat");

		/*
		 * Read in Mq_Mq representor norm data
		 */
		int Q_m_hat = fQm * (fQm + 1) / 2;
		Mq_Mq_representor_norms = new double[Q_m_hat][][];
		for (int q1 = 0; q1 < fQm; q1++) {
			for (int q2 = 0; q2 < fQm - q1; q2++) {
				filename = "Mq_Mq_" + String.format("%03d", q1) + "_"
						+ String.format("%03d", q2) + ".bin";
				Mq_Mq_representor_norms[q2 + q1 * fQm] = mr
						.readRawDoubleMatrix(m.getInStream(filename));
			}
		}
		Log.d(DEBUG_TAG, "Finished reading Mq_Mq_norms.dat");

		/*
		 * Read in Aq_M representor norm data
		 */
		Aq_Mq_representor_norms = new double[getQa()][fQm][][];
		for (int i = 0; i < getQa(); i++) {
			for (int j = 0; j < fQm; j++) {
				filename = "Aq_Mq_" + String.format("%03d", i) + "_"
						+ String.format("%03d", j) + "_norms.bin";
				Aq_Mq_representor_norms[i][j] = mr.readRawDoubleMatrix(m
						.getInStream(filename));
			}
		}
		Log.d(DEBUG_TAG, "Finished reading Aq_Mq_norms.dat");
	}

	/**
	 * Perform online solve with the N RB basis functions, for the set of
	 * parameters in current_params, where 1 <= N <= RB_size.
	 */
	@Override
	public double RB_solve(int N) {
		current_N = N;

		if (N > getNBF()) {
			throw new RuntimeException(
					"ERROR: N cannot be larger than the number "
							+ "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException(
					"ERROR: N must be greater than 0 in RB_solve");
		}

		RealMatrix RB_mass_matrix_N = null;
		RealMatrix RB_LHS_matrix = null;
		RealMatrix RB_RHS_matrix = null;

		boolean tdM = ((ITransient) affineFunctionsInstance).isTimeDependentM();
		// Have time-dependent RB LHS/RHS matrices also if only theta^m are
		// time-dependent
		boolean tdA = affineFunctionsInstance.isTimeDependentA() || tdM;

		boolean LHSMatrixIsID = false;
		/*
		 * Initialize constant RB matrices for time-independent theta
		 * coefficient functions. This saves time during the simulations.
		 */
		if (!tdM) {
			// First assemble the mass matrix
			RB_mass_matrix_N = new Array2DRowRealMatrix(N, N);
			for (int q_m = 0; q_m < getQm(); q_m++) {
				RB_mass_matrix_N = RB_mass_matrix_N.add(RB_M_q_matrix[q_m]
						.getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
								thetaQm(q_m)));
			}
		}
		if (!tdA) {
			// No need to copy, "add" returns a new matrix anyways
			RB_LHS_matrix = RB_mass_matrix_N;
			RB_RHS_matrix = RB_mass_matrix_N;
			// RB_LHS_matrix = new Array2DRowRealMatrix(N, N);
			// RB_RHS_matrix = new Array2DRowRealMatrix(N, N);
			//
			// RB_LHS_matrix = RB_LHS_matrix.add(RB_mass_matrix_N
			// .scalarMultiply(1. / getdt()));
			// RB_RHS_matrix = RB_RHS_matrix.add(RB_mass_matrix_N
			// .scalarMultiply(1. / getdt()));

			for (int q_a = 0; q_a < getQa(); q_a++) {
				RB_LHS_matrix = RB_LHS_matrix.add(RB_A_q_vector[q_a]
						.getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
								getEulerTheta() * getdt() * thetaQa(q_a)));
				RB_RHS_matrix = RB_RHS_matrix.add(RB_A_q_vector[q_a]
						.getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
								-(1. - getEulerTheta()) * getdt()
										* thetaQa(q_a)));
			}
			LHSMatrixIsID = isIdentityMatrix(RB_LHS_matrix);
		}

		setTimeStep(0); // Sets the member variable timestep to zero

		// Get initial conditions
		RB_solution = getInitialCoefficients(N);
		RealVector RB_solution_N = RB_solution.copy();
		RealVector old_RB_solution_N = RB_solution_N.copy();

		// Initialize rhs
		RealVector RB_rhs_N = new ArrayRealVector(N);

		// Load the initial condition into RB_temporal_solution_data
		timestepRBSolutions[timestep] = RB_solution_N;

		double error_bound_sum = 0.;

		// Set error bound at _k=0
		error_bound_all_k[timestep] = Math.sqrt(error_bound_sum);

		// Compute the outputs and associated error bounds at _k=0
		for (int k = 0; k < getNumOutputs(); k++) {
			RB_outputs_all_k[k][timestep] = 0.;
			RB_output_error_bounds_all_k[k][timestep] = 0.;
			for (int q_l = 0; q_l < getQl(k); q_l++) {
				RB_outputs_all_k[k][timestep] += thetaQl(k, q_l, 0)
						* (RB_output_vectors[k][q_l].getSubVector(0, N)
								.dotProduct(RB_solution_N));
			}
			RB_output_error_bounds_all_k[k][timestep] = compute_output_dual_norm(
					k, 0) * error_bound_all_k[timestep];
		}

		double alpha_LB = get_SCM_lower_bound();

		// Precompute time-invariant parts of the dual norm of the residual.
		cache_online_residual_terms(N);

		/*
		 * Main time step loop
		 */
		for (int time_level = 1; time_level <= fTotalTimesteps; time_level++) {
			// The current time
			double t = getdt() * time_level;

			/*
			 * Initialize constant RB matrices for time-independent theta
			 * coefficient functions. This saves time during the simulations.
			 */
			if (tdM) {
				// First assemble the mass matrix
				RB_mass_matrix_N = new Array2DRowRealMatrix(N, N);
				for (int q_m = 0; q_m < getQm(); q_m++) {
					RB_mass_matrix_N = RB_mass_matrix_N.add(RB_M_q_matrix[q_m]
							.getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
									thetaQm(q_m, t)));
				}
			}
			if (tdA) {
				// No need to copy, "add" returns a new matrix anyways
				RB_LHS_matrix = RB_mass_matrix_N;
				RB_RHS_matrix = RB_mass_matrix_N;
				// RB_LHS_matrix = new Array2DRowRealMatrix(N, N);
				// RB_RHS_matrix = new Array2DRowRealMatrix(N, N);
				//
				// RB_LHS_matrix = RB_LHS_matrix.add(RB_mass_matrix_N
				// .scalarMultiply(1. / getdt()));
				// RB_RHS_matrix = RB_RHS_matrix.add(RB_mass_matrix_N
				// .scalarMultiply(1. / getdt()));

				for (int q_a = 0; q_a < getQa(); q_a++) {
					RB_LHS_matrix = RB_LHS_matrix
							.add(RB_A_q_vector[q_a].getSubMatrix(0, N - 1, 0,
									N - 1)
									.scalarMultiply(
											getEulerTheta() * getdt()
													* thetaQa(q_a, t)));
					RB_RHS_matrix = RB_RHS_matrix.add(RB_A_q_vector[q_a]
							.getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
									-(1. - getEulerTheta()) * getdt()
											* thetaQa(q_a, t)));
				}
				LHSMatrixIsID = isIdentityMatrix(RB_LHS_matrix);
			}

			setTimeStep(time_level); // This updates the member variable
										// timestep
			old_RB_solution_N = RB_solution_N;

			// Compute RB_rhs, as RB_LHS_matrix x old_RB_solution
			RB_rhs_N = RB_RHS_matrix.operate(old_RB_solution_N);

			// Add forcing terms
			RealVector force = new ArrayRealVector(N);
			for (int q_f = 0; q_f < getQf(); q_f++) {
				force = force.add(RB_F_q_vector[q_f].getSubVector(0, N)
						.mapMultiplyToSelf(thetaQf(q_f, t)));
			}
			RB_rhs_N = RB_rhs_N.add(force.mapMultiplyToSelf(getdt()));

			if (!LHSMatrixIsID) {
				// Solve the linear system
				RB_solution_N = new LUDecompositionImpl(RB_LHS_matrix)
						.getSolver().solve(RB_rhs_N);
			} else {
				RB_solution_N = RB_rhs_N;
			}

			double[] sol = RB_solution_N.getData();
			String sol_str = "[";
			for (int i = 0; i < sol.length; i++) {
				sol_str += String.format("%1.15e  ", sol[i]);
			}
			Log.d("TransientRBSystem",
					"RB_solution at t="
							+ String.format("%5f", time_level * getdt())
							+ ": " + sol_str + "]");

			// Save RB_solution for current time level
			timestepRBSolutions[timestep] = RB_solution_N;

			// Evaluate the dual norm of the residual for RB_solution_vector
			RB_solution = RB_solution_N;
			old_RB_solution = old_RB_solution_N;
			double epsilon_N = compute_residual_dual_norm(N);

			error_bound_sum += residual_scaling_numer(alpha_LB)
					* Math.pow(epsilon_N, 2.);

			// store error bound at time-level _k
			error_bound_all_k[timestep] = Math.sqrt(error_bound_sum
					/ residual_scaling_denom(alpha_LB));

			// Now compute the outputs and associated errors
			for (int i = 0; i < getNumOutputs(); i++) {
				RB_outputs_all_k[i][timestep] = 0.;
				RB_output_error_bounds_all_k[i][timestep] = 0.;
				for (int q_l = 0; q_l < getQl(i); q_l++) {
					RB_outputs_all_k[i][timestep] += thetaQl(i, q_l, t)
							* (RB_output_vectors[i][q_l].getSubVector(0, N)
									.dotProduct(RB_solution_N));
				}
				RB_output_error_bounds_all_k[i][timestep] = compute_output_dual_norm(
						i, t) * error_bound_all_k[timestep];
			}
		}

		// double[] sol = RB_outputs_all_k[0];
		// String sol_str = "[";
		// for (int i = 0; i < sol.length; i++) {
		// sol_str += String.format("%1.5e %1.5e\n", sol[i],
		// RB_output_error_bounds_all_k[0][i]);
		// }
		// Log.d("TransientRBSystem", "RB_outputs_all_k: " + sol_str + "]");

		// Now compute the L2 norm of the RB solution at time-level _K
		// to normalize the error bound
		// We reuse RB_rhs here
		RealMatrix RB_L2_matrix_N = RB_L2_matrix.getSubMatrix(0, N - 1, 0,
				N - 1);
		double final_RB_L2_norm = Math.sqrt(RB_L2_matrix_N.operate(
				RB_solution_N).dotProduct(RB_solution_N));

		if (apply_temporal_filter_flag) {
			apply_temporal_filter();
		}

		return (return_rel_error_bound ? error_bound_all_k[fTotalTimesteps]
				/ final_RB_L2_norm : error_bound_all_k[fTotalTimesteps]);
	}

	private boolean isIdentityMatrix(RealMatrix m) {
		if (!m.isSquare())
			return false;
		for (int i = 0; i < m.getRowDimension(); i++) {
			for (int j = 0; j < m.getColumnDimension(); j++) {
				if ((i == j && m.getEntry(i, j) != 1)
						|| (i != j && m.getEntry(i, j) != 0))
					return false;
			}
		}
		return true;
	}

	@Override
	protected void readConfigurationJRB(AModelManager m) {
		super.readConfigurationJRB(m);
		fQm = ((ITransient) affineFunctionsInstance).getQm();
		dt = Double.parseDouble(m.getModelXMLTagValue("rb_model.timeinfo.dt"));
		euler_theta = Double.parseDouble(m
				.getModelXMLTagValue("rb_model.timeinfo.euler_theta"));
		fTotalTimesteps = Integer.parseInt(m
				.getModelXMLTagValue("rb_model.timeinfo.K"));

		n_plotting_steps = Integer.parseInt(m.getModelXMLTagValue(
				"model.visual.plotSteps", "" + (fTotalTimesteps + 1)));
	}

	/**
	 * 
	 * @see rb.java.RBSystem#readConfigurationRBAppMIT(rb.java.GetPot)
	 */
	@Override
	protected void readConfigurationRBAppMIT(GetPot infile) {
		super.readConfigurationRBAppMIT(infile);

		dt = infile.call("dt", 0.);
		fTotalTimesteps = infile.call("K", 0);
		euler_theta = infile.call("euler_theta", 1.);

		int apply_temporal_filter_flag_in = infile.call(
				"apply_temporal_filter_flag", 0);
		apply_temporal_filter_flag = (apply_temporal_filter_flag_in != 0);

		double filter_width_in = infile.call("filter_width", 2.);
		filter_width = filter_width_in;

		int n_plotting_steps_in = infile.call("n_plotting_steps",
				getTotalTimesteps() + 1);
		n_plotting_steps = n_plotting_steps_in;

		Log.d(DEBUG_TAG, "TransientRBSystem parameters from "
				+ Const.parameters_filename + ":");
		Log.d(DEBUG_TAG, "dt: " + getdt());
		Log.d(DEBUG_TAG, "Number of time steps: " + getTotalTimesteps());
		Log.d(DEBUG_TAG, "euler_theta (for generalized Euler): "
				+ getEulerTheta());
		Log.d(DEBUG_TAG, "Apply a temporal filter? "
				+ apply_temporal_filter_flag);
		if (apply_temporal_filter_flag) {
			Log.d(DEBUG_TAG, "Temporal filter std. dev. " + filter_width);
			Log.d(DEBUG_TAG, "Number of timesteps to be plotted"
					+ n_plotting_steps);
		}

		fQm = ((ITransient) affineFunctionsInstance).getQm();
		Log.d(DEBUG_TAG, "Q_m = " + fQm);
	}

	/**
	 * Specifies the residual scaling on the denominator to be used in the a
	 * posteriori error bound. Overload in subclass in order to obtain the
	 * desired error bound.
	 */
	@Override
	protected double residual_scaling_denom(double alpha_LB) {
		return alpha_LB;
	}

	/**
	 * Specifies the residual scaling on the numerator to be used in the a
	 * posteriori error bound. Overload in subclass in order to obtain the
	 * desired error bound.
	 */
	protected double residual_scaling_numer(double alpha_LB) {
		return getdt();
	}

	// PROTECTED FUNCTIONS

	public void setTimeStep(int k_in) {
		this.timestep = k_in;
	}

	/**
	 * Evaluate theta_q_m (for the q^th mass matrix term) at the current
	 * parameter.
	 * 
	 * @param i
	 * @return
	 */
	public double thetaQm(int i) {
		return thetaQm(i, 0);
	}

	// public void set_dt(double dt_in) {
	// this.dt = dt_in;
	// }
	//
	// public void set_euler_theta(double euler_theta_in) {
	// this.euler_theta = euler_theta_in;
	// }
	//
	// public void set_K(int K_in) {
	// this._K = K_in;
	// }

	/**
	 * Evaluate theta_q_m (for the q^th mass matrix term) at the current
	 * parameter.
	 * 
	 * @param i
	 * @param t
	 * @return
	 */
	public double thetaQm(int i, double t) {
		return ((ITransient) affineFunctionsInstance).thetaQm(i, getParams()
				.getCurrent(), t);
	}

	// /**
	// * Set the secondary SCM system
	// */
	// public void setSecondarySCM(RBSCMSystem second_scm_system) {
	// mSecondRbScmSystem = second_scm_system;
	// }

	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector. This does not assume cached data hence works for
	 * parameters that change as a function of time.
	 */
	protected double uncached_compute_residual_dual_norm(int N) {

		RealVector RB_u_euler_theta = RB_solution.mapMultiply(getEulerTheta())
				.add(old_RB_solution.mapMultiply(1. - getEulerTheta()));
		RealVector mass_coeffs = RB_solution.subtract(old_RB_solution)
				.mapMultiply(-1. / getdt());

		double residual_norm_sq = 0.;

		int q = 0;
		for (int q_f1 = 0; q_f1 < getQf(); q_f1++) {
			double cached_theta_q_f1 = thetaQf(q_f1);
			for (int q_f2 = q_f1; q_f2 < getQf(); q_f2++) {
				double delta = (q_f1 == q_f2) ? 1. : 2.;
				residual_norm_sq += delta * cached_theta_q_f1 * thetaQf(q_f2)
						* Fq_representor_norms[q];

				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			double cached_theta_q_f = thetaQf(q_f);
			for (int q_a = 0; q_a < getQa(); q_a++) {
				double cached_theta_q_a = thetaQa(q_a);
				for (int i = 0; i < N; i++) {
					residual_norm_sq += 2. * RB_u_euler_theta.getEntry(i)
							* cached_theta_q_f * cached_theta_q_a
							* Fq_Aq_representor_norms[q_f][q_a][i];
				}
			}
		}

		q = 0;
		for (int q_a1 = 0; q_a1 < getQa(); q_a1++) {
			double cached_theta_q_a1 = thetaQa(q_a1);
			for (int q_a2 = q_a1; q_a2 < getQa(); q_a2++) {
				double cached_theta_q_a2 = thetaQa(q_a2);
				double delta = (q_a1 == q_a2) ? 1. : 2.;

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						residual_norm_sq += delta
								* RB_u_euler_theta.getEntry(i)
								* RB_u_euler_theta.getEntry(j)
								* cached_theta_q_a1 * cached_theta_q_a2
								* Aq_Aq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

		// Now add the terms due to the time-derivative
		q = 0;
		for (int q_m1 = 0; q_m1 < getQm(); q_m1++) {
			double cached_theta_q_m1 = thetaQm(q_m1);
			for (int q_m2 = q_m1; q_m2 < getQm(); q_m2++) {
				double cached_theta_q_m2 = thetaQm(q_m2);
				double delta = (q_m1 == q_m2) ? 1. : 2.;

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						residual_norm_sq += delta * mass_coeffs.getEntry(i)
								* mass_coeffs.getEntry(j) * cached_theta_q_m1
								* cached_theta_q_m2
								* Mq_Mq_representor_norms[q][i][j];
					}
				}
				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			double cached_theta_q_f = thetaQf(q_f);

			for (int q_m = 0; q_m < getQm(); q_m++) {
				double cached_theta_q_m = thetaQm(q_m);

				for (int i = 0; i < N; i++) {
					residual_norm_sq += 2. * mass_coeffs.getEntry(i)
							* cached_theta_q_f * cached_theta_q_m
							* Fq_Mq_representor_norms[q_f][q_m][i];
				}
			}
		}

		for (int q_a = 0; q_a < getQa(); q_a++) {
			double cached_theta_q_a = thetaQa(q_a);
			for (int q_m = 0; q_m < getQm(); q_m++) {
				double cached_theta_q_m = thetaQm(q_m);

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						residual_norm_sq += 2. * RB_u_euler_theta.getEntry(i)
								* mass_coeffs.getEntry(j) * cached_theta_q_a
								* cached_theta_q_m
								* Aq_Mq_representor_norms[q_a][q_m][i][j];
					}
				}
			}
		}

		if (residual_norm_sq < 0) {
			Log.d(DEBUG_TAG, "Warning: Square of residual norm is negative "
					+ "in TransientRBSystem::compute_residual_dual_norm()");

			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect result
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}
}
