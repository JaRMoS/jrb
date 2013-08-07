package rb;

import jarmos.Log;
import jarmos.io.AModelManager;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

/**
 * This class provides the Online reduced basis functionality for quadratically nonlinear time-dependent problems. This
 * class is modeled on the QNTransientRBSystem class in rbappmit
 * 
 * This class has been taken from the original @ref rbappmit package and modified to fit into the current JaRMoS
 * framework.
 * 
 * @author Daniel Wirtz
 * @date Aug 29, 2011
 * 
 */
public class QNTransientRBSystem extends TransientRBSystem {

	// Logging tag
	private static final String DEBUG_TAG = "QNTransientRBSystem";

	// The tolerance for the Newton solver
	private double nonlinear_tolerance;

	// The maximum number of Newton iterations
	private int n_newton_steps;

	// The nominal lower bound for the stability factor
	private double nominal_rho_min;
	private double nominal_rho_max;

	// The RB data for the trilinear form
	protected double[][][] RB_trilinear_form;

	/**
	 * Vectors storing the residual representor inner products to be used in computing the residuals online.
	 */
	double[][][] Fq_C_representor_norms;
	double[][][][] Mq_C_representor_norms;
	double[][][][] Aq_C_representor_norms;
	double[][][][] C_C_representor_norms;

	/**
	 * Member variable that stores the exponential growth factor of the error bound
	 */
	private double tau_prod_k;

	// Private member that we need in calling the SCM
	private int _N_in_RB_solve;

	/**
	 * Perform online solve for current_params with the N basis functions. Overload this to solve the nonlinear RB
	 * system using Newton's method.
	 * 
	 * TODO: Use the current time information when the theta's are time-dependent, as in TransientRBSystem.
	 */
	@Override
	public double RB_solve(int N) {

		current_N = N;
		_N_in_RB_solve = N;

		// Initialize tau_prod_k
		tau_prod_k = 1.;

		if (N > getNBF()) {
			throw new RuntimeException("ERROR: N cannot be larger than the number " + "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException("ERROR: N must be greater than 0 in RB_solve");
		}

		double dt = getdt();

		// First assemble the mass matrix
		RealMatrix RB_mass_matrix_N = new Array2DRowRealMatrix(N, N);
		for (int q_m = 0; q_m < getQm(); q_m++) {
			RB_mass_matrix_N = RB_mass_matrix_N.add(RB_M_q_matrix[q_m].getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
					thetaQm(q_m)));
		}

		RealMatrix RB_RHS_Aq_matrix = new Array2DRowRealMatrix(N, N);

		RealMatrix Base_RB_LHS_matrix = RB_mass_matrix_N.scalarMultiply(1. / dt);

		for (int q_a = 0; q_a < getQa(); q_a++) {
			double cached_theta_q_a = thetaQa(q_a);
			Base_RB_LHS_matrix = Base_RB_LHS_matrix.add(RB_A_q_vector[q_a].getSubMatrix(0, N - 1, 0, N - 1)
					.scalarMultiply(getEulerTheta() * cached_theta_q_a));
			RB_RHS_Aq_matrix = RB_RHS_Aq_matrix.add(RB_A_q_vector[q_a].getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(
					-cached_theta_q_a));
		}

		// Set system time level to 0
		setTimeStep(0);

		// Initialize a vector to store our current Newton iterate
		RealVector RB_u_bar = new ArrayRealVector(N);

		// and load the _k=0 data
		RB_solution = RB_u_bar;
		timestepRBSolutions[timestep] = RB_u_bar; // Should use .copy()
													// here!

		double error_bound_sum = 0.;

		// Set error bound at _k=0
		error_bound_all_k[timestep] = Math.sqrt(error_bound_sum);

		// Compute the outputs and associated error bounds at _k=0
		for (int i = 0; i < getNumOutputs(); i++) {
			RB_outputs_all_k[i][timestep] = 0.;
			RB_output_error_bounds_all_k[i][timestep] = 0.;
			for (int q_l = 0; q_l < getQl(i); q_l++) {
				RB_outputs_all_k[i][timestep] += thetaQl(i, q_l, 0)
						* (RB_output_vectors[i][q_l].getSubVector(0, N).dotProduct(RB_solution));
			}
			RB_output_error_bounds_all_k[i][timestep] = compute_output_dual_norm(i, 0) * error_bound_all_k[timestep];
		}

		// Initialize a vector to store the solution from the old time-step
		RealVector RB_u_old = new ArrayRealVector(N);

		// Initialize a vector to store the Newton increment, RB_delta_u
		RealVector RB_delta_u = new ArrayRealVector(N);

		// Pre-compute eval_theta_c()
		double cached_theta_c = eval_theta_c();

		for (int time_level = 1; time_level <= fTotalTimesteps; time_level++) {

			double t = time_level * getdt();

			setTimeStep(time_level); // update the member variable _k

			// Set RB_u_old to be the result of the previous Newton loop
			RB_u_old = RB_u_bar.copy();

			// Now we begin the nonlinear loop
			for (int l = 0; l < n_newton_steps; ++l) {
				// Get u_euler_theta = euler_theta*RB_u_bar +
				// (1-euler_theta)*RB_u_old
				RealVector RB_u_euler_theta = RB_u_bar.mapMultiply(getEulerTheta()).add(
						RB_u_old.mapMultiply(1. - getEulerTheta()));

				// Assemble the left-hand side for the RB linear system
				RealMatrix RB_LHS_matrix = Base_RB_LHS_matrix.copy();

				// Add the trilinear term
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						double new_entry = RB_LHS_matrix.getEntry(i, j);
						for (int n = 0; n < N; n++) {
							new_entry += cached_theta_c * getEulerTheta() * RB_u_euler_theta.getEntry(n)
									* (RB_trilinear_form[n][i][j] + RB_trilinear_form[j][i][n]);
						}
						RB_LHS_matrix.setEntry(i, j, new_entry);
					}
				}

				// Assemble the right-hand side for the RB linear system (the
				// residual)
				// First add forcing terms
				RealVector RB_rhs = new ArrayRealVector(N);

				for (int q_f = 0; q_f < getQf(); q_f++) {
					RB_rhs = RB_rhs.add(RB_F_q_vector[q_f].getSubVector(0, N).mapMultiply(thetaQf(q_f)));
				}

				// Now add -1./dt * M * (RB_u_bar - RB_u_old)
				RB_rhs = RB_rhs.add(RB_mass_matrix_N.operate(RB_u_bar).mapMultiply(-1. / dt));
				RB_rhs = RB_rhs.add(RB_mass_matrix_N.operate(RB_u_old).mapMultiply(1. / dt));

				// Now add the Aq stuff
				RB_rhs = RB_rhs.add(RB_RHS_Aq_matrix.operate(RB_u_euler_theta));

				// Finally add the trilinear term
				for (int i = 0; i < N; i++) {
					double new_entry = RB_rhs.getEntry(i);

					for (int j = 0; j < N; j++) {
						double RB_u_euler_theta_j = RB_u_euler_theta.getEntry(j);

						for (int n = 0; n < N; n++) {
							new_entry -= cached_theta_c * RB_u_euler_theta.getEntry(n) * RB_u_euler_theta_j
									* RB_trilinear_form[n][i][j];
						}
					}
					RB_rhs.setEntry(i, new_entry);
				}

				DecompositionSolver solver = new LUDecompositionImpl(RB_LHS_matrix).getSolver();
				RB_delta_u = solver.solve(RB_rhs);

				// update the Newton iterate
				RB_u_bar = RB_u_bar.add(RB_delta_u);

				// Compute the l2 norm of RB_delta_u
				double RB_delta_u_norm = RB_delta_u.getNorm();

				if (RB_delta_u_norm < nonlinear_tolerance) {
					break;
				}

				if ((l == (n_newton_steps - 1)) && (RB_delta_u_norm > nonlinear_tolerance)) {
					throw new RuntimeException("RB Newton loop did not converge");
				}
			}

			// Load RB_solution into RB_solution_vector for residual computation
			RB_solution = RB_u_bar;
			old_RB_solution = RB_u_old;
			timestepRBSolutions[timestep] = RB_u_bar; // should use copy
														// here!

			double rho_LB = (mRbScmSystem == null) ? get_nominal_rho_LB() : get_SCM_lower_bound();

			// Evaluate the dual norm of the residual for RB_solution_vector
			double epsilon_N = compute_residual_dual_norm(N);

			error_bound_sum += residual_scaling_numer(rho_LB) * Math.pow(epsilon_N, 2.);

			// store error bound at time-level _k
			error_bound_all_k[timestep] = Math.sqrt(error_bound_sum / residual_scaling_denom(rho_LB));

			// Now compute the outputs and associated error bounds
			for (int i = 0; i < getNumOutputs(); i++) {
				RB_outputs_all_k[i][timestep] = 0.;
				RB_output_error_bounds_all_k[i][timestep] = 0.;
				for (int q_l = 0; q_l < getQl(i); q_l++) {
					RB_outputs_all_k[i][timestep] += thetaQl(i, q_l)
							* (RB_output_vectors[i][q_l].getSubVector(0, N).dotProduct(RB_solution));
				}
				RB_output_error_bounds_all_k[i][timestep] = compute_output_dual_norm(i, t)
						* error_bound_all_k[timestep];
			}
			Log.d(DEBUG_TAG, "output = " + RB_outputs_all_k[0][timestep] + ", bound="
					+ RB_output_error_bounds_all_k[0][timestep]);
		}

		// Now compute the L2 norm of the RB solution at time-level _K
		// to normalize the error bound
		// We reuse RB_rhs here
		double final_RB_L2_norm = Math.sqrt(RB_mass_matrix_N.operate(RB_solution).dotProduct(RB_solution));

		return (return_rel_error_bound ? error_bound_all_k[fTotalTimesteps] / final_RB_L2_norm
				: error_bound_all_k[fTotalTimesteps]);
	}

	/**
	 * Override appropriate for quadratically nonlinear error bounds.
	 */
	@Override
	protected double residual_scaling_numer(double rho_LB) {
		double tau_LB = (0.5 * rho_LB < 0.) ? 0.5 * rho_LB : 0.;

		return getdt() * tau_prod_k / (1. - tau_LB * getdt());
	}

	/**
	 * Override appropriate for quadratically nonlinear error bounds.
	 */
	@Override
	protected double residual_scaling_denom(double rho_LB) {
		double tau_LB = (0.5 * rho_LB < 0.) ? 0.5 * rho_LB : 0.;

		// Update tau_prod_k
		tau_prod_k *= (1. + tau_LB * getdt()) / (1. - tau_LB * getdt());

		// and return it
		return tau_prod_k;
	}

	/**
	 * Set the nonlinear tolerance for Newton's method for both the truth and RB solves.
	 */
	public void set_nonlinear_tolerance(double nonlinear_tolerance_in) {
		nonlinear_tolerance = nonlinear_tolerance_in;
	}

	/**
	 * Get the nonlinear tolerance for Newton's method.
	 */
	public double get_nonlinear_tolerance() {
		return nonlinear_tolerance;
	}

	/**
	 * Set the maximum number of Newton steps for both the truth and RB solves.
	 */
	public void set_n_newton_steps(int n_newton_steps_in) {
		n_newton_steps = n_newton_steps_in;
	}

	/**
	 * Set the maximum number of Newton steps for both the truth and RB solves.
	 */
	public int get_n_newton_steps() {
		return n_newton_steps;
	}

	/**
	 * Evaluate theta_c (for the quadratic nonlinearity) at the current parameter.
	 */
	public double eval_theta_c() {

		Method meth;

		try {
			Class<?> partypes[] = new Class[1];
			partypes[0] = double[].class;

			meth = oldAffFcnCl.getMethod("evaluateC", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateC failed", nsme);
		}

		Double theta_val;
		try {
			Object arglist[] = new Object[1];
			arglist[0] = getParams().getCurrent();

			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			theta_val = (Double) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return theta_val.doubleValue();
	}

	/**
	 * Compute the dual norm of the residual for the solution saved in RB_solution_vector.
	 */
	@Override
	protected double compute_residual_dual_norm(int N) {
		// Use the stored representor inner product values
		// to evaluate the residual dual norm
		double dt = getdt();

		double residual_norm_sq = 0.;

		// Use TransientRBSystem to compute all the linear terms
		residual_norm_sq += Math.pow(super.uncached_compute_residual_dual_norm(N), 2.);

		// Now just need to add the terms involving the nonlinearity
		RealVector RB_u_euler_theta = RB_solution.mapMultiply(getEulerTheta()).add(
				old_RB_solution.mapMultiply(1. - getEulerTheta()));
		RealVector mass_coeffs = RB_solution.subtract(old_RB_solution).mapMultiply(-1. / dt);

		// Pre-compute eval_theta_c()
		double cached_theta_c = eval_theta_c();

		// All residual terms can be treated as positive quantities...
		for (int q_f = 0; q_f < getQf(); q_f++) {
			double cached_theta_q_f = thetaQf(q_f);
			for (int n1 = 0; n1 < N; n1++) {
				for (int j1 = 0; j1 < N; j1++) {
					residual_norm_sq += 2. * cached_theta_q_f * cached_theta_c * RB_u_euler_theta.getEntry(n1)
							* RB_u_euler_theta.getEntry(j1) * Fq_C_representor_norms[q_f][n1][j1];
				}
			}
		}

		for (int q_m = 0; q_m < getQm(); q_m++) {
			double cached_theta_q_m = thetaQm(q_m);
			for (int i = 0; i < N; i++) {
				for (int n1 = 0; n1 < N; n1++) {
					for (int j1 = 0; j1 < N; j1++) {
						residual_norm_sq += 2. * cached_theta_q_m * cached_theta_c * mass_coeffs.getEntry(i)
								* RB_u_euler_theta.getEntry(n1) * RB_u_euler_theta.getEntry(j1)
								* Mq_C_representor_norms[q_m][i][n1][j1];
					}
				}
			}
		}

		for (int q_a = 0; q_a < getQa(); q_a++) {
			double cached_theta_q_a = thetaQa(q_a);
			for (int i = 0; i < N; i++) {
				for (int n1 = 0; n1 < N; n1++) {
					for (int j1 = 0; j1 < N; j1++) {
						residual_norm_sq += 2. * cached_theta_q_a * cached_theta_c * RB_u_euler_theta.getEntry(i)
								* RB_u_euler_theta.getEntry(n1) * RB_u_euler_theta.getEntry(j1)
								* Aq_C_representor_norms[q_a][i][n1][j1];
					}
				}
			}
		}

		for (int n1 = 0; n1 < N; n1++) {
			for (int j1 = 0; j1 < N; j1++) {
				double RB_u_euler_theta_1 = RB_u_euler_theta.getEntry(n1) * RB_u_euler_theta.getEntry(j1);

				for (int n2 = n1; n2 < N; n2++) {
					int init_j2_index = (n2 == n1) ? j1 : 0;
					for (int j2 = init_j2_index; j2 < N; j2++) {
						double RB_u_euler_theta_2 = RB_u_euler_theta.getEntry(n2) * RB_u_euler_theta.getEntry(j2);

						double delta = ((n2 == n1) && (j2 == j1)) ? 1. : 2.;

						residual_norm_sq += delta * cached_theta_c * cached_theta_c * RB_u_euler_theta_1
								* RB_u_euler_theta_2 * C_C_representor_norms[n1][j1][n2][j2];
					}
				}

			}
		}

		if (residual_norm_sq < 0) {
			// Sometimes this is negative due to rounding error,
			// but this error shouldn't affect the error bound
			// too much...
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}

	// Get/set the nominal rho min/max values, we typically read these
	// in from the input file
	public void set_nominal_rho_min(double nominal_rho_LB_in) {
		nominal_rho_min = nominal_rho_LB_in;
	}

	public void set_nominal_rho_max(double nominal_rho_LB_in) {
		nominal_rho_max = nominal_rho_LB_in;
	}

	public double get_nominal_rho_min() {
		return nominal_rho_min;
	}

	public double get_nominal_rho_max() {
		return nominal_rho_max;
	}

	/**
	 * Get the nominal stability factor lower bound. By default, this is a linear function of parameter 0.
	 */
	public double get_nominal_rho_LB() {
		double mu_min = getParams().getMinValue(0);
		double mu_max = getParams().getMaxValue(0);
		double current_mu = getParams().getCurrent()[0];
		return (nominal_rho_max * (mu_max - current_mu) + nominal_rho_min * (current_mu - mu_min)) / (mu_max - mu_min);
	}

	/**
	 * @return the SCM lower bound for current_parameters
	 */
	public double get_SCM_lower_bound() {

		if (mRbScmSystem != null) {
			// Cast to a QNTransientSCMSystem
			QNTransientSCMSystem qnScmSystem = (QNTransientSCMSystem) mRbScmSystem;

			// Tell the SCM system the number of basis functions
			qnScmSystem.set_n_basis_functions(_N_in_RB_solve);

			// Create a parameter vector in which the current time-level
			// is appended to current_parameters.
			int np = getParams().getNumParams();
			double[] params = new double[np + 1];
			for (int i = 0; i < np; i++) {
				params[i] = getParams().getCurrent()[i];
			}
			params[np] = timestep;

			// Set the parameter
			qnScmSystem.sys.getParams().setCurrent(params);

			// Also, construct a vector storing the RB coefficients
			RealVector RB_u_euler_theta = RB_solution.mapMultiply(getEulerTheta()).add(
					old_RB_solution.mapMultiply(1. - getEulerTheta()));

			// Pass params and RB_u_euler_theta to the associated SCM system
			qnScmSystem.set_current_RB_coeffs(RB_u_euler_theta);

			return mRbScmSystem.get_SCM_LB();
		} else {
			return get_SCM_from_AffineFunction();
		}
	}

	/**
	 * 
	 * @see rb.TransientRBSystem#readConfiguration(jarmos.io.AModelManager)
	 */
	@Override
	public void readConfigurationRBAppMIT(GetPot infile) {
		super.readConfigurationRBAppMIT(infile);

		double nonlinear_tolerance_in = infile.call("nonlinear_tolerance", 1.e-8);
		set_nonlinear_tolerance(nonlinear_tolerance_in);

		int n_newton_steps_in = infile.call("n_newton_steps", 15);
		set_n_newton_steps(n_newton_steps_in);

		double nominal_rho_min_in = infile.call("nominal_rho_min", 0);
		set_nominal_rho_min(nominal_rho_min_in);

		double nominal_rho_max_in = infile.call("nominal_rho_max", 0);
		set_nominal_rho_max(nominal_rho_max_in);

		Log.d(DEBUG_TAG, "QNTransientRBSystem parameters from " + Const.parameters_filename + ":");
		Log.d(DEBUG_TAG, "Tolerance for Newton's method: " + get_nonlinear_tolerance());
		Log.d(DEBUG_TAG, "Max number of Newton steps: " + get_n_newton_steps());
		Log.d(DEBUG_TAG, "Nominal rho min: " + get_nominal_rho_min());
		Log.d(DEBUG_TAG, "Nominal rho max: " + get_nominal_rho_max());
	}

	/**
	 * Override read_offline_data_from_files in order to read in the data for the nonlinear problem.
	 */
	@Override
	public void loadOfflineData_rbappmit(AModelManager m) throws IOException {

		super.loadOfflineData_rbappmit(m);

		int n_bfs = getNBF();

		// Read in the trlinear form

		{
			Log.d(DEBUG_TAG, "Starting read RB_trilinear_form.dat");

			String[] tokens;
			{
				BufferedReader reader = m.getBufReader("RB_trilinear_form.dat");

				String line = reader.readLine();
				reader.close();
				tokens = line.split(" ");
			}

			// Set the size of the inner product matrix
			RB_trilinear_form = new double[n_bfs][n_bfs][n_bfs];

			// Fill the array
			int count = 0;
			for (int i = 0; i < n_bfs; i++) {
				for (int j = 0; j < n_bfs; j++) {
					for (int l = 0; l < n_bfs; l++) {
						RB_trilinear_form[i][j][l] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			}

			Log.d(DEBUG_TAG, "Finished reading RB_trilinear_form.dat");
		}

		// Read in Fq_C representor norm data
		{
			String[] tokens;

			{
				BufferedReader reader = m.getBufReader("Fq_C_norms.dat");
				String line = reader.readLine();
				reader.close();
				tokens = line.split(" ");
			}

			// Declare the array
			Fq_C_representor_norms = new double[getQf()][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < getQf(); q_f++)
				for (int i = 0; i < n_bfs; i++) {
					for (int j = 0; j < n_bfs; j++) {
						Fq_C_representor_norms[q_f][i][j] = Double.parseDouble(tokens[count]);
						count++;
					}
				}

			Log.d(DEBUG_TAG, "Finished reading Fq_C_norms.dat");
		}

		// Read in M_M representor norm data
		{
			String[] tokens;
			{
				BufferedReader reader = m.getBufReader("Mq_C_norms.dat");

				tokens = reader.readLine().split(" ");
				reader.close();
			}

			// Declare the array
			Mq_C_representor_norms = new double[getQm()][n_bfs][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int q_m = 0; q_m < getQm(); q_m++)
				for (int i = 0; i < n_bfs; i++)
					for (int j = 0; j < n_bfs; j++) {
						for (int k = 0; k < n_bfs; k++) {
							Mq_C_representor_norms[q_m][i][j][k] = Double.parseDouble(tokens[count]);
							count++;
						}
					}

			Log.d(DEBUG_TAG, "Finished reading Mq_C_norms.dat");
		}

		// Read in Aq_C representor norm data
		{
			BufferedReader reader = m.getBufReader("Aq_C_norms.dat");

			String[] tokens = reader.readLine().split(" ");
			reader.close();
			reader = null;

			// Declare the array
			Aq_C_representor_norms = new double[getQa()][n_bfs][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int q_a = 0; q_a < getQa(); q_a++) {
				for (int i = 0; i < n_bfs; i++) {
					for (int n1 = 0; n1 < n_bfs; n1++) {
						for (int j1 = 0; j1 < n_bfs; j1++) {
							Aq_C_representor_norms[q_a][i][n1][j1] = Double.parseDouble(tokens[count]);
							count++;
						}
					}
				}
			}
			Log.d(DEBUG_TAG, "Finished reading Aq_C_norms.dat");
		}

		// Read in C_C representor norm data
		{
			BufferedReader reader = m.getBufReader("C_C_norms.dat");

			String[] tokens = reader.readLine().split(" ");
			reader.close();
			reader = null;

			// Declare the array
			C_C_representor_norms = new double[n_bfs][n_bfs][n_bfs][n_bfs];

			// Fill it
			int count = 0;
			for (int ii = 0; ii < n_bfs; ii++) {
				for (int i = 0; i < n_bfs; i++) {
					for (int n1 = 0; n1 < n_bfs; n1++) {
						for (int j1 = 0; j1 < n_bfs; j1++) {
							C_C_representor_norms[ii][i][n1][j1] = Double.parseDouble(tokens[count]);
							count++;
						}
					}
				}
			}
			Log.d(DEBUG_TAG, "Finished reading C_C_norms.dat");
		}
	}
}
