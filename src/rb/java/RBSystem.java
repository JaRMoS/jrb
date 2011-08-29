package rb.java;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import rb.java.affinefcn.IWithuL;
import rmcommon.Log;
import rmcommon.io.AModelManager;
import rmcommon.io.MathObjectReader;

/**
 * This class provides the Online stage for the reduced basis method for
 * elliptic steady state problems.
 * 
 * This class is modeled on RBSystem from rbOOmit
 * 
 * Changes by
 * 
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 * 
 */
public class RBSystem extends RBBase {

	// Logging tag
	private static final String DEBUG_TAG = "RBSystem";

	// PUBLIC MEMBER VARIABLES

	protected double[][][] Aq_Aq_representor_norms;

	/*
	 * The number of nodes. TODO: This value should be obtainable via a
	 * reference to the system's GeometryData object after it has been loaded..
	 */
	private int calN;

	// current N
	protected int current_N;

	// The number of output functionals
	private int fNumOutputs;
	// The affine function expansion size for all outputs
	private int[] Ql_values;

	// PRIVATE MEMBER VARIABLES

	protected double[][][] Fq_Aq_representor_norms;

	/**
	 * Arrays storing the residual representor inner products to be used in
	 * computing the residuals in the Online stage. These are resized by reading
	 * in the Offline data written out by rbOOmit.
	 */
	protected double[] Fq_representor_norms;

	/* The number of fields: 0 means no visualization */
	private int fNumFields;

	/**
	 * The number of terms in the affine expansion of the rhs
	 */
	private int fQf;

	/**
	 * The number of uL functions.
	 * 
	 * Zero per default.
	 */
	private int fQuL = 0;

	/**
	 * A reference to the SCM system.
	 */
	public RBSCMSystem mRbScmSystem;

	/**
	 * The number of basis functions in the RB space.
	 */
	private int n_bfs;

	/**
	 * This array stores the dual norms for each output. Row n stores the Q_l
	 * dual norms for the expansion of the n^th output.
	 */
	public double[][] output_dual_norms;

	/**
	 * Dense matrices for the RB computations.
	 */
	protected RealMatrix[] RB_A_q_vector;

	/**
	 * Dense vector for the RHS.
	 */
	protected RealVector[] RB_F_q_vector;

	/**
	 * The vector storing the RB error bounds in the steady-state case.
	 */
	public double[] RB_output_error_bounds;

	/**
	 * The output error bounds at all time levels.
	 */
	public double[][] RB_output_error_bounds_all_k;

	/**
	 * The vectors storing the RB output vectors
	 */
	protected RealVector[][] RB_output_vectors;

	/**
	 * The vector storing the RB output values in the steady-state case.
	 */
	public double[] RB_outputs;
	// (those two below are moved from TransientRBSystem
	/**
	 * The outputs at all time levels.
	 */
	public double[][] RB_outputs_all_k;
	/**
	 * The RB solution vector. Stored as a Vector so that we can easily resize
	 * it during an RB_solve.
	 */
	protected RealVector RB_solution;

	protected double[][][] RB_sweep_solution;
	/**
	 * Boolean flag to indicate whether RB_solve returns an absolute or relative
	 * error bound. True => relative, false => absolute.
	 */
	public boolean return_rel_error_bound;

	protected float[][] uL_vector;

	// PUBLIC FUNCTIONS

	protected float[][][] Z_vector;

	/**
	 * Constructor.
	 */
	public RBSystem() {

		// Initialize n_bfs to 0
		n_bfs = 0;
	}

	/**
	 * Compute the dual norm of the i^th output function at the current
	 * parameter value
	 */
	protected double compute_output_dual_norm(int i, double t) {

		// Use the stored representor inner product values
		// to evaluate the output dual norm
		double output_norm_sq = 0.;

		int q = 0;
		for (int q_l1 = 0; q_l1 < getQl(i); q_l1++) {
			for (int q_l2 = q_l1; q_l2 < getQl(i); q_l2++) {
				double delta = (q_l1 == q_l2) ? 1. : 2.;
				output_norm_sq += delta * eval_theta_q_l(i, q_l1, t)
						* eval_theta_q_l(i, q_l2, t) * output_dual_norms[i][q];
				q++;
			}
		}

		return Math.sqrt(output_norm_sq);
	}
	
	/**
	 * 
	 * @param i
	 * @return
	 */
	public double eval_theta_q_f(int i) {
		return eval_theta_q_f(i, 0);
	}
	
	/**
	 * 
	 * @param i
	 * @param t
	 * @return
	 */
	public double eval_theta_q_f(int i, double t) {
		return affineFunctionsInstance.thetaQf(i, getParams().getCurrent(), t);
	}
	
	/**
	 * 
	 * @param k 
	 * @param i
	 * @return
	 */
	public double eval_theta_q_l(int k, int i) {
		return eval_theta_q_l(k, i, 0);
	}
	
	/**
	 * 
	 * @param k 
	 * @param i
	 * @param t
	 * @return
	 */
	public double eval_theta_q_l(int k, int i, double t) {
		return affineFunctionsInstance.thetaQl(k, i, getParams().getCurrent(), t);
	}

	/**
	 * Compute the dual norm of the residual for the solution saved in
	 * RB_solution_vector.
	 */
	protected double compute_residual_dual_norm(int N) {

		// Use the stored representor inner product values
		// to evaluate the residual norm
		double residual_norm_sq = 0.;

		int q = 0;
		for (int q_f1 = 0; q_f1 < getQf(); q_f1++) {
			for (int q_f2 = q_f1; q_f2 < getQf(); q_f2++) {
				double delta = (q_f1 == q_f2) ? 1. : 2.;
				residual_norm_sq += delta * eval_theta_q_f(q_f1)
						* eval_theta_q_f(q_f2) * Fq_representor_norms[q];

				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			for (int q_a = 0; q_a < getQa(); q_a++) {
				for (int i = 0; i < N; i++) {
					double delta = 2.;
					residual_norm_sq += get_soln_coeff(i) * delta
							* eval_theta_q_f(q_f) * eval_theta_q_a(q_a)
							* Fq_Aq_representor_norms[q_f][q_a][i];
				}
			}
		}

		q = 0;
		for (int q_a1 = 0; q_a1 < getQa(); q_a1++) {
			for (int q_a2 = q_a1; q_a2 < getQa(); q_a2++) {
				double delta = (q_a1 == q_a2) ? 1. : 2.;

				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						residual_norm_sq += get_soln_coeff(i)
								* get_soln_coeff(j) * delta
								* eval_theta_q_a(q_a1) * eval_theta_q_a(q_a2)
								* Aq_Aq_representor_norms[q][i][j];
					}
				}

				q++;
			}
		}

		if (residual_norm_sq < 0.) {
			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect error bound much...
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}

	// Return calN
	public int get_calN() {
		return calN;
	}

//	public double get_dt() {
//		return 0.;
//	}

	public int get_K() {
		return 1;
	}

	/**
	 * Number of output value fields.
	 * 
	 * Zero for no visualization.
	 * @return The number of output fields
	 */
	public int getNumFields() {
		return fNumFields;
	}

	public int get_N() {
		return current_N;
	}

	public int get_nt() {
		return 1;
	}

	/**
	 * TODO: if all affine_functions implement the IAffineFunctions interface,
	 * just call the getQf method of the local instance.
	 * 
	 * @return Q_f, the number of term in the affine expansion of the right-hand
	 *         side
	 */
	public int getQf() {
		return fQf;
	}

	/**
	 * TODO: if all affine_functions implement the IAffineFunctions interface,
	 * just call the getQl method of the local instance.
	 * 
	 * @param output_index
	 *            The index of the output we are interested in
	 * @return the number of terms in the affine expansion of the specified
	 *         output
	 */
	protected int getQl(int output_index) {
		return Ql_values[output_index];
	}

	public int getQuL() {
		return fQuL;
	}

	public double get_RB_output(int n_output, boolean Rpart) {
		return RB_outputs[n_output];
	}

	public double get_RB_output_error_bound(int n_output, boolean Rpart) {
		return RB_output_error_bounds[n_output];
	}

	public double[][] get_RBsolution() {
		double[][] RBsol = new double[1][];
		RBsol[0] = RB_solution.toArray();
		return RBsol;
	}

	/**
	 * A private helper function to get the SCM from AffineFunctions in the case
	 * that SCM_TYPE = NONE
	 */
	protected double get_SCM_from_AffineFunction() {
		// we assume that an SCM LB function has been specified
		// in AffineFunctions.jar
		Method meth;

		try {
			Class<?> partypes[] = new Class[1];
			partypes[0] = double[].class;

			meth = oldAffFcnCl.getMethod("get_SCM_LB", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_SCM_LB failed", nsme);
		}

		Double SCM_val;
		try {
			Object arglist[] = new Object[1];
			arglist[0] = getParams().getCurrent();

			Object SCM_obj = meth.invoke(oldAffFcnObj, arglist);
			SCM_val = (Double) SCM_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return SCM_val.doubleValue();
	}

	/**
	 * @return the SCM lower bound for current_parameters
	 */
	public double get_SCM_lower_bound() {
		if (mRbScmSystem != null) {
			mRbScmSystem.setParams(getParams());
			return mRbScmSystem.get_SCM_LB();
		} else {
			return get_SCM_from_AffineFunction();
		}
	}

	/**
	 * @return the SCM upper bound for current_parameters
	 */
	double get_SCM_upper_bound() {

		if (mRbScmSystem != null) {
			mRbScmSystem.setParams(getParams());
			return mRbScmSystem.get_SCM_UB();
		} else {
			return get_SCM_from_AffineFunction();
		}
	}

	/**
	 * @return coefficient i of RB_solution from the most recent RB_solve.
	 */
	double get_soln_coeff(int i) {
		return RB_solution.getEntry(i);
	}

	public float[][][] get_sweep_truth_sol() {
		int N = RB_sweep_solution[0][0].length;
		int numSweep = RB_sweep_solution.length;
		float[][][] truth_sol = new float[getNumFields()][1][calN * numSweep];
		for (int ifn = 0; ifn < getNumFields(); ifn++) {
			double tmpval;
			for (int iSweep = 0; iSweep < numSweep; iSweep++)
				for (int i = 0; i < calN; i++) {
					tmpval = 0;
					for (int j = 0; j < N; j++)
						tmpval += Z_vector[ifn][j][i]
								* RB_sweep_solution[iSweep][0][j];
					truth_sol[ifn][0][iSweep * calN + i] = (float) tmpval;
				}
		}
		return truth_sol;
	}

	/**
	 * Returns a float array of transformation data for each node.
	 * 
	 * Only to be called for models who have a parameterized geometry.
	 * 
	 * The get_local_transformation method of the AffineFunctions class is
	 * called (using the current parameter \mu) and the linear transform
	 * function returned.
	 * 
	 * @return
	 */
	public float[][] get_tranformation_data() {
		Method meth = null;
		boolean isOK = true;
		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[1];
			partypes[0] = double[].class;

			meth = oldAffFcnCl.getMethod("get_local_transformation", partypes);
		} catch (NoSuchMethodException nsme) {
			// throw new RuntimeException("getMethod for evaluateF failed",
			// nsme);
			isOK = false;
		}

		float[][] T_vector;
		if (isOK) {
			try {
				Object arglist[] = new Object[1];
				arglist[0] = getParams().getCurrent();

				Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
				T_vector = (float[][]) theta_obj;
			} catch (IllegalAccessException iae) {
				throw new RuntimeException(iae);
			} catch (InvocationTargetException ite) {
				throw new RuntimeException(ite.getCause());
			}
		} else {
			/*
			 * Return a fake transformation vector if method invocation went
			 * wrong. (Should better throw an exception here?)
			 */
			T_vector = new float[1][12];
			T_vector[0][0] = 1f;
			T_vector[0][1] = 0f;
			T_vector[0][2] = 0f;
			T_vector[0][3] = 0f;
			T_vector[0][4] = 1f;
			T_vector[0][5] = 0f;
			T_vector[0][6] = 0f;
			T_vector[0][7] = 0f;
			T_vector[0][8] = 1f;
			T_vector[0][9] = 0f;
			T_vector[0][10] = 0f;
			T_vector[0][11] = 0f;
		}

		return T_vector;
	}

	/**
	 * Constructs the true solution from the most recent RB reduced solution
	 * using the full Z data and current coefficients.
	 * 
	 * The second dimension stands for real and complex data, so
	 * sol[field_nr][1] is the complex part of field field_nr.
	 * 
	 * @see ComplexRBSystem for an overridden version of this method.
	 * 
	 * @return a float array with the true solution, with entries
	 *         [field_nr][0][node_nr]
	 */
	public float[][][] get_truth_sol() {
		int N = RB_solution.getDimension();
		float[][][] truth_sol = new float[getNumFields()][1][calN];
		/*
		 * Assign solutions of each field
		 */
		for (int ifn = 0; ifn < getNumFields(); ifn++) {
			double tmpval;
			/*
			 * i is the node number of the grid
			 */
			for (int i = 0; i < calN; i++) {
				tmpval = 0;
				/*
				 * j is the current RB dimension (with coeffs)
				 */
				for (int j = 0; j < N; j++)
					tmpval += Z_vector[ifn][j][i] * get_soln_coeff(j);
				truth_sol[ifn][0][i] = (float) tmpval;
			}
		}
		return truth_sol;
	}

	/**
	 * @return The number of basis functions in the system.
	 */
	public int getNBF() {
		return n_bfs;
	}

	/**
	 * @return the number of output functionals
	 */
	public int getNumOutputs() {
		return fNumOutputs;
	}

	/**
	 * Resize the vectors that store solution data and output data.
	 */
	protected void initialize_data_vectors() {
		// Also, resize RB_outputs and RB_output_error_error_bounds arrays
		RB_outputs = new double[getNumOutputs()];
		RB_output_error_bounds = new double[getNumOutputs()];
	}

	/**
	 * Loads the offline data for the RBSystem
	 * 
	 * @param m
	 * @throws IOException
	 */
	public final void loadOfflineData(AModelManager m) throws IOException {
		if ("rbappmit".equals(m.getModelType())) {
			loadOfflineData_rbappmit(m);
		} else {
			loadOfflineDataJRB(m);
		}
		initialize_data_vectors();
	}

	protected void loadOfflineDataJRB(AModelManager m) throws IOException {

		MathObjectReader mr = m.getMathObjReader();
		String filename;

		/*
		 * Get output dual norms
		 */
		RB_output_vectors = new RealVector[fNumOutputs][];
		output_dual_norms = new double[fNumOutputs][];
		for (int i = 0; i < fNumOutputs; i++) {
			filename = "output_" + String.format("%03d", i) + "_dual_norms.bin";

			output_dual_norms[i] = mr.readRawDoubleVector(m.getInStream(filename));
			// int Q_l_hat = getQl(i) * (getQl(i) + 1) / 2;
			// output_dual_norms[i] = new double[Q_l_hat];
			// for (int q = 0; q < Q_l_hat; q++) {
			// output_dual_norms[i][q] =
			// Double.parseDouble(dual_norms_tokens[q]);
			// }

			RB_output_vectors[i] = new RealVector[Ql_values[i]];
			for (int q_l = 0; q_l < Ql_values[i]; q_l++) {
				filename = "output_" + String.format("%03d", i) + "_"
						+ String.format("%03d", q_l) + ".bin";

				RB_output_vectors[i][q_l] = mr.readVector(m.getInStream(filename));
			}
		}
		Log.d(DEBUG_TAG, "Finished reading output data");

		/*
		 * Read in the F_q vectors
		 */
		RB_F_q_vector = new RealVector[fQf];
		for (int q_f = 0; q_f < fQf; q_f++) {
			filename = "RB_F_" + String.format("%03d", q_f) + ".bin";
			RB_F_q_vector[q_f] = mr.readVector(m.getInStream(filename));
		}
		Log.d(DEBUG_TAG, "Finished reading RB_F_q data");

		/*
		 * Read in the A_q matrices
		 */
		RB_A_q_vector = new RealMatrix[getQa()];
		for (int q_a = 0; q_a < getQa(); q_a++) {
			filename = "RB_A_" + String.format("%03d", q_a) + ".bin";
			RB_A_q_vector[q_a] = mr.readMatrix(m.getInStream(filename));
		}
		Log.d(DEBUG_TAG, "Finished reading RB_A_q data");

		/*
		 * Read in F_q representor norm data Contains the upper triangular
		 * entries of the pairwise norm matrix
		 */
		Fq_representor_norms = mr.readRawDoubleVector(m.getInStream("Fq_norms.bin"));

		// Read in Fq_Aq representor norm data

		// Declare the array
		Fq_Aq_representor_norms = new double[fQf][getQa()][n_bfs];
		for (int q_a = 0; q_a < getQa(); q_a++) {
			for (int q_f = 0; q_f < fQf; q_f++) {
				filename = "Fq_Aq_" + String.format("%03d", q_f) + "_"
						+ String.format("%03d", q_a) + ".bin";
				Fq_Aq_representor_norms[q_f][q_a] = mr.readRawDoubleVector(m.getInStream(filename));
			}
		}
		Log.d(DEBUG_TAG, "Finished reading Fq_Aq_norms.dat");

		/*
		 * Read in Aq_Aq representor norm data
		 */
		int triuQa = getQa() * (getQa() + 1) / 2;
		Aq_Aq_representor_norms = new double[triuQa][][];
		for (int i = 0; i < getQa(); i++) {
			for (int j = 0; j < getQa() - i; j++) {
				filename = "Aq_Aq_" + String.format("%03d", i) + "_"
						+ String.format("%03d", j) + "_norms.bin";
				Aq_Aq_representor_norms[j + getQa() * i] = mr.readRawDoubleMatrix(m.getInStream(filename));
			}
		}
		Log.d(DEBUG_TAG, "Finished reading Aq_Aq_norms.dat");

		/*
		 * Reading uL data, if some is present
		 */
		if (fQuL > 0) {
			uL_vector = new float[fQuL][];
			for (int q_uL = 0; q_uL < fQuL; q_uL++) {
				filename = "uL_" + String.format("%03d", q_uL) + ".bin";
				uL_vector[q_uL] = mr.readRawFloatVector(m.getInStream(filename));
			}
			Log.d(DEBUG_TAG, "Finished reading uL.dat");
		}

		/*
		 *  Read in Z data
		 */
		if (fNumFields > 0) {
			Z_vector = new float[fNumFields][n_bfs][];
			for (int imf = 0; imf < fNumFields; imf++)
				for (int inbfs = 0; inbfs < n_bfs; inbfs++) {
					filename = "Z_" + String.format("%03d", imf) + "_"
							+ String.format("%03d", inbfs) + ".bin";
					Z_vector[imf][inbfs] = mr.readRawFloatVector(m.getInStream(filename));
				}
			Log.d(DEBUG_TAG, "Finished reading Z data");
		}
	}

	/**
	 * 
	 * @param m
	 * @throws IOException
	 */
	protected void loadOfflineData_rbappmit(AModelManager m) throws IOException {

		{
			BufferedReader reader = m.getBufReader("n_bfs.dat");

			String line = reader.readLine();
			n_bfs = Integer.parseInt(line);
			reader.close();
			reader = null;

			Log.d(DEBUG_TAG, "Finished reading n_bfs.dat");
		}

		{
			if (getNumFields() > 0) {
				BufferedReader reader = m.getBufReader("calN.dat");

				String line = reader.readLine();

				calN = Integer.parseInt(line);
				reader.close();
				reader = null;
			}
			Log.d(DEBUG_TAG, "Finished reading calN.dat");
		}

		// Read in output data
		if (getNumOutputs() > 0) {
			// Get output dual norms
			{
				RB_output_vectors = new RealVector[getNumOutputs()][];
				output_dual_norms = new double[getNumOutputs()][];
				String[] dual_norms_tokens;
				for (int i = 0; i < getNumOutputs(); i++) {
					{
						BufferedReader reader = m.getBufReader("output_"
								+ String.format("%03d", i) + "_dual_norms.dat");

						String line1 = reader.readLine();
						reader.close();
						reader = null;
						dual_norms_tokens = line1.split(" ");
					}

					{
						int Q_l_hat = getQl(i) * (getQl(i) + 1) / 2;
						output_dual_norms[i] = new double[Q_l_hat];
						for (int q = 0; q < Q_l_hat; q++) {
							output_dual_norms[i][q] = Double.parseDouble(dual_norms_tokens[q]);
						}
					}

					RB_output_vectors[i] = new RealVector[getQl(i)];
					String[] output_i_tokens;
					for (int q_l = 0; q_l < getQl(i); q_l++) {
						// Now read in the RB output vectors
						{
							BufferedReader reader_i = m.getBufReader("output_"
									+ String.format("%03d", i) + "_"
									+ String.format("%03d", q_l) + ".dat");

							String line_i = reader_i.readLine();
							reader_i.close();
							output_i_tokens = line_i.split(" ");
						}

						RB_output_vectors[i][q_l] = new ArrayRealVector(n_bfs);
						for (int j = 0; j < n_bfs; j++) {
							RB_output_vectors[i][q_l].setEntry(j, Double.parseDouble(output_i_tokens[j]));
						}

					}
				}
			}
		}

		Log.d(DEBUG_TAG, "Finished reading output data");

		// Read in the F_q vectors
		{
			RB_F_q_vector = new RealVector[getQf()];
			String[] tokens;
			for (int q_f = 0; q_f < getQf(); q_f++) {
				{
					BufferedReader reader = m.getBufReader("RB_F_"
							+ String.format("%03d", q_f) + ".dat");
					String line = reader.readLine();
					reader.close();
					reader = null;
					tokens = line.split(" ");
				}

				// Set the size of the inner product matrix
				RB_F_q_vector[q_f] = new ArrayRealVector(n_bfs);

				// Fill the vector
				for (int i = 0; i < n_bfs; i++) {
					RB_F_q_vector[q_f].setEntry(i, Double.parseDouble(tokens[i]));
				}

			}
			Log.d(DEBUG_TAG, "Finished reading RB_F_q data");
		}

		// Read in the A_q matrices
		{
			RB_A_q_vector = new RealMatrix[getQa()];
			String[] tokens;
			for (int q_a = 0; q_a < getQa(); q_a++) {
				{
					BufferedReader reader = m.getBufReader("RB_A_"
							+ String.format("%03d", q_a) + ".dat");
					String line = reader.readLine();
					reader.close();
					reader = null;
					tokens = line.split(" ");
				}

				// Set the size of the inner product matrix
				RB_A_q_vector[q_a] = new Array2DRowRealMatrix(n_bfs, n_bfs);

				// Fill the vector
				int count = 0;
				for (int i = 0; i < n_bfs; i++)
					for (int j = 0; j < n_bfs; j++) {
						RB_A_q_vector[q_a].setEntry(i, j, Double.parseDouble(tokens[count]));
						count++;
					}

			}
			Log.d(DEBUG_TAG, "Finished reading RB_A_q data");
		}

		// Read in F_q representor norm data
		{
			BufferedReader reader = m.getBufReader("Fq_norms.dat");

			String line = reader.readLine();
			reader.close();
			reader = null;
			String[] tokens = line.split(" ");

			// Declare the array
			int Q_f_hat = getQf() * (getQf() + 1) / 2;
			Fq_representor_norms = new double[Q_f_hat];

			// Fill it
			for (int i = 0; i < Q_f_hat; i++) {
				Fq_representor_norms[i] = Double.parseDouble(tokens[i]);
			}

			Log.d(DEBUG_TAG, "Finished reading Fq_norms.dat");
		}

		// Read in Fq_Aq representor norm data
		{
			BufferedReader reader = m.getBufReader("Fq_Aq_norms.dat");
			String line = reader.readLine();
			reader.close();
			reader = null;
			String[] tokens = line.split(" ");

			// Declare the array
			Fq_Aq_representor_norms = new double[getQf()][getQa()][n_bfs];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < getQf(); q_f++)
				for (int q_a = 0; q_a < getQa(); q_a++)
					for (int i = 0; i < n_bfs; i++) {
						Fq_Aq_representor_norms[q_f][q_a][i] = Double.parseDouble(tokens[count]);
						count++;
					}

			Log.d(DEBUG_TAG, "Finished reading Fq_Aq_norms.dat");
		}

		MathObjectReader mr = m.getMathObjReader();

		// Read in Aq_Aq representor norm data
		{
			// Declare the array
			int Q_a_hat = getQa() * (getQa() + 1) / 2;
			Aq_Aq_representor_norms = new double[Q_a_hat][][];

			int count = 0;
			String file;
			for (int i = 0; i < getQa(); i++) {
				for (int j = i; j < getQa(); j++) {
					file = "Aq_Aq_" + String.format("%03d", i) + "_"
							+ String.format("%03d", j) + "_norms.bin";
					InputStream in = m.getInStream(file);
					try {
						Aq_Aq_representor_norms[count] = mr.readRawDoubleMatrix(in, n_bfs, n_bfs);
					} finally {
						in.close();
					}
					count++;
				}
			}
			Log.d(DEBUG_TAG, "Finished reading Aq_Aq_norms.dat");
		}

		// Reading uL data
		{
			if (getQuL() > 0) {
				uL_vector = new float[getQuL()][];
				InputStream in;
				for (int q_uL = 0; q_uL < getQuL(); q_uL++) {
					in = m.getInStream("uL_" + String.format("%03d", q_uL)
							+ ".bin");
					try {
						uL_vector[q_uL] = mr.readRawFloatVector(in, get_calN());
					} finally {
						in.close();
					}
				}
			}
			Log.d(DEBUG_TAG, "Finished reading uL.dat");
		}

		// Read in Z data
		{
			int mf = getNumFields();
			if (mf > 0) {
				Z_vector = new float[mf][n_bfs][];
				InputStream in;
				for (int imf = 0; imf < mf; imf++)
					for (int inbfs = 0; inbfs < n_bfs; inbfs++) {
						in = m.getInStream("Z_" + String.format("%03d", imf)
								+ "_" + String.format("%03d", inbfs) + ".bin");
						try {
							Z_vector[imf][inbfs] = mr.readRawFloatVector(in, calN);
						} finally {
							in.close();
						}
					}
			}

			Log.d(DEBUG_TAG, "Finished reading Z.dat");
		}

	}

	/**
	 * Perform online solve with the N RB basis functions, for the set of
	 * parameters in current_params, where 1 <= N <= RB_size.
	 */
	public double RB_solve(int N) {

		current_N = N;

		if (N > getNBF()) {
			throw new RuntimeException("ERROR: N cannot be larger than the number "
					+ "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException("ERROR: N must be greater than 0 in RB_solve");
		}

		// Assemble the RB system
		RealMatrix RB_system_matrix_N = new Array2DRowRealMatrix(N, N);

		for (int q_a = 0; q_a < getQa(); q_a++) {
			RB_system_matrix_N = RB_system_matrix_N.add(RB_A_q_vector[q_a].getSubMatrix(0, N - 1, 0, N - 1).scalarMultiply(eval_theta_q_a(q_a)));
		}

		// Assemble the RB rhs
		RealVector RB_rhs_N = new ArrayRealVector(N);

		for (int q_f = 0; q_f < fQf; q_f++) {
			// Note getSubVector takes an initial index and the number of
			// entries
			// i.e. the interface is a bit different to getSubMatrix
			RB_rhs_N = RB_rhs_N.add(RB_F_q_vector[q_f].getSubVector(0, N).mapMultiply(eval_theta_q_f(q_f)));
		}

		// Solve the linear system
		DecompositionSolver solver = new LUDecompositionImpl(RB_system_matrix_N).getSolver();
		RB_solution = solver.solve(RB_rhs_N);

		// Evaluate the dual norm of the residual for RB_solution_vector
		double epsilon_N = compute_residual_dual_norm(N);

		// Get lower bound for coercivity constant
		double alpha_LB = get_SCM_lower_bound();

		// If SCM lower bound is negative
		if (alpha_LB <= 0) { // Get an upper bound instead
			alpha_LB = get_SCM_upper_bound();
		}

		// Store (absolute) error bound
		double abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

		// Compute the norm of RB_solution
		double RB_solution_norm = RB_solution.getNorm();

		// Now compute the outputs and associated errors
		RealVector RB_output_vector_N = new ArrayRealVector(N);
		for (int i = 0; i < getNumOutputs(); i++) {
			RB_outputs[i] = 0.;

			for (int q_l = 0; q_l < getQl(i); q_l++) {
				RB_output_vector_N = RB_output_vectors[i][q_l].getSubVector(0, N);
				RB_outputs[i] += eval_theta_q_l(i, q_l)
						* RB_solution.dotProduct(RB_output_vector_N);
			}
			RB_output_error_bounds[i] = compute_output_dual_norm(i, 0) // Zero is current time, not in RBSystem
					* abs_error_bound;
		}

		return (return_rel_error_bound ? abs_error_bound / RB_solution_norm : abs_error_bound);
	}

	@Override
	protected void readConfigurationJRB(AModelManager m) {
		super.readConfigurationJRB(m);
		// Number of basis functions
		n_bfs = Integer.parseInt(m.getModelXMLAttribute("num_basisfcn", "rb_model"));
		
		fNumFields = Integer.parseInt(m.getModelXMLAttribute("fields", "rb_model"));

		calN = Integer.parseInt(m.getModelXMLTagValue("geometry.nodes"));
		
		// Number of output functionals
		fNumOutputs = affineFunctionsInstance.getNumOutputs();

		Ql_values = affineFunctionsInstance.getQl();
		fQf = affineFunctionsInstance.getQf();
		if (affineFunctionsInstance instanceof IWithuL) {
			fQuL = ((IWithuL) affineFunctionsInstance).getQuL();
		}
	}

	@Override
	protected void readConfigurationRBAppMIT(GetPot infile) {
		super.readConfigurationRBAppMIT(infile);
		/*
		 * Read config values from the input.in file
		 */
		Log.d(DEBUG_TAG, "Entered parse_parameters_file, filename = "
				+ Const.parameters_filename);

		fNumFields = infile.call("n_field", 1);
		Log.d(DEBUG_TAG, "n_field = " + fNumFields);

		int return_rel_error_bound_in = infile.call("return_rel_error_bound", 1);
		return_rel_error_bound = (return_rel_error_bound_in != 0);

		Log.d(DEBUG_TAG, "return a relative error bound from RB_solve? "
				+ return_rel_error_bound);

		/*
		 * Read values using the new IAffineFunctions wrapper
		 */
		Ql_values = affineFunctionsInstance.getQl();
		fQf = affineFunctionsInstance.getQf();
		Log.d(DEBUG_TAG, "Q_f = " + fQf);
		if (affineFunctionsInstance instanceof IWithuL) {
			fQuL = ((IWithuL) affineFunctionsInstance).getQuL();
		}
		
		fNumOutputs = affineFunctionsInstance.getNumOutputs();
		Log.d(DEBUG_TAG, "n_outputs = " + fNumOutputs);
	}

	/**
	 * Specifies the residual scaling on the denominator to be used in the a
	 * posteriori error bound. Overload in subclass in order to obtain the
	 * desired error bound.
	 */
	protected double residual_scaling_denom(double alpha_LB) {
		return Math.sqrt(alpha_LB);
	}

	// Set calN
	public void set_calN(int _calN) {
		calN = _calN;
	}

	public void set_n_basis_functions(int _N) {
		n_bfs = _N;
	}

	public void set_sweep_sol(double[][][] _sweep_sol) {
		RB_sweep_solution = _sweep_sol;
	}

	/**
	 * Set the primary SCM
	 */
	public void setPrimarySCM(RBSCMSystem scm_system) {
		mRbScmSystem = scm_system;
	}

}
