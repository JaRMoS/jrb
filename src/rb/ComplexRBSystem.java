package rb;

// rbAPPmit: An Android front-end for the Certified Reduced Basis Method
// Copyright (C) 2010 David J. Knezevic and Phuong Huynh
//
// This file is part of rbAPPmit
//
// @ref rbappmit is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// @ref rbappmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with rbAPPmit. If not, see <http://www.gnu.org/licenses/>.

import jarmos.ComplexSolutionField;
import jarmos.FieldDescriptor;
import jarmos.Log;
import jarmos.SimulationResult;
import jarmos.geometry.MeshTransform;
import jarmos.io.AModelManager;
import jarmos.io.MathObjectReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Array;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.linear.Array2DRowFieldMatrix;
import org.apache.commons.math.linear.ArrayFieldVector;
import org.apache.commons.math.linear.FieldMatrix;
import org.apache.commons.math.linear.FieldVector;

/**
 * RB system class for complex-valued fields.
 * 
 * @author Daniel Wirtz @date 07.08.2013
 * 
 */
public class ComplexRBSystem extends RBSystem {

	// Logging tag
	private static final String DEBUG_TAG = "ComplexRBSystem";

	public Complex[] RB_outputs;
	public Complex[] RB_output_error_bounds;

	protected FieldVector<Complex> RB_solution;

	public Complex[][] output_dual_norms;
	protected FieldVector<Complex>[][] RB_output_vectors;
	protected FieldMatrix<Complex>[] RB_A_q_vector;
	protected FieldVector<Complex>[] RB_F_q_vector;

	protected FieldVector<Complex> theta_a;

	protected Complex[] Fq_representor_norms;
	protected Complex[][][] Fq_Aq_representor_norms;
	protected Complex[][][] Aq_Aq_representor_norms;

	protected Complex[][][] Z_vector;
	protected Complex[][] uL_vector;

	@Override
	public void loadOfflineData_rbappmit(AModelManager m) throws IOException {

		isReal = false;

		// Find out dimension of the RB space
		{
			BufferedReader reader = m.getBufReader("n_bfs.dat");

			set_n_basis_functions(Integer.parseInt(reader.readLine()));

			reader.close();
		}

		Log.d(DEBUG_TAG, "Finished reading n_bfs.dat");

		// Read in output data
		if (getNumOutputs() > 0) {
			// Get output dual norms
			{
				RB_output_vectors = (FieldVector<Complex>[][]) new FieldVector<?>[getNumOutputs()][];
				output_dual_norms = new Complex[getNumOutputs()][];
				for (int i = 0; i < getNumOutputs(); i++) {

					String[] dual_norms_tokens;
					{
						BufferedReader reader = m
								.getBufReader("output_" + String.format("%03d", i) + "_dual_norms.dat");

						dual_norms_tokens = reader.readLine().split(" ");
						reader.close();
					}

					{
						int Q_l_hat = getQl(i) * (getQl(i) + 1) / 2;
						output_dual_norms[i] = new Complex[Q_l_hat];
						for (int q = 0; q < Q_l_hat; q++) {
							output_dual_norms[i][q] = new Complex(Double.parseDouble(dual_norms_tokens[q]),
									Double.parseDouble(dual_norms_tokens[Q_l_hat + q]));
						}
					}

					{
						RB_output_vectors[i] = (FieldVector<Complex>[]) new FieldVector<?>[getQl(i)];
						String[] output_i_tokens;
						for (int q_l = 0; q_l < getQl(i); q_l++) {
							// Now read in the RB output vectors
							{
								BufferedReader reader_i = m.getBufReader("output_" + String.format("%03d", i) + "_"
										+ String.format("%03d", q_l) + ".dat");

								output_i_tokens = reader_i.readLine().split(" ");
								reader_i.close();
							}

							RB_output_vectors[i][q_l] = new ArrayFieldVector<Complex>(getNBF(), new Complex(0d, 0d));
							for (int j = 0; j < getNBF(); j++) {
								RB_output_vectors[i][q_l].setEntry(
										j,
										new Complex(Double.parseDouble(output_i_tokens[j]), Double
												.parseDouble(output_i_tokens[getNBF() + j])));
							}
						}
					}
				}
			}
			Log.d(DEBUG_TAG, "Finished reading output data");
		} else
			Log.d(DEBUG_TAG, "No output data set. (get_n_outputs() == 0)");

		/*
		 * // Read in the inner product matrix { InputStreamReader isr; String
		 * dataString = directory_name + "/RB_inner_product_matrix.dat";
		 * 
		 * if(!isAssetFile) { HttpGet request = new HttpGet(dataString);
		 * HttpResponse response = client.execute(request); isr = new
		 * InputStreamReader(response.getEntity() .getContent()); } else { //
		 * Read from assets isr = new InputStreamReader(
		 * context.getAssets().open(dataString)); } BufferedReader
		 * BufferedReader reader = new BufferedReader(isr,buffer_size);
		 * 
		 * String String line = reader.readLine(); String[] tokens =
		 * line.split(" ");
		 * 
		 * // Set the size of the inner product matrix RB_inner_product_matrix =
		 * new Array2DRowRealMatrix(n_bfs, n_bfs);
		 * 
		 * // Fill the matrix int count = 0; for (int i = 0; i < n_bfs; i++) for
		 * (int j = 0; j < n_bfs; j++) { RB_inner_product_matrix.setEntry(i, j,
		 * Double .parseDouble(tokens[count])); count++; } reader.close(); }
		 * 
		 * Log.d(DEBUG_TAG, "Finished reading RB_inner_product_matrix.dat");
		 */

		// Read in the F_q vectors
		{
			RB_F_q_vector = (FieldVector<Complex>[]) new FieldVector<?>[getQf()];
			String[] tokens;
			for (int q_f = 0; q_f < getQf(); q_f++) {
				{
					BufferedReader reader = m.getBufReader("RB_F_" + String.format("%03d", q_f) + ".dat");

					tokens = reader.readLine().split(" ");
					reader.close();
				}

				// Set the size of the inner product matrix
				RB_F_q_vector[q_f] = new ArrayFieldVector<Complex>(getNBF(), new Complex(0d, 0d));

				// Fill the vector
				for (int i = 0; i < getNBF(); i++) {
					RB_F_q_vector[q_f].setEntry(i,
							new Complex(Double.parseDouble(tokens[i]), Double.parseDouble(tokens[getNBF() + i])));
				}
			}
			Log.d(DEBUG_TAG, "Finished reading RB_F_q data");
		}

		// Read in the A_q matrices
		{
			RB_A_q_vector = (FieldMatrix<Complex>[]) Array.newInstance(FieldMatrix.class, getQa());// (FieldMatrix<Complex>[])
																									// new
																									// FieldMatrix<?>[get_Q_a()];
			String[] tokens;
			for (int q_a = 0; q_a < getQa(); q_a++) {
				{
					BufferedReader reader = m.getBufReader("RB_A_" + String.format("%03d", q_a) + ".dat");
					tokens = reader.readLine().split(" ");
					reader.close();
				}

				// Set the size of the inner product matrix
				RB_A_q_vector[q_a] = new Array2DRowFieldMatrix<Complex>((new Complex(0, 0)).getField(), getNBF(),
						getNBF());

				// Fill the vector
				int count = 0;
				for (int i = 0; i < getNBF(); i++)
					for (int j = 0; j < getNBF(); j++) {
						RB_A_q_vector[q_a].setEntry(
								i,
								j,
								new Complex(Double.parseDouble(tokens[count]), Double.parseDouble(tokens[count
										+ getNBF() * getNBF()])));
						count++;
					}
			}
			Log.d(DEBUG_TAG, "Finished reading RB_A_q data");
		}

		// Read in F_q representor norm data
		{
			BufferedReader reader = m.getBufReader("Fq_norms.dat");

			String[] tokens = reader.readLine().split(" ");
			reader.close();

			// Declare the array
			int Q_f_hat = getQf() * (getQf() + 1) / 2;
			Fq_representor_norms = new Complex[Q_f_hat];

			// Fill it
			for (int i = 0; i < Q_f_hat; i++) {
				Fq_representor_norms[i] = new Complex(Double.parseDouble(tokens[i * 2 + 0]),
						Double.parseDouble(tokens[i * 2 + 1]));
			}

			Log.d(DEBUG_TAG, "Finished reading Fq_norms.dat");
		}

		// Read in Fq_Aq representor norm data
		{
			BufferedReader reader = m.getBufReader("Fq_Aq_norms.dat");

			String[] tokens = reader.readLine().split(" ");
			reader.close();
			reader = null;

			// Declare the array
			Fq_Aq_representor_norms = new Complex[getQf()][getQa()][getNBF()];

			double[][][] Rdata = new double[getQf()][getQa()][getNBF()];
			double[][][] Idata = new double[getQf()][getQa()][getNBF()];
			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < getQf(); q_f++)
				for (int q_a = 0; q_a < getQa(); q_a++) {
					for (int i = 0; i < getNBF(); i++) {
						Rdata[q_f][q_a][i] = Double.parseDouble(tokens[count]);
						count++;
					}
					for (int i = 0; i < getNBF(); i++) {
						Idata[q_f][q_a][i] = Double.parseDouble(tokens[count]);
						count++;
					}
				}

			for (int q_f = 0; q_f < getQf(); q_f++)
				for (int q_a = 0; q_a < getQa(); q_a++)
					for (int i = 0; i < getNBF(); i++)
						Fq_Aq_representor_norms[q_f][q_a][i] = new Complex(Rdata[q_f][q_a][i], Idata[q_f][q_a][i]);

			Log.d(DEBUG_TAG, "Finished reading Fq_Aq_norms.dat");
		}

		MathObjectReader mr = m.getMathObjReader();
		// Read in Aq_Aq representor norm data
		{
			// Declare the array
			int Q_a_hat = getQa() * (getQa() + 1) / 2;
			Aq_Aq_representor_norms = new Complex[Q_a_hat][getNBF()][getNBF()];

			int count = 0;
			double[][] Rdata2 = null, Idata2 = null;
			for (int i = 0; i < getQa(); i++)
				for (int j = i; j < getQa(); j++) {
					String file = "Aq_Aq_" + String.format("%03d", i) + "_" + String.format("%03d", j) + "_norms.bin";

					int n = getNBF();
					InputStream in = m.getInStream(file);
					try {
						Rdata2 = mr.readRawDoubleMatrix(in, n, n);
						Idata2 = mr.readRawDoubleMatrix(in, n, n);
					} catch (IOException io) {
						Log.e("ComplexRBSystem", "IOException with file " + file, io);
					} finally {
						in.close();
						in = null;
					}
					for (int k = 0; k < n; k++)
						for (int l = 0; l < n; l++)
							Aq_Aq_representor_norms[count][k][l] = new Complex(Rdata2[k][l], Idata2[k][l]);

					count++;
				}
			Log.d(DEBUG_TAG, "Finished reading Aq_Aq_norms.dat");
		}

		// // Read calN number
		// {
		// if (getNumFields() > 0) {
		// BufferedReader reader = m.getBufReader("calN.dat");
		//
		// String line = reader.readLine();
		//
		// set_calN(Integer.parseInt(line));
		// reader.close();
		// }
		//
		// Log.d(DEBUG_TAG, "Finished reading calN.dat");
		// }

		int n = getGeometry().getNumVertices();
		// Reading uL data
		{
			if (getQuL() > 0) {
				uL_vector = new Complex[getQuL()][n];
				float[] Rdata3, Idata3;
				for (int q_uL = 0; q_uL < getQuL(); q_uL++) {
					InputStream in = m.getInStream("uL_" + String.format("%03d", q_uL) + ".bin");
					try {
						Rdata3 = mr.readRawFloatVector(in, n);
						Idata3 = mr.readRawFloatVector(in, n);
					} finally {
						in.close();
					}
					for (int i = 0; i < n; i++)
						uL_vector[q_uL][i] = new Complex(Rdata3[i], Idata3[i]);
				}
			}
			Log.d(DEBUG_TAG, "Finished reading uL.dat");
		}

		// Read in Z data
		{
			if (getNumDoFFields() > 0) {
				Z_vector = new Complex[getNumDoFFields()][getNBF()][n];
				float[] Rdata3, Idata3;
				for (int imf = 0; imf < getNumDoFFields(); imf++)
					for (int inbfs = 0; inbfs < getNBF(); inbfs++) {
						InputStream in = m.getInStream("Z_" + String.format("%03d", imf) + "_"
								+ String.format("%03d", inbfs) + ".bin");
						try {
							Rdata3 = mr.readRawFloatVector(in, n);
							Idata3 = mr.readRawFloatVector(in, n);
						} finally {
							in.close();
						}
						for (int i = 0; i < n; i++)
							Z_vector[imf][inbfs][i] = new Complex(Rdata3[i], Idata3[i]);
					}
			}
			Log.d(DEBUG_TAG, "Finished reading Z.dat");
		}

		initialize_data_vectors();
	}

	protected void initialize_data_vectors() {
		// Also, resize RB_outputs and RB_output_error_error_bounds arrays
		RB_outputs = new Complex[getNumOutputs()];
		RB_output_error_bounds = new Complex[getNumOutputs()];
	}

	@Override
	public double RB_solve(int N) {

		current_N = N;

		theta_a = complex_eval_theta_q_a();

		if (N > getNBF()) {
			throw new RuntimeException("ERROR: N cannot be larger than the number " + "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException("ERROR: N must be greater than 0 in RB_solve");
		}

		// Assemble the RB system
		FieldMatrix<Complex> RB_system_matrix_N = new Array2DRowFieldMatrix<Complex>((new Complex(0, 0)).getField(), N,
				N);

		for (int q_a = 0; q_a < getQa(); q_a++) {
			RB_system_matrix_N = RB_system_matrix_N.add(RB_A_q_vector[q_a].getSubMatrix(0, N - 1, 0, N - 1)
					.scalarMultiply(theta_a.getEntry(q_a)));
			// scalarMultiply(complex_eval_theta_q_a(q_a) ) );
		}

		// Assemble the RB rhs
		FieldVector<Complex> RB_rhs_N = new ArrayFieldVector<Complex>(N, new Complex(0d, 0d));

		for (int q_f = 0; q_f < getQf(); q_f++) {
			// Note getSubVector takes an initial index and the number of
			// entries
			// i.e. the interface is a bit different to getSubMatrix
			RB_rhs_N = RB_rhs_N.add(RB_F_q_vector[q_f].getSubVector(0, N).mapMultiply(complex_eval_theta_q_f(q_f)));
		}

		// Solve the linear system by Gaussian elimination

		RB_solution = new ArrayFieldVector<Complex>(N, new Complex(0., 0.));
		for (int j = 1; j < N; j++)
			for (int i = j; i < N; i++) {
				Complex m = RB_system_matrix_N.getEntry(i, j - 1).divide(RB_system_matrix_N.getEntry(j - 1, j - 1));
				for (int k = 0; k < N; k++)
					RB_system_matrix_N.setEntry(
							i,
							k,
							RB_system_matrix_N.getEntry(i, k).subtract(
									RB_system_matrix_N.getEntry(j - 1, k).multiply(m)));
				RB_rhs_N.setEntry(i, RB_rhs_N.getEntry(i).subtract(m.multiply(RB_rhs_N.getEntry(j - 1))));
			}
		RB_solution.setEntry(N - 1, RB_rhs_N.getEntry(N - 1).divide(RB_system_matrix_N.getEntry(N - 1, N - 1)));
		for (int j = N - 2; j >= 0; j--) {
			Complex m = new Complex(0., 0.);
			for (int i = j + 1; i < N; i++)
				m = m.add(RB_system_matrix_N.getEntry(j, i).multiply(RB_solution.getEntry(i)));
			RB_solution.setEntry(j, (RB_rhs_N.getEntry(j).subtract(m)).divide(RB_system_matrix_N.getEntry(j, j)));
		}

		// Evaluate the dual norm of the residual for RB_solution_vector
		double epsilon_N = compute_residual_dual_norm(N);

		// Get lower bound for coercivity constant
		double alpha_LB = get_SCM_lower_bound();

		// If SCM lower bound is negative
		if (alpha_LB < 0) { // Get an upper bound instead
			alpha_LB = get_SCM_upper_bound();
		}

		// Store (absolute) error bound
		double abs_error_bound = epsilon_N / residual_scaling_denom(alpha_LB);

		// Compute the norm of RB_solution
		/*
		 * RealMatrix RB_inner_product_matrix_N =
		 * RB_inner_product_matrix.getSubMatrix(0, N-1, 0, N-1);
		 */

		double RB_solution_norm = 0.0d;
		for (int i = 0; i < N; i++)
			RB_solution_norm += ((RB_solution.getEntry(i)).multiply((RB_solution.getEntry(i)).conjugate())).getReal();
		RB_solution_norm = Math.sqrt(RB_solution_norm);

		// Now compute the outputs and associated errors
		FieldVector<Complex> RB_output_vector_N = new ArrayFieldVector<Complex>(N, new Complex(0d, 0d));
		for (int i = 0; i < getNumOutputs(); i++) {
			RB_outputs[i] = new Complex(0., 0.);

			RB_output_vector_N = (RB_output_vectors[i][0].getSubVector(0, N)).mapMultiply(complex_eval_theta_q_l(i, 0));
			for (int q_l = 1; q_l < getQl(i); q_l++)
				RB_output_vector_N = RB_output_vector_N.add((RB_output_vectors[i][q_l].getSubVector(0, N))
						.mapMultiply(complex_eval_theta_q_l(i, q_l)));
			for (int j = 0; j < N; j++)
				RB_outputs[i] = RB_outputs[i].add((RB_output_vector_N.getEntry(j).conjugate()).multiply((RB_solution
						.getEntry(j))));

			RB_output_error_bounds[i] = new Complex(compute_output_dual_norm(i, 0) // Zero
																					// means
																					// no
																					// time
																					// used
																					// here
					* abs_error_bound, compute_output_dual_norm(i, 0) // Zero
																		// means
																		// no
																		// time
																		// used
																		// here
					* abs_error_bound);
		}

		cal_derived_output();

		return (return_rel_error_bound ? abs_error_bound / RB_solution_norm : abs_error_bound);
	}

	@Override
	protected double compute_residual_dual_norm(int N) {

		// Use the stored representor inner product values
		// to evaluate the residual norm
		double res_ff = 0;
		double res_af = 0;
		double res_aa = 0;

		int q = 0;
		for (int q_f1 = 0; q_f1 < getQf(); q_f1++) {
			for (int q_f2 = q_f1; q_f2 < getQf(); q_f2++) {
				double delta = (q_f1 == q_f2) ? 1. : 2.;
				res_ff += delta
						* ((complex_eval_theta_q_f(q_f1).multiply(complex_eval_theta_q_f(q_f2).conjugate()))
								.multiply(Fq_representor_norms[q])).getReal();
				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			for (int q_a = 0; q_a < getQa(); q_a++) {
				for (int i = 0; i < N; i++) {
					res_af += 2. * (get_complex_soln_coeff(i).conjugate().multiply(
					// complex_eval_theta_q_f(q_f).multiply(complex_eval_theta_q_a(q_a).conjugate())
							complex_eval_theta_q_f(q_f).multiply(theta_a.getEntry(q_a).conjugate()))
							.multiply(Fq_Aq_representor_norms[q_f][q_a][i])).getReal();
				}
			}
		}

		q = 0;
		for (int q_a1 = 0; q_a1 < getQa(); q_a1++) {
			for (int q_a2 = q_a1; q_a2 < getQa(); q_a2++) {
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < N; j++) {
						double delta = (q_a1 == q_a2) ? 1. : 2.;
						res_aa += delta
								* ((get_complex_soln_coeff(i).conjugate().multiply(get_complex_soln_coeff(j)))
										.multiply(
										// (complex_eval_theta_q_a(q_a1).conjugate().multiply(complex_eval_theta_q_a(q_a2)))
										(theta_a.getEntry(q_a1).conjugate().multiply(theta_a.getEntry(q_a2))))
										.multiply(Aq_Aq_representor_norms[q][i][j])).getReal();
					}
				}
				q++;
			}
		}

		double residual_norm_sq = res_ff + res_af + res_aa;

		if (residual_norm_sq < 0.) {
			// Sometimes this is negative due to rounding error,
			// but error is on the order of 1.e-10, so shouldn't
			// affect error bound much...
			residual_norm_sq = Math.abs(residual_norm_sq);
		}

		return Math.sqrt(residual_norm_sq);
	}

	/**
	 * @see rb.RBSystem#compute_output_dual_norm(int, double) TODO: make use of the t variable and change the
	 * complex_eval_theta etc
	 */
	@Override
	protected double compute_output_dual_norm(int i, double t) {

		// Use the stored representor inner product values
		// to evaluate the output dual norm
		double output_norm_sq = 0.;

		int q = 0;
		for (int q_l1 = 0; q_l1 < getQl(i); q_l1++) {
			for (int q_l2 = q_l1; q_l2 < getQl(i); q_l2++) {
				if (q_l1 == q_l2)
					output_norm_sq += 1. * ((complex_eval_theta_q_l(i, q_l1).multiply(complex_eval_theta_q_l(i, q_l2)
							.conjugate())).multiply(output_dual_norms[i][q])).getReal();
				else
					output_norm_sq += 2. * ((complex_eval_theta_q_l(i, q_l1).multiply(complex_eval_theta_q_l(i, q_l2)
							.conjugate())).multiply(output_dual_norms[i][q])).getReal();
				q++;
			}
		}

		return Math.sqrt(output_norm_sq);
	}

	Complex get_complex_soln_coeff(int i) {
		return RB_solution.getEntry(i);
	}

	@Override
	public double[][] get_RBsolution() {
		Complex[] c = RB_solution.toArray();
		double[][] d = new double[2][c.length];
		for (int i = 0; i < c.length; i++) {
			d[0][i] = c[i].getReal();
			d[1][i] = c[i].getImaginary();
		}
		return d;
	}

	@Override
	public double get_RB_output(int n_output, boolean Rpart) {
		if (Rpart)
			return RB_outputs[n_output].getReal();
		else
			return RB_outputs[n_output].getImaginary();
	}

	@Override
	public double get_RB_output_error_bound(int n_output, boolean Rpart) {
		if (Rpart)
			return RB_output_error_bounds[n_output].getReal();
		else
			return RB_output_error_bounds[n_output].getImaginary();
	}

	/**
	 * Returns the results of the complex RB simulation.
	 * 
	 * Formerly, the returned float[][] array was again sorted into visualization data in RBVisualization.onCreate. As
	 * there have only been @ref rbappmit models with complex field data, which has been split into real, imaginary and
	 * norm field values for visualization, this new method is the native way of transforming the results into
	 * ComplexSolutionFields, however many there may be. For each field the visualization will automatically create
	 * views for real, imaginary and norm value fields.
	 */
	@Override
	public SimulationResult getSimulationResults() {
		int N = RB_solution.getDimension();
		SimulationResult res = new SimulationResult(1);

		int fnumcnt = 0;
		for (FieldDescriptor sftype : logicalFieldTypes) {
			if (fnumcnt + sftype.Type.requiredDoFFields > getNumDoFFields()) {
				throw new RuntimeException("Too many output fields used by current "
						+ "SolutionFieldTypes set in RBSystem. Check your model.xml.");
			}
			int numDoF = Z_vector[fnumcnt][0].length;
			switch (sftype.Type) {
			case ComplexValue: {
				ComplexSolutionField f = new ComplexSolutionField(sftype, numDoF);
				Complex tmpval;
				for (int nodenr = 0; nodenr < numDoF; nodenr++) {
					tmpval = new Complex(0., 0.);
					for (int rbdim = 0; rbdim < N; rbdim++) {
						tmpval = tmpval.add(Z_vector[fnumcnt][rbdim][nodenr].multiply(get_complex_soln_coeff(rbdim)));
					}
					f.setComplexValue(nodenr, tmpval);
				}
				if (getQuL() > 0) {
					for (int q_uL = 0; q_uL < getQuL(); q_uL++)
						for (int i = 0; i < numDoF; i++) {
							f.addComplexValue(i, (float) uL_vector[q_uL][i].getReal(),
									(float) uL_vector[q_uL][i].getImaginary());
						}
				}
				res.addField(f);
				break;
			}
			default:
				throw new RuntimeException("Invalid/unimplemented solution field type '" + sftype.Type
						+ "' for complex RB system");
			}
			/*
			 * Increase field counter by the amount the current solution field
			 * used
			 */
			fnumcnt += sftype.Type.requiredDoFFields;
		}
		for (MeshTransform m : transforms) {
			res.addTransform(m);
		}
		return res;
	}

	@Override
	public SimulationResult getSweepSimResults() {
		int N = RB_sweep_solution[0][0].length;
		int numSweep = RB_sweep_solution.length;

		/*
		 * Preparations
		 */
		Complex[][] RB_sweep_sol = new Complex[numSweep][N];
		for (int i = 0; i < numSweep; i++)
			for (int j = 0; j < N; j++)
				RB_sweep_sol[i][j] = new Complex(RB_sweep_solution[i][0][j], RB_sweep_solution[i][1][j]);
		SimulationResult res = new SimulationResult(numSweep);
		int fnumcnt = 0;
		for (FieldDescriptor sftype : logicalFieldTypes) {
			if (fnumcnt + sftype.Type.requiredDoFFields > getNumDoFFields()) {
				throw new RuntimeException("Too many output fields used by current "
						+ "SolutionFieldTypes set in RBSystem. Check your model.xml.");
			}
			int numDoF = Z_vector[fnumcnt][0].length;
			switch (sftype.Type) {
			case ComplexValue: {
				ComplexSolutionField f = new ComplexSolutionField(sftype, numDoF * numSweep);
				Complex tmpval;
				for (int iSweep = 0; iSweep < numSweep; iSweep++) {
					for (int i = 0; i < numDoF; i++) {
						tmpval = new Complex(0., 0.);
						for (int rbdim = 0; rbdim < N; rbdim++) {
							tmpval = tmpval.add(Z_vector[fnumcnt][rbdim][i].multiply(RB_sweep_sol[iSweep][rbdim]));
						}
						f.setComplexValue(iSweep * numDoF + i, tmpval);
					}
				}
				if (getQuL() > 0) {
					for (int q_uL = 0; q_uL < getQuL(); q_uL++)
						for (int iSweep = 0; iSweep < numSweep; iSweep++)
							for (int nodenr = 0; nodenr < numDoF; nodenr++) {
								f.addComplexValue(iSweep * numDoF + nodenr, uL_vector[q_uL][nodenr]);
							}
				}
				res.addField(f);
				break;
			}
			default:
				throw new RuntimeException("Invalid/unimplemented solution field type '" + sftype.Type
						+ "' for complex RB system sweep");
			}
			/*
			 * Increase field counter by the amount the current solution field
			 * used
			 */
			fnumcnt += sftype.Type.requiredDoFFields;
		}
		for (MeshTransform m : transforms) {
			res.addTransform(m);
		}
		return res;
	}

	public boolean is_derived_output() {
		Method meth;

		try {
			Class<?> partypes[] = null;
			meth = oldAffFcnCl.getMethod("is_derived_output", partypes);
		} catch (NoSuchMethodException nsme) {
			return false;
		}

		try {
			Object arglist[] = null;
			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			boolean val = (Boolean) theta_obj;
			return val;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
	}

	public void cal_derived_output() {
		if (is_derived_output()) {
			Method meth;

			try {

				Class<?> partypes[] = new Class[1];
				partypes[0] = double[].class;

				meth = oldAffFcnCl.getMethod("cal_derived_output", partypes);
			} catch (NoSuchMethodException nsme) {
				throw new RuntimeException("getMethod for cal_derived_output failed", nsme);
			}

			try {

				Object arglist[] = new Object[1];
				double[] input = new double[4];
				for (int i = 0; i < getNumOutputs(); i++) {
					input[0] = RB_outputs[i].getReal();
					input[1] = RB_outputs[i].getImaginary();
					input[2] = RB_output_error_bounds[i].getReal();
					input[3] = RB_output_error_bounds[i].getImaginary();

					arglist[0] = input;

					Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
					double[] output = (double[]) theta_obj;
					RB_outputs[i] = new Complex(output[0], output[1]);
					RB_output_error_bounds[i] = new Complex(output[2], output[3]);
				}
			} catch (IllegalAccessException iae) {
				throw new RuntimeException(iae);
			} catch (InvocationTargetException ite) {
				throw new RuntimeException(ite.getCause());
			}

		}
	}
}
