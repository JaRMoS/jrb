package rb.java;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayFieldVector;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.DecompositionSolver;
import org.apache.commons.math.linear.FieldVector;
import org.apache.commons.math.linear.LUDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import rb.java.affinefcn.IAffineFunctions;
import rb.java.affinefcn.IAffineInitials;
import rb.java.affinefcn.IWithuL;
import rmcommon.FieldDescriptor;
import rmcommon.Log;
import rmcommon.ModelBase;
import rmcommon.ModelType;
import rmcommon.Parameters;
import rmcommon.DefaultSolutionField;
import rmcommon.SimulationResult;
import rmcommon.geometry.AffineLinearMeshTransform;
import rmcommon.geometry.DefaultTransform;
import rmcommon.geometry.DisplacementField;
import rmcommon.geometry.GeometryData;
import rmcommon.geometry.MeshTransform;
import rmcommon.io.AModelManager;
import rmcommon.io.MathObjectReader;
import rmcommon.io.AModelManager.ModelManagerException;
import rmcommon.io.MathObjectReader.MathReaderException;

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
 *       TODO: import & use implicit operators LL_I from rbmatlab TODO: report
 *       progress for transient RB systems and/or normal systems
 * 
 */
public class RBSystem extends ModelBase {

	// Logging tag
	private static final String DEBUG_TAG = "RBSystem";

	/**
	 * Maximum number of sweep points
	 */
	private static final int MAX_SWEEP_POINTS = 10;

	/**
	 * The Class object for the affine functions class.
	 * 
	 * Implements the interface import rb.java.affinefcn.IAffineFunctions.
	 */
	public Class<IAffineFunctions> affineFunctionsClass;
	/**
	 * The member object that defines the parameter-dependent functions for the
	 * affine expansion of the left-hand, right-hand side and outputs. We need
	 * to access this object at runtime, hence we use a ClassLoader and
	 * reflection in order to call the relevant functions.
	 */
	public IAffineFunctions affineFunctionsInstance;

	protected double[][][] Aq_Aq_representor_norms;

	// current N
	protected int current_N;

	// The number of output functionals
	private int fNumOutputs;

	protected double[][][] Fq_Aq_representor_norms;

	/**
	 * Arrays storing the residual representor inner products to be used in
	 * computing the residuals in the Online stage. These are resized by reading
	 * in the Offline data written out by rbOOmit.
	 */
	protected double[] Fq_representor_norms;

	// PUBLIC MEMBER VARIABLES

	/**
	 * The number of terms in the affine expansion of the bilinear form
	 */
	private int fQa;

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

	// PRIVATE MEMBER VARIABLES

	protected float[][][] fullBasisVectors;

	/**
	 * 
	 */
	public boolean isReal = true;

	/**
	 * A reference to the SCM system.
	 */
	public RBSCMSystem mRbScmSystem;

	/**
	 * The parameters used for a parameter sweep.
	 */
	private double[][] mSweepParam;

	/**
	 * The number of basis functions in the RB space.
	 */
	private int numBasisFuncs;

	/**
	 * @deprecated Here for rbappmit compatibility.
	 */
	public Class<?> oldAffFcnCl;

	/**
	 * @deprecated Here for rbappmit compatibility.
	 */
	public Object oldAffFcnObj;

	/**
	 * This array stores the dual norms for each output. Row n stores the Q_l
	 * dual norms for the expansion of the n^th output.
	 */
	public double[][] output_dual_norms;

	/**
	 * The system's parameters object containing the parameter values and
	 * descriptions
	 */
	private Parameters params;

	// The affine function expansion size for all outputs
	private int[] Ql_values;

	/**
	 * Dense matrices for the RB computations.
	 */
	protected RealMatrix[] RB_A_q_vector;

	/**
	 * Dense vector for the RHS.
	 */
	protected RealVector[] RB_F_q_vector;

	/**
	 * The initial value solution coefficients
	 */
	private RealVector[] RB_initial_coeffs;

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

	private double sweepIncrement = 0;

	private double[][] sweepOutputBounds;
	private double[][] sweepOutputs;

	private Method transformationMethod = null;

	protected float[][] uL_vector;

	private List<MeshTransform> transforms;

	/**
	 * Constructor.
	 */
	public RBSystem() {

		// Initialize n_bfs to 0
		numBasisFuncs = 0;
		transforms = new ArrayList<MeshTransform>();
	}

	// PUBLIC FUNCTIONS

	/**
	 * 
	 */

	public FieldVector<Complex> complex_eval_theta_q_a() {
		Method meth;
		// Complex c;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[1];
			partypes[0] = double[].class;

			meth = oldAffFcnCl.getMethod("evaluateA_array", partypes);
		} catch (NoSuchMethodException nsme) {
			FieldVector<Complex> c = new ArrayFieldVector<Complex>(fQa, new Complex(0d, 0d));
			for (int i = 0; i < fQa; i++)
				c.setEntry(i, complex_eval_theta_q_a(i));
			return c;
			// throw new RuntimeException("getMethod for evaluateA failed",
			// nsme);
		}

		double[][] theta_val;
		try {
			Object arglist[] = new Object[1];
			arglist[0] = params.getCurrent();
			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			theta_val = (double[][]) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
		FieldVector<Complex> c = new ArrayFieldVector<Complex>(fQa, new Complex(0d, 0d));
		for (int i = 0; i < fQa; i++)
			c.setEntry(i, new Complex(theta_val[i][0], theta_val[i][1]));

		return c;
	}

	public Complex complex_eval_theta_q_a(int q) {
		Method meth;
		// Complex c;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[3];
			partypes[0] = Integer.TYPE;
			partypes[1] = double[].class;
			partypes[2] = boolean.class;

			meth = oldAffFcnCl.getMethod("evaluateA", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateA failed", nsme);
		}

		Double theta_val_r, theta_val_i;
		try {
			Object arglist[] = new Object[3];
			arglist[0] = new Integer(q);
			arglist[1] = params.getCurrent();
			arglist[2] = true;

			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			theta_val_r = (Double) theta_obj;

			arglist[2] = false;
			theta_val_i = (Double) meth.invoke(oldAffFcnObj, arglist);
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
		Complex c = new Complex(theta_val_r.doubleValue(), theta_val_i.doubleValue());

		return c;
	}

	public Complex complex_eval_theta_q_f(int q) {
		Method meth;
		// Complex c;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[3];
			partypes[0] = Integer.TYPE;
			partypes[1] = double[].class;
			partypes[2] = boolean.class;

			meth = oldAffFcnCl.getMethod("evaluateF", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateF failed", nsme);
		}

		Double theta_val_r, theta_val_i;
		try {
			Object arglist[] = new Object[3];
			arglist[0] = new Integer(q);
			arglist[1] = params.getCurrent();
			arglist[2] = true;

			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			theta_val_r = (Double) theta_obj;

			arglist[2] = false;
			theta_val_i = (Double) meth.invoke(oldAffFcnObj, arglist);
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
		Complex c = new Complex(theta_val_r.doubleValue(), theta_val_i.doubleValue());

		return c;
	}

	public Complex complex_eval_theta_q_l(int n, int q) {
		Method meth;
		// Complex c;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[4];
			partypes[0] = Integer.TYPE;
			partypes[1] = Integer.TYPE;
			partypes[2] = double[].class;
			partypes[3] = boolean.class;

			meth = oldAffFcnCl.getMethod("evaluateL", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateL failed", nsme);
		}

		Double theta_val_r, theta_val_i;
		try {
			Object arglist[] = new Object[4];
			arglist[0] = new Integer(n);
			arglist[1] = new Integer(q);
			arglist[2] = params.getCurrent();
			arglist[3] = true;

			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			theta_val_r = (Double) theta_obj;

			arglist[3] = false;
			theta_val_i = (Double) meth.invoke(oldAffFcnObj, arglist);
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
		Complex c = new Complex(theta_val_r.doubleValue(), theta_val_i.doubleValue());

		return c;
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
				output_norm_sq += delta * thetaQl(i, q_l1, t) * thetaQl(i, q_l2, t) * output_dual_norms[i][q];
				q++;
			}
		}

		return Math.sqrt(output_norm_sq);
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
				residual_norm_sq += delta * thetaQf(q_f1) * thetaQf(q_f2) * Fq_representor_norms[q];

				q++;
			}
		}

		for (int q_f = 0; q_f < getQf(); q_f++) {
			for (int q_a = 0; q_a < getQa(); q_a++) {
				for (int i = 0; i < N; i++) {
					double delta = 2.;
					residual_norm_sq += get_soln_coeff(i) * delta * thetaQf(q_f) * thetaQa(q_a)
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
						residual_norm_sq += get_soln_coeff(i) * get_soln_coeff(j) * delta * thetaQa(q_a1)
								* thetaQa(q_a2) * Aq_Aq_representor_norms[q][i][j];
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

	// /**
	// * @param mu
	// * @param gData
	// */
	// public void customMeshTransform(double[][] mu, GeometryData gData) {
	// int numT = mu.length;
	//
	// gData.vertices = new float[numT * gData.getNumVertices() * 3];
	// for (int i = 0; i < numT; i++) {
	// // get current nodal data
	// float[] tmpnode = rbappmitCustomMeshTransform(mu[i],
	// gData.originalVertices.clone());
	// // Log.d("GLRenderer", mu[i][0] + " " + mu[i][1]);
	// // Log.d("GLRenderer", tmpnode[4] + " " + data.vertices[4]);
	// // copy current nodal data into animation list
	// for (int j = 0; j < gData.getNumVertices(); j++) {
	// for (int k = 0; k < 3; k++) {
	// gData.vertices[i * gData.getNumVertices() * 3 + j * 3 + k] = tmpnode[j *
	// 3 + k];
	// }
	// }
	// }
	// }

	public int get_N() {
		return current_N;
	}

	public double get_RB_output(int n_output, boolean Rpart) {
		return RB_outputs[n_output];
	}

	public double get_RB_output_error_bound(int n_output, boolean Rpart) {
		return RB_output_error_bounds[n_output];
	}

	protected double[][] get_RBsolution() {
		return new double[][] { RB_solution.toArray() };
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
			// mRbScmSystem.setParams(getParams());
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
			// mRbScmSystem.setParams(getParams());
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

	/**
	 * Returns the initial conditions solution coefficients for RB size N
	 * 
	 * TODO: Use this method in any RB_Solve member in subclasses!
	 * 
	 * @param N
	 *            The RB size
	 * @return The initial solution coefficient vector of size N
	 */
	public RealVector getInitialCoefficients(int N) {
		if (affineFunctionsInstance instanceof IAffineInitials) {
			RealVector res = new ArrayRealVector(N);
			IAffineInitials a = (IAffineInitials) affineFunctionsInstance;
			for (int i = 0; i < a.getQu0(); i++) {
				res = res.add(RB_initial_coeffs[i].getSubVector(0, N).mapMultiply(
						a.thetaQu0(i, getParams().getCurrent())));
			}
			return res;
		} else {
			return RB_initial_coeffs[0].getSubVector(0, N);
		}
	}

	/**
	 * @return The number of basis functions in the system.
	 */
	public int getNBF() {
		return numBasisFuncs;
	}

	/**
	 * @return the number of output functionals
	 */
	public int getNumOutputs() {
		return fNumOutputs;
	}

	/**
	 * Returns the system's parameters object
	 * 
	 * @return The current parameters
	 */
	public Parameters getParams() {
		return params;
	}

	/**
	 * @return Q_a, the number of term in the affine expansion of the bilinear
	 *         form
	 */
	public int getQa() {
		return fQa;
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

	/**
	 * Constructs the true solution from the most recent RB reduced solution
	 * using the full Z data and current coefficients.
	 * 
	 * @attention This is a "hack" as the former rbappmit had only a limited
	 *            number of models, of which all visualization has been
	 *            organized a certain way (previously coded inside the
	 *            RBVisualization.onCreate method). This is as follows: One
	 *            field: Only field values Two fields: If not transient, it is
	 *            2D displacement data Three fields: If not transient, it is 3D
	 *            displacement data Four fields: If not transient, it is 3D
	 *            displacement and a fourth field with normal values.
	 * 
	 * @see ComplexRBSystem for an overridden version of this method.
	 * 
	 * @return
	 */
	public SimulationResult getSimulationResults() {
		int N = RB_solution.getDimension();
		SimulationResult res = new SimulationResult(1);
		double co;
		int fnumcnt = 0;
		for (FieldDescriptor sftype : logicalFieldTypes) {
			if (fnumcnt + sftype.Type.requiredDoFFields > getNumDoFFields()) {
				throw new RuntimeException("Too many output fields used by current "
						+ "SolutionFieldTypes set in RBSystem. Check your model.xml.");
			}
			int numDoF = fullBasisVectors[fnumcnt][0].length;
			switch (sftype.Type) {
			case Displacement3D: {
				DisplacementField d = new DisplacementField(sftype, numDoF);
				for (int nodenr = 0; nodenr < numDoF; nodenr++) {
					double x = 0, y = 0, z = 0;
					for (int rbdim = 0; rbdim < N; rbdim++) {
						co = get_soln_coeff(rbdim);
						x += fullBasisVectors[fnumcnt][rbdim][nodenr] * co;
						y += fullBasisVectors[fnumcnt + 1][rbdim][nodenr] * co;
						z += fullBasisVectors[fnumcnt + 2][rbdim][nodenr] * co;
					}
					d.setDisplacement(nodenr, (float) x, (float) y, (float) z);
				}
				res.addField(d);
				break;
			}
			case Displacement2D: {
				DisplacementField d = new DisplacementField(sftype, numDoF);
				for (int nodenr = 0; nodenr < numDoF; nodenr++) {
					double x = 0, y = 0;
					for (int rbdim = 0; rbdim < N; rbdim++) {
						co = get_soln_coeff(rbdim);
						x += fullBasisVectors[fnumcnt][rbdim][nodenr] * co;
						y += fullBasisVectors[fnumcnt + 1][rbdim][nodenr] * co;
					}
					d.setDisplacement(nodenr, (float) x, (float) y);
				}
				res.addField(d);
				break;
			}
			case RealValue: {
				DefaultSolutionField f = new DefaultSolutionField(sftype, numDoF);
				for (int nodenr = 0; nodenr < numDoF; nodenr++) {
					double x = 0;
					for (int rbdim = 0; rbdim < N; rbdim++) {
						x += fullBasisVectors[fnumcnt][rbdim][nodenr] * get_soln_coeff(rbdim);
					}
					f.setValue(nodenr, (float) x);
				}
				res.addField(f);
				break;
			}
			default: {
				throw new RuntimeException("Invalid/unimplemented solution field type '" + sftype + "'");
			}
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

	/**
	 * Returns the increment used during the last parameter sweep or zero if
	 * there hasn't been any.
	 * 
	 * @return
	 */
	public double getSweepIncrement() {
		return sweepIncrement;
	}

	public double[][] getSweepOutputBounds() {
		return sweepOutputBounds;
	}

	public double[][] getSweepOutputs() {
		return sweepOutputs;
	}

	public SimulationResult getSweepSimResults() {
		int N = RB_sweep_solution[0][0].length;
		int numSweeps = RB_sweep_solution.length;

		double tmpval;
		SimulationResult res = new SimulationResult(numSweeps);

		int currentDoFField = 0;
		for (FieldDescriptor sftype : logicalFieldTypes) {
			if (currentDoFField + sftype.Type.requiredDoFFields > getNumDoFFields()) {
				throw new RuntimeException("Too many output fields used by current "
						+ "SolutionFieldTypes set in RBSystem. Check your model.xml.");
			}
			int numDoFs = fullBasisVectors[currentDoFField][0].length;
			Log.d("RBSystem", "Creating sweep solution field of type " + sftype.Type + ", sweeps:" + numSweeps
					+ ", Dofs: " + numDoFs);
			switch (sftype.Type) {
			case Displacement3D: {
				DisplacementField d = new DisplacementField(sftype, numDoFs * numSweeps);
				for (int iSweep = 0; iSweep < numSweeps; iSweep++) {
					for (int nodenr = 0; nodenr < numDoFs; nodenr++) {
						double x = 0, y = 0, z = 0;
						for (int rbdim = 0; rbdim < N; rbdim++) {
							x += fullBasisVectors[currentDoFField][rbdim][nodenr] * RB_sweep_solution[iSweep][0][rbdim];
							y += fullBasisVectors[currentDoFField + 1][rbdim][nodenr]
									* RB_sweep_solution[iSweep][0][rbdim];
							z += fullBasisVectors[currentDoFField + 2][rbdim][nodenr]
									* RB_sweep_solution[iSweep][0][rbdim];
						}
						d.setDisplacement(iSweep * numDoFs + nodenr, (float) x, (float) y, (float) z);
					}
				}
				res.addField(d);
				break;
			}
			case Displacement2D: {
				DisplacementField d = new DisplacementField(sftype, numDoFs * numSweeps);
				for (int iSweep = 0; iSweep < numSweeps; iSweep++) {
					for (int nodenr = 0; nodenr < numDoFs; nodenr++) {
						double x = 0, y = 0;
						for (int rbdim = 0; rbdim < N; rbdim++) {
							x += fullBasisVectors[currentDoFField][rbdim][nodenr] * RB_sweep_solution[iSweep][0][rbdim];
							y += fullBasisVectors[currentDoFField + 1][rbdim][nodenr]
									* RB_sweep_solution[iSweep][0][rbdim];
						}
						d.setDisplacement(iSweep * numDoFs + nodenr, (float) x, (float) y);
					}
				}
				res.addField(d);
				break;
			}
			case RealValue: {
				DefaultSolutionField f = new DefaultSolutionField(sftype, numSweeps * numDoFs);
				for (int sweepNr = 0; sweepNr < numSweeps; sweepNr++) {
					for (int i = 0; i < numDoFs; i++) {
						tmpval = 0;
						for (int j = 0; j < N; j++) {
							tmpval += fullBasisVectors[currentDoFField][j][i] * RB_sweep_solution[sweepNr][0][j];
						}
						f.setValue(sweepNr * numDoFs + i, (float) tmpval);
					}
				}
				res.addField(f);
				break;
			}
			default:
				throw new RuntimeException("Invalid/unimplemented solution field type '" + sftype.Type + "' for sweep.");
			}
		}
		for (MeshTransform m : transforms) {
			res.addTransform(m);
		}
		return res;
	}

	/**
	 * TODO check where this is used
	 * 
	 * @return
	 */
	public int getTotalTimesteps() {
		return 1;
	}

	/**
	 * Returns a float array of transformation data for each node, using the
	 * specified parameter mu.
	 * 
	 * Only to be called for models who have a parameterized geometry.
	 * 
	 * The get_local_transformation method of the AffineFunctions class is
	 * called (using the current parameter \mu) and the linear transform
	 * function returned.
	 * 
	 * @return
	 */
	protected float[][] getAffLinTransformationData(double[] mu) {
		assert transformationMethod != null;
		if (transformationMethod != null) {
			try {
				return (float[][]) transformationMethod.invoke(oldAffFcnObj, new Object[] { mu });
			} catch (IllegalAccessException iae) {
				throw new RuntimeException(iae);
			} catch (InvocationTargetException ite) {
				throw new RuntimeException(ite.getCause());
			}
		}
		return null;
	}

	public int getVisualNumTimesteps() {
		return 1;
	}

	/**
	 * @return
	 */
	protected boolean hasCustomMeshTransform() {
		Method meth;

		try {
			Class<?> partypes[] = null;
			meth = oldAffFcnCl.getMethod("is_custom_mesh_transform", partypes);
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

	protected boolean hasAffLinTransformationData() {
		return transformationMethod != null;
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
	@Override
	public void loadOfflineData(AModelManager m) throws MathReaderException, ModelManagerException, IOException {
		super.loadOfflineData(m);
		if (m.getModelType() == ModelType.rbappmit) {
			loadOfflineData_rbappmit(m);
		} else {
			loadOfflineDataJRB(m);
		}

		initialize_data_vectors();
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
			numBasisFuncs = Integer.parseInt(line);
			reader.close();
			reader = null;

			Log.d(DEBUG_TAG, "Finished reading n_bfs.dat");
		}

		/*
		 * Load zero initial conditions for old rbappmit models
		 */
		RB_initial_coeffs = new RealVector[1];
		RB_initial_coeffs[0] = new ArrayRealVector(numBasisFuncs);

		// Read in output data
		if (getNumOutputs() > 0) {
			// Get output dual norms
			{
				RB_output_vectors = new RealVector[getNumOutputs()][];
				output_dual_norms = new double[getNumOutputs()][];
				String[] dual_norms_tokens;
				for (int i = 0; i < getNumOutputs(); i++) {
					{
						BufferedReader reader = m
								.getBufReader("output_" + String.format("%03d", i) + "_dual_norms.dat");

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
							BufferedReader reader_i = m.getBufReader("output_" + String.format("%03d", i) + "_"
									+ String.format("%03d", q_l) + ".dat");

							String line_i = reader_i.readLine();
							reader_i.close();
							output_i_tokens = line_i.split(" ");
						}

						RB_output_vectors[i][q_l] = new ArrayRealVector(numBasisFuncs);
						for (int j = 0; j < numBasisFuncs; j++) {
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
					BufferedReader reader = m.getBufReader("RB_F_" + String.format("%03d", q_f) + ".dat");
					String line = reader.readLine();
					reader.close();
					reader = null;
					tokens = line.split(" ");
				}

				// Set the size of the inner product matrix
				RB_F_q_vector[q_f] = new ArrayRealVector(numBasisFuncs);

				// Fill the vector
				for (int i = 0; i < numBasisFuncs; i++) {
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
					BufferedReader reader = m.getBufReader("RB_A_" + String.format("%03d", q_a) + ".dat");
					String line = reader.readLine();
					reader.close();
					reader = null;
					tokens = line.split(" ");
				}

				// Set the size of the inner product matrix
				RB_A_q_vector[q_a] = new Array2DRowRealMatrix(numBasisFuncs, numBasisFuncs);

				// Fill the vector
				int count = 0;
				for (int i = 0; i < numBasisFuncs; i++)
					for (int j = 0; j < numBasisFuncs; j++) {
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
			Fq_Aq_representor_norms = new double[getQf()][getQa()][numBasisFuncs];

			// Fill it
			int count = 0;
			for (int q_f = 0; q_f < getQf(); q_f++)
				for (int q_a = 0; q_a < getQa(); q_a++)
					for (int i = 0; i < numBasisFuncs; i++) {
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
					file = "Aq_Aq_" + String.format("%03d", i) + "_" + String.format("%03d", j) + "_norms.bin";
					InputStream in = m.getInStream(file);
					try {
						Aq_Aq_representor_norms[count] = mr.readRawDoubleMatrix(in, numBasisFuncs, numBasisFuncs);
					} finally {
						in.close();
					}
					count++;
				}
			}
			Log.d(DEBUG_TAG, "Finished reading Aq_Aq_norms.dat");
		}

		// Read in Z data
		{
			int mf = getNumDoFFields();
			if (mf > 0) {
				fullBasisVectors = new float[mf][numBasisFuncs][];
				InputStream in;
				for (int imf = 0; imf < mf; imf++)
					for (int inbfs = 0; inbfs < numBasisFuncs; inbfs++) {
						in = m.getInStream("Z_" + String.format("%03d", imf) + "_" + String.format("%03d", inbfs)
								+ ".bin");
						try {
							fullBasisVectors[imf][inbfs] = mr.readRawFloatVector(in, getGeometry().getNumVertices());
						} finally {
							in.close();
						}
					}
			}

			Log.d(DEBUG_TAG, "Finished reading Z.dat");
		}

		// Reading uL data
		{
			if (getQuL() > 0) {
				uL_vector = new float[getQuL()][];
				InputStream in;
				for (int q_uL = 0; q_uL < getQuL(); q_uL++) {
					in = m.getInStream("uL_" + String.format("%03d", q_uL) + ".bin");
					try {
						uL_vector[q_uL] = mr.readRawFloatVector(in, getGeometry().getNumVertices());
					} finally {
						in.close();
					}
				}
			}
			Log.d(DEBUG_TAG, "Finished reading uL.dat");
		}

		// See if theres a transformation
		try {
			// Get a reference to get_n_L_functions, which does not
			// take any argument
			transformationMethod = oldAffFcnCl.getMethod("get_local_transformation", new Class[] { double[].class });
		} catch (NoSuchMethodException nsme) {
			transformationMethod = null;
		}

	}

	protected void loadOfflineDataJRB(AModelManager m) throws IOException {

		MathObjectReader mr = m.getMathObjReader();
		String filename;

		/*
		 * Get initial value coefficient vectors
		 */
		int Qu0 = 1;
		if (affineFunctionsInstance instanceof IAffineInitials) {
			Qu0 = ((IAffineInitials) affineFunctionsInstance).getQu0();
		}
		RB_initial_coeffs = new RealVector[Qu0];
		for (int i = 0; i < Qu0; i++) {
			filename = "RB_initial_" + String.format("%03d", i) + ".bin";
			RB_initial_coeffs[i] = mr.readVector(m.getInStream(filename));
		}
		Log.d(DEBUG_TAG, "Finished initial value data");

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
				filename = "output_" + String.format("%03d", i) + "_" + String.format("%03d", q_l) + ".bin";

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
		Fq_Aq_representor_norms = new double[fQf][getQa()][numBasisFuncs];
		for (int q_a = 0; q_a < getQa(); q_a++) {
			for (int q_f = 0; q_f < fQf; q_f++) {
				filename = "Fq_Aq_" + String.format("%03d", q_f) + "_" + String.format("%03d", q_a) + ".bin";
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
				filename = "Aq_Aq_" + String.format("%03d", i) + "_" + String.format("%03d", j) + "_norms.bin";
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
		 * Read in Z data
		 */
		if (getNumDoFFields() > 0) {
			fullBasisVectors = new float[getNumDoFFields()][numBasisFuncs][];
			for (int fieldNr = 0; fieldNr < getNumDoFFields(); fieldNr++)
				for (int bfNr = 0; bfNr < numBasisFuncs; bfNr++) {
					filename = "Z_" + String.format("%03d", fieldNr) + "_" + String.format("%03d", bfNr) + ".bin";
					fullBasisVectors[fieldNr][bfNr] = mr.readRawFloatVector(m.getInStream(filename));
				}
			Log.d(DEBUG_TAG, "Finished reading Z data");
		}
	}

	/**
	 * Custom mesh transformation method of old rbappmit-models
	 * 
	 * @param mu
	 * @param x
	 * @return
	 */
	public float[] rbappmitCustomMeshTransform(double[] mu, float[] x) {
		Method meth;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[2];
			partypes[0] = double[].class;
			partypes[1] = float[].class;

			meth = oldAffFcnCl.getMethod("mesh_transform", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for mesh_transform failed", nsme);
		}

		float[] xt;
		try {
			Object arglist[] = new Object[2];
			arglist[0] = mu;
			arglist[1] = x;

			Object theta_obj = meth.invoke(oldAffFcnObj, arglist);
			xt = (float[]) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return xt;
	}

	/**
	 * Performs a parameter sweep with automatic guess of the number of
	 * parameter sweeps.
	 * 
	 * @param sweepIndex
	 * @param N
	 * @return The number of sweeps performed
	 */
	public int performSweep(int sweepIndex, int N) {
		Log.d("RBSystem", "Starting sweep: sweepIdx:" + sweepIndex + ", N:" + N + ", NumDoFs:" + getNumDoFFields()
				+ ", vertices:" + getGeometry().getNumVertices());
		return performSweep(sweepIndex, N, Math.round(100000 / (getNumDoFFields() * getGeometry().getNumVertices())));
	}

	/**
	 * Performs a parameter sweep for the current RB system.
	 * 
	 * @param sweepIndex
	 *            The index of the parameter to sweep
	 * @param N
	 *            The reduced basis size for sweeping
	 * @param numSweepPts
	 *            The desired number of sweep points
	 * @return The number of sweeps effectively performed
	 */
	public int performSweep(int sweepIndex, int N, int numSweepPts) {
		Parameters p = getParams();
		if (!isReal) {
			Log.d("RBSystem", "Dividing numSweepPts by three.");
			numSweepPts /= 3;
		}
		if (numSweepPts > MAX_SWEEP_POINTS) {
			Log.e("RBSystem", "Too large number of sweep points, allowed: " + MAX_SWEEP_POINTS + ", requested: "
					+ numSweepPts);
			numSweepPts = MAX_SWEEP_POINTS;
		}
		int n_outputs = getNumOutputs();

		mSweepParam = new double[numSweepPts][];

		RB_sweep_solution = null;
		if (isReal) {
			RB_sweep_solution = new double[numSweepPts][1][current_N];
		} else {
			RB_sweep_solution = new double[numSweepPts][2][current_N];
			n_outputs *= 2;
		}

		sweepOutputs = new double[n_outputs][numSweepPts];
		sweepOutputBounds = new double[n_outputs][numSweepPts];

		double sweepParamRange = p.getMaxValue(sweepIndex) - p.getMinValue(sweepIndex);
		sweepIncrement = sweepParamRange / (numSweepPts - 1);

		transforms.clear();
		for (int sweepNr = 0; sweepNr < numSweepPts; sweepNr++) {
			double new_param = p.getMinValue(sweepIndex) + sweepNr * sweepIncrement;
			mSweepParam[sweepNr] = p.getCurrent().clone();
			Log.d(DEBUG_TAG, "ParamSweep: New param " + Arrays.toString(p.getCurrent()));
			p.setCurrent(sweepIndex, new_param);
			RB_solve(N);

			if (isReal) {
				for (int n = 0; n < n_outputs; n++) {
					sweepOutputs[n][sweepNr] = RB_outputs[n];
					sweepOutputBounds[n][sweepNr] = RB_output_error_bounds[n];
				}
			} else {
				for (int n = 0; n < n_outputs / 2; n++) {
					sweepOutputs[n][sweepNr] = get_RB_output(n, true);
					sweepOutputs[n + n_outputs / 2][sweepNr] = get_RB_output(n, false);
					sweepOutputBounds[n][sweepNr] = get_RB_output_error_bound(n, true);
					sweepOutputBounds[n + n_outputs / 2][sweepNr] = get_RB_output_error_bound(n, false);
				}
			}

			if (getNumDoFFields() > 0) {
				RB_sweep_solution[sweepNr] = get_RBsolution();
				if (!hasCustomMeshTransform()) {
					if (hasAffLinTransformationData()) {
						float[][] tmp = getAffLinTransformationData(p.getCurrent());
						Log.d("RBSystem", "Storing AffLinTrafo: " + Log.dumpArr(tmp));
						transforms.add(new AffineLinearMeshTransform(tmp, getGeometry().vertexLTFuncNr));
					}
				} else {
					transforms.add(new rbappmitCustomMeshTransform(mSweepParam[sweepNr], this));
				}
			}
		}
		return numSweepPts;
	}

	/**
	 * Perform online solve with the N RB basis functions, for the set of
	 * parameters in current_params, where 1 <= N <= RB_size.
	 */
	protected double RB_solve(int N) {

		current_N = N;

		if (N > getNBF()) {
			throw new RuntimeException("ERROR: N cannot be larger than the number " + "of basis functions in RB_solve");
		}
		if (N == 0) {
			throw new RuntimeException("ERROR: N must be greater than 0 in RB_solve");
		}

		// Assemble the RB system
		RealMatrix RB_system_matrix_N = new Array2DRowRealMatrix(N, N);

		for (int q_a = 0; q_a < getQa(); q_a++) {
			RB_system_matrix_N = RB_system_matrix_N.add(RB_A_q_vector[q_a].getSubMatrix(0, N - 1, 0, N - 1)
					.scalarMultiply(thetaQa(q_a)));
		}

		// Assemble the RB rhs
		RealVector RB_rhs_N = getInitialCoefficients(N);

		for (int q_f = 0; q_f < fQf; q_f++) {
			// Note getSubVector takes an initial index and the number of
			// entries
			// i.e. the interface is a bit different to getSubMatrix
			RB_rhs_N = RB_rhs_N.add(RB_F_q_vector[q_f].getSubVector(0, N).mapMultiplyToSelf(thetaQf(q_f)));
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
				RB_outputs[i] += thetaQl(i, q_l) * RB_solution.dotProduct(RB_output_vector_N);
			}
			RB_output_error_bounds[i] = compute_output_dual_norm(i, 0) // Zero
																		// is
																		// current
																		// time,
																		// not
																		// in
																		// RBSystem
					* abs_error_bound;
		}

		return (return_rel_error_bound ? abs_error_bound / RB_solution_norm : abs_error_bound);
	}

	/**
	 * 
	 * @param m
	 * @return
	 */
	public final boolean readConfiguration(AModelManager m) {
		fQa = affineFunctionsInstance.getQa();
		Log.d("RBBase", "Q_a = " + fQa);

		if (m.getModelType() == ModelType.rbappmit) {
			try {
				GetPot infile = new GetPot(m.getInStream(Const.parameters_filename), Const.parameters_filename);
				Log.d("RBBase", "Created GetPot object");
				readConfigurationRBAppMIT(infile);
			} catch (IOException e) {
				Log.e("RBBase", "Error opening input.in", e);
				return false;
			}
		} else {
			readConfigurationJRB(m);
		}
		return true;
	}

	protected void readConfigurationJRB(AModelManager m) {
		params = m.getParameters();

		// Number of basis functions
		numBasisFuncs = Integer.parseInt(m.getModelXMLAttribute("num_basisfcn", "rb_model"));

		// Number of output functionals
		fNumOutputs = affineFunctionsInstance.getNumOutputs();

		Ql_values = affineFunctionsInstance.getQl();
		fQf = affineFunctionsInstance.getQf();
		if (affineFunctionsInstance instanceof IWithuL) {
			fQuL = ((IWithuL) affineFunctionsInstance).getQuL();
		}
	}

	protected void readConfigurationRBAppMIT(GetPot infile) {

		int n_parameters = infile.call("n_parameters", 1);
		Log.d("RBBase", "n_parameters = " + n_parameters);
		if (n_parameters > 0) {
			params = new Parameters();
			for (int i = 0; i < n_parameters; i++) {
				String label = infile.call("param" + Integer.toString(i) + "_label", "mu_" + i);
				// Read in the min/max for the i^th parameter
				String min_string = new String("mu" + i + "_min");
				double mu_i_min = infile.call(min_string, 0.);
				String max_string = new String("mu" + i + "_max");
				double mu_i_max = infile.call(max_string, 0.);
				params.addParam(label, mu_i_min, mu_i_max);
				Log.d("RBBase", "Parameter " + i + ": Min = " + mu_i_min + ", Max = " + mu_i_max);
			}
		}

		/*
		 * Read config values from the input.in file
		 */
		Log.d(DEBUG_TAG, "Entered parse_parameters_file, filename = " + Const.parameters_filename);

		int return_rel_error_bound_in = infile.call("return_rel_error_bound", 1);
		return_rel_error_bound = (return_rel_error_bound_in != 0);

		Log.d(DEBUG_TAG, "return a relative error bound from RB_solve? " + return_rel_error_bound);

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

	public void set_n_basis_functions(int _N) {
		numBasisFuncs = _N;
	}

	/**
	 * Only used to set the current parameters object for an SCMSystem, as this
	 * also inherits from the RBBase class.
	 * 
	 * @param p
	 */
	protected void setParams(Parameters p) {
		params = p;
	}

	/**
	 * Set the primary SCM
	 */
	public void setPrimarySCM(RBSCMSystem scm_system) {
		mRbScmSystem = scm_system;
	}

	public double solveRB(int N) {
		/**
		 * If this system is unsteady dont perform geometric transformations
		 * (not checked yet/no models there yet) Comes from the LINEAR_STEADY /
		 * LINEAR_COMPLEX_STEADY enums in ROMSim/rbappmit, where the check was
		 * made originally.
		 */
		transforms.clear();
		if (!(this instanceof TransientRBSystem)) {
			if (!hasCustomMeshTransform()) {
				if (hasAffLinTransformationData()) {
					float[][] tmp = getAffLinTransformationData(getParams().getCurrent());
					Log.d("RBSystem", "Storing AffLinTrafo: " + Log.dumpArr(tmp));
					transforms.add(new AffineLinearMeshTransform(tmp, getGeometry().vertexLTFuncNr));
				}
			} else {
				transforms.add(new rbappmitCustomMeshTransform(getParams().getCurrent(), this));
			}
		} else {
			transforms.add(new DefaultTransform());
		}
		return RB_solve(N);
	}

	// private void sweepUpdateGeometry(float[][][] RB_sweep_LTfuncs, int
	// numSweepPts) {
	// GeometryData geoData = getGeometry();
	// if (!hasCustomMeshTransform()) {
	// /*
	// * Apply affine-linear transformation according to sweeped parameter
	// * values if present
	// */
	// if (RB_sweep_LTfuncs != null) {
	// /*
	// * The number of sweeps os encoded in the length of the first
	// * dimension of the RB_sweep_LTfuncs array as there is a set for
	// * each sweep.
	// */
	// geoData.applyAffLinVertexTransformation(RB_sweep_LTfuncs);
	// /*
	// * otherwise just copy the nodes as often as there have been
	// * parameter sweeps
	// *
	// * @TODO separate geometry / field value frames for less memory
	// * usage!
	// */
	// } else {
	// geoData.vertices = new float[numSweepPts * geoData.getNumVertices() * 3];
	// for (int sweep = 0; sweep < numSweepPts; sweep++) {
	// System.arraycopy(geoData.originalVertices, 0, geoData.vertices, sweep *
	// geoData.getNumVertices()
	// * 3, geoData.getNumVertices());
	// }
	// }
	// } else {
	// customMeshTransform(mSweepParam, geoData);
	// }
	// }

	/**
	 * Evaluate theta_q_a (for the q^th bilinear form) at the current parameter.
	 * 
	 * @param q
	 * @return
	 */
	public double thetaQa(int q) {
		return thetaQa(q, 0);
	}

	/**
	 * 
	 * @param q
	 * @param t
	 * @return
	 */
	public double thetaQa(int q, double t) {
		return affineFunctionsInstance.thetaQa(q, params.getCurrent(), t);
	}

	/**
	 * 
	 * @param i
	 * @return
	 */
	public double thetaQf(int i) {
		return thetaQf(i, 0);
	}

	/**
	 * 
	 * @param i
	 * @param t
	 * @return
	 */
	public double thetaQf(int i, double t) {
		return affineFunctionsInstance.thetaQf(i, params.getCurrent(), t);
	}

	/**
	 * 
	 * @param k
	 * @param i
	 * @return
	 */
	public double thetaQl(int k, int i) {
		return thetaQl(k, i, 0);
	}

	/**
	 * 
	 * @param k
	 * @param i
	 * @param t
	 * @return
	 */
	public double thetaQl(int k, int i, double t) {
		return affineFunctionsInstance.thetaQl(k, i, getParams().getCurrent(), t);
	}
}
