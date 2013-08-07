package rb;

import jarmos.Log;
import jarmos.io.AModelManager;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.linear.RealVector;

/**
 * Base class for quadratically nonlinear RB systems including the SCM method for error bound computation.
 * 
 * This class has been taken from the original @ref rbappmit package and modified to fit into the current JaRMoS
 * framework.
 * 
 * @author Daniel Wirtz @date 07.08.2013
 * 
 */
public class QNTransientSCMSystem extends RBSCMSystem {

	public QNTransientSCMSystem(RBSystem sys) {
		super(sys);
	}

	// Logging tag
	private static final String DEBUG_TAG = "QNTransientSCMSystem";

	/**
	 * The current number of basis functions.
	 */
	private int n_bfs;

	/**
	 * The RB coefficients at the C_J parameters.
	 */
	private double[][] C_J_RB_coeffs;

	/**
	 * The current RB coefficientss.
	 */
	private RealVector current_RB_coeffs;

	// We also need to save the RB coefficients during LB calculations
	private RealVector saved_RB_coeffs;

	// We may need to scale the theta_c function for the sake of the SCM!
	private double SCM_theta_c_scaling;

	/**
	 * Get/set the number of basis functions
	 */
	int get_n_basis_functions() {
		return n_bfs;
	}

	void set_n_basis_functions(int n_bfs_in) {
		n_bfs = n_bfs_in;
	}

	/**
	 * Set the current RB coefficients.
	 */
	void set_current_RB_coeffs(RealVector RB_coeffs) {
		current_RB_coeffs = RB_coeffs;
	}

	/**
	 * Evaluate theta_c (for the quadratic nonlinearity) at the current parameter.
	 */
	public double eval_theta_c() {

		Method meth;

		try {
			Class<?> partypes[] = new Class[1];
			partypes[0] = double[].class;

			meth = sys.oldAffFcnCl.getMethod("evaluateC", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateC failed", nsme);
		}

		Double theta_val;
		try {
			Object arglist[] = new Object[1];
			arglist[0] = sys.getParams().getCurrent();

			Object theta_obj = meth.invoke(sys.oldAffFcnObj, arglist);
			theta_val = (Double) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return SCM_theta_c_scaling * theta_val.doubleValue();
	}

	/**
	 * Override get_Q_a since we have n_bfs extra terms
	 */
	@Override
	public int getQa() {
		return super.getQa() + get_n_basis_functions();
	}

	/**
	 * Override eval_theta_q_a in order to account for the affine terms related to basis functions
	 */
	@Override
	public double thetaQa(int q) {

		if (q < get_n_basis_functions()) {
			double theta_c_value = eval_theta_c();
			return current_RB_coeffs.getEntry(q) * theta_c_value;
		} else {

			Method meth;

			try {
				Class<?> partypes[] = new Class[2];
				partypes[0] = Integer.TYPE;
				partypes[1] = double[].class;

				meth = sys.oldAffFcnCl.getMethod("evaluateA", partypes);
			} catch (NoSuchMethodException nsme) {
				throw new RuntimeException("getMethod for evaluateA failed", nsme);
			}

			Double theta_val;
			try {
				Object arglist[] = new Object[2];
				arglist[0] = new Integer(q - get_n_basis_functions());
				arglist[1] = sys.getParams().getCurrent();

				Object theta_obj = meth.invoke(sys.oldAffFcnObj, arglist);
				theta_val = (Double) theta_obj;
			} catch (IllegalAccessException iae) {
				throw new RuntimeException(iae);
			} catch (InvocationTargetException ite) {
				throw new RuntimeException(ite.getCause());
			}

			return theta_val.doubleValue();
		}
	}

	/**
	 * Override to also load the RB coefficients.
	 */
	@Override
	protected void get_current_parameters_from_C_J(int index) {
		super.get_current_parameters_from_C_J(index);

		for (int i = 0; i < get_n_basis_functions(); i++)
			current_RB_coeffs.setEntry(i, C_J_RB_coeffs[index][i]);
	}

	/**
	 * Override to also save the RB coefficients.
	 */
	@Override
	protected void save_current_parameters() {
		super.save_current_parameters();

		saved_RB_coeffs = current_RB_coeffs.copy();
	}

	/**
	 * Override to also load the RB coefficients.
	 */
	@Override
	protected void reload_current_parameters() {
		super.reload_current_parameters();

		set_current_RB_coeffs(saved_RB_coeffs);
	}

	/**
	 * Override read_offline_data in order to read in the extra data in the QNTransient case.
	 */
	@Override
	public void loadOfflineData(AModelManager m) throws IOException, InconsistentStateException {

		// Initially set number of basis functions from n_bfs.dat
		{
			BufferedReader reader = m.getBufReader("n_bfs.dat");

			String line = reader.readLine();
			reader.close();
			reader = null;

			n_bfs = Integer.parseInt(line);

			Log.d(DEBUG_TAG, "Finished reading n_bfs.dat");
		}

		super.loadOfflineData(m);

		// Read C_J_RB_coeffs
		{
			C_J_RB_coeffs = new double[C_J_stability_vector.length][get_n_basis_functions()];
			if (C_J_stability_vector != null) {
				String[] tokens;
				{
					BufferedReader reader = m.getBufReader("C_J_RB_coeffs.dat");
					String line = reader.readLine();
					reader.close();
					reader = null;
					tokens = line.split(" ");
				}

				int count = 0;
				for (int i = 0; i < C_J_stability_vector.length; i++) {
					for (int j = 0; j < get_n_basis_functions(); j++) {
						C_J_RB_coeffs[i][j] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			}
			Log.d(DEBUG_TAG, "Finished reading C_J_RB_coeffs.dat");
		}
	}

	/**
	 * @param parameters_filename
	 * The name of the file to parse Parse the input file to initialize this RBSCMSystem.
	 * @throws IOException
	 */
	@Override
	public void readConfiguration(AModelManager m) throws IOException {
		super.readConfiguration(m);
		GetPot infile = new GetPot(m.getInStream(Const.parameters_filename), Const.parameters_filename);
		SCM_theta_c_scaling = infile.call("SCM_theta_c_scaling", 1.);
		Log.d(DEBUG_TAG, "SCM_theta_c_scaling = " + SCM_theta_c_scaling);
	}

}
