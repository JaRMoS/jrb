package rb.java;

import jarmos.Log;
import jarmos.io.AModelManager;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.linear.LinearConstraint;
import org.apache.commons.math.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math.optimization.linear.Relationship;
import org.apache.commons.math.optimization.linear.SimplexSolver;



// This class implements the Online stage
// of the Successive Constraint Method
// for coercive problems.
// This class is modeled on the RBSCMSystem
// class in rbOOmit

/**
 * Changes made by
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 *
 */
public class RBSCMSystem {
	
	protected RBSystem sys;

	// Logging tag
	private static final String DEBUG_TAG = "RBSCMSystem";

	/**
	 * The maximum number of constraints we impose.
	 */
	private int SCM_M;

	/**
	 * The bounding box values
	 */
	private double[] B_min;
	private double[] B_max;

	/**
	 * The Greedily selected parameters.
	 */
	private Vector<double[]> C_J;

	/**
	 * The values of the stability factor at the greedily selected parameters.
	 */
	protected double[] C_J_stability_vector;

	/**
	 * This matrix stores the infimizing vectors y_1(\mu),...,y_Q_a(\mu), for
	 * each \mu in C_J, which are used in computing the SCM upper bounds.
	 */
	private double[][] SCM_UB_vectors;

	/**
	 * A private Parameter used to temporarily store current_parameters during
	 * the SCM calculation
	 */
	private double[] saved_parameters;
	
	public RBSCMSystem(RBSystem sys) {
		this.sys = sys;
	}
	
	protected double thetaQa(int q) {
		return sys.thetaQa(q);
	}
	
	protected int getQa() {
		return sys.getQa();
	}

	/**
	 * @return the SCM lower bound for the current parameters.
	 */
	public double get_SCM_LB() {
		double min_J_obj = 0.;

		try {

			// First, declare the constraints
			Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();

			// Add bounding box constraints for the get_Q_a() variables
			for (int q = 0; q < getQa(); q++) {
				double[] index = new double[getQa()];
				index[q] = 1.;

				constraints.add(new LinearConstraint(index, Relationship.GEQ,
						B_min[q]));
				constraints.add(new LinearConstraint(index, Relationship.LEQ,
						B_max[q]));
			}

			// Sort the indices of C_J based on distance from current_parameters
			List<Integer> sortedIndices = getSorted_CJ_Indices();

			// Save the current_parameters since we'll change them in the loop
			// below
			save_current_parameters();

			// Add the constraint rows
			int n_rows = Math.min(SCM_M, C_J.size());
			int count = 1;

			if (n_rows > 0) {
				for (Integer mu_index : sortedIndices) {
					get_current_parameters_from_C_J(mu_index);
					// current_parameters = C_J.get(mu_index);

					double[] constraint_row = new double[getQa()];
					for (int q = 0; q < getQa(); q++) {
						constraint_row[q] = sys.thetaQa(q);
					}

					constraints.add(new LinearConstraint(constraint_row,
							Relationship.GEQ, C_J_stability_vector[mu_index]));

					if (count >= n_rows)
						break;

					count++;
				}
			}

			// Now load the original parameters back into current_parameters
			// in order to set the coefficients of the objective function
			reload_current_parameters();

			// Create objective function object
			double[] objectiveFn = new double[getQa()];
			for (int q = 0; q < getQa(); q++) {
				objectiveFn[q] = sys.thetaQa(q);
			}
			LinearObjectiveFunction f = new LinearObjectiveFunction(
					objectiveFn, 0.);

			SimplexSolver solver = new SimplexSolver();
			RealPointValuePair opt_pair = solver.optimize(f, constraints,
					GoalType.MINIMIZE, false);
			min_J_obj = opt_pair.getValue();
		} catch (OptimizationException e) {
			Log.e("DEBUG_TAG", "Optimal solution not found");
			e.printStackTrace();
		} catch (Exception e) {
			Log.e("DEBUG_TAG", "Exception occurred during SCM_LB calculation");
			e.printStackTrace();
		}

		Log.d(DEBUG_TAG, "SCM val = " + min_J_obj);
		return min_J_obj;
	}

	/**
	 * Evaluate the SCM upper bound for current_parameters.
	 */
	public double get_SCM_UB() {

		// Sort the indices of C_J based on distance from current_parameters
		List<Integer> sortedIndices = getSorted_CJ_Indices();

		// For each mu, we just find the minimum of J_obj over
		// the subset of vectors in SCM_UB_vectors corresponding
		// to C_J_M (SCM_UB_vectors contains vectors for all of
		// C_J).
		double min_J_obj = 0.;
		int n_rows = Math.min(SCM_M, C_J.size());
		int count = 1;
		for (Iterator<Integer> it = sortedIndices.iterator(); it.hasNext();) {
			Integer mu_index = (Integer) it.next();

			get_current_parameters_from_C_J(mu_index);

			double[] UB_vector = SCM_UB_vectors[mu_index];

			double J_obj = 0.;
			for (int q = 0; q < getQa(); q++) {
				J_obj += sys.thetaQa(q) * UB_vector[q];
			}

			if ((count == 1) || (J_obj < min_J_obj)) {
				min_J_obj = J_obj;
			}

			if (count >= n_rows)
				break;

			count++;
		}

		return min_J_obj;
	}

	/**
	 * @return Euclidean distance between two parameters
	 */
	public static double param_dist(double[] mu_1, double[] mu_2) {
		// Default distance is Euclidean norm
		double sum = 0.;

		for (int i = 0; i < mu_1.length; i++) {
			sum += Math.pow(mu_1[i] - mu_2[i], 2.);
		}

		return Math.sqrt(sum);
	}

	/**
	 * Load the current_parameters from the set C_J.
	 */
	protected void get_current_parameters_from_C_J(int index) {
		sys.getParams().setCurrent(C_J.get(index));
	}

	/**
	 * Save current_parameters in saved_parameters.
	 */
	protected void save_current_parameters() {
		saved_parameters = sys.getParams().getCurrent().clone();
	}

	/**
	 * Reload from saved_parameters
	 */
	protected void reload_current_parameters() {
		sys.getParams().setCurrent(saved_parameters);
	}

	/**
	 * @return the current parameters
	 */
	public double[] get_current_parameters() {
		return sys.getParams().getCurrent();
	}

//	@Override
//	protected void readConfigurationJRB(AModelManager m) {
//		super.readConfigurationJRB(m);
//	}
	
	/**
	 * 
	 * @param m
	 * @return
	 * @throws IOException 
	 */
	public void readConfiguration(AModelManager m) throws IOException {
		GetPot infile = new GetPot(m.getInStream(Const.parameters_filename), Const.parameters_filename);
		// int n_SCM_parameters = infile.call("n_SCM_parameters",1);
		int n_SCM_parameters = infile.call("n_SCM_parameters",
				infile.call("n_parameters", 1));
		Log.d(DEBUG_TAG, "n_parameters = " + n_SCM_parameters);

		SCM_M = infile.call("SCM_M", 0);

		Log.d(DEBUG_TAG, "RBSCMSystem parameters from " + Const.parameters_filename
				+ ":");
		Log.d(DEBUG_TAG, "SCM_M: " + SCM_M);
	}

	/**
	 * Read in the stored data from the specified URL in order to initialize the
	 * SCM.
	 * 
	 * @param directory_name
	 *            The URL of the directory containing the Offline data Read in
	 *            the Offline data to initialize this RBSystem.
	 */
	public void loadOfflineData(AModelManager m) throws IOException,
			InconsistentStateException {

		// Read in the bounding box minimum values
		{
			BufferedReader reader = m.getBufReader("B_min.dat");

			String line = reader.readLine();
			String[] tokens = line.split(" ");
			reader.close(); reader = null;

			B_min = new double[getQa()];
			for (int i = 0; i < B_min.length; i++) {
				B_min[i] = Double.parseDouble(tokens[i]);
			}

			Log.d(DEBUG_TAG, "Finished reading B_min.dat");
		}

		// Read in the bounding box maximum values
		{
			BufferedReader reader = m.getBufReader("B_max.dat");

			String line = reader.readLine();
			String[] tokens = line.split(" ");

			B_max = new double[getQa()];
			for (int i = 0; i < B_max.length; i++) {
				B_max[i] = Double.parseDouble(tokens[i]);
			}
			reader.close(); reader = null;

			Log.d(DEBUG_TAG, "Finished reading B_max.dat");
		}

		// Read in the stability constant values
		{
			BufferedReader reader = m.getBufReader("C_J_stability_vector.dat");

			String line = reader.readLine();
			reader.close(); reader = null;

			try {
				String[] tokens = line.split(" ");

				if ((tokens.length == 1) && (tokens[0] == "")) {
					C_J_stability_vector = null;
				} else {
					C_J_stability_vector = new double[tokens.length];
					for (int i = 0; i < C_J_stability_vector.length; i++) {
						C_J_stability_vector[i] = Double.parseDouble(tokens[i]);
					}
				}
			} catch (Exception e) {
				Log.d(DEBUG_TAG, "Exception occurred when splitting string, "
						+ "setting C_J_stability_vector to null");
				C_J_stability_vector = null;
			}

			Log.d(DEBUG_TAG, "Finished reading C_J_stability_vector.dat");
		}

		// Read in C_J, the Greedily selected parameters
		{
			BufferedReader reader = m.getBufReader("C_J.dat");

			C_J = new Vector<double[]>(0);
			if (C_J_stability_vector != null) {

				String line = reader.readLine();
				reader.close(); reader = null;
				String[] tokens = line.split(" ");

				int count = 0;
				int np = sys.getParams().getNumParams();
				for (int i = 0; i < C_J_stability_vector.length; i++) {
					C_J.add(new double[np]);
					for (int j = 0; j < np; j++) {
						C_J.get(i)[j] = Double.parseDouble(tokens[count]);
						count++;
					}
				}
			}

			Log.d(DEBUG_TAG, "Finished reading C_J.dat");
		}

		// Read in SCM_UB_vectors
		{
			BufferedReader reader = m.getBufReader("SCM_UB_vectors.dat");

			if (C_J_stability_vector != null) {

				String line = reader.readLine();
				reader.close(); reader = null;
				String[] tokens = line.split(" ");

				int count = 0;
				// Resize SCM_UB_vectors based on C_J_stability_vector and Q_a
				SCM_UB_vectors = new double[C_J_stability_vector.length][getQa()];
				for (int i = 0; i < SCM_UB_vectors.length; i++) {
					for (int j = 0; j < getQa(); j++) {
						SCM_UB_vectors[i][j] = Double
								.parseDouble(tokens[count]);
						count++;
					}
				}
			}

			Log.d(DEBUG_TAG, "Finished reading SCM_UB_vectors.dat");
		}

	}

	/**
	 * Private helper function to sort the indices of C_J based on distance from
	 * current_parameters
	 */
	private List<Integer> getSorted_CJ_Indices() {

		int J = C_J.size();

		LinkedHashMap<Double, Integer> dist_from_mu = new LinkedHashMap<Double, Integer>(
				J);

		for (int j = 0; j < J; j++) {
			double dist = param_dist(get_current_parameters(), C_J.get(j));
			dist_from_mu.put(dist, j);
		}

		List<Map.Entry<Double, Integer>> list = new LinkedList<Map.Entry<Double, Integer>>(
				dist_from_mu.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<Double, Integer>>() {
			public int compare(Map.Entry<Double, Integer> o1,
					Map.Entry<Double, Integer> o2) {
				return o1.getKey().compareTo(o2.getKey());
				/*
				 * return ((Comparable<?>) ((Map.Entry) (o1)).getKey())
				 * .compareTo(((Map.Entry) (o2)).getKey());
				 */
			}
		});

		// Create a sorted list of values to return
		List<Integer> result = new LinkedList<Integer>();
		for (Map.Entry<Double, Integer> e : list) {
			result.add(e.getValue());
		}
		return result;
	}

}
