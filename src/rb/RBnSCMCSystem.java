package rb;

import jarmos.Log;
import jarmos.io.AModelManager;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.linear.LinearConstraint;
import org.apache.commons.math.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math.optimization.linear.Relationship;
import org.apache.commons.math.optimization.linear.SimplexSolver;

/**
 * Reduced basis SCM system
 * 
 * This class has been taken from the original @ref rbappmit package and modified to fit into the current JaRMoS
 * framework.
 * 
 * @author Daniel Wirtz @date 07.08.2013
 * 
 */
public class RBnSCMCSystem extends RBSCMSystem {

	public RBnSCMCSystem(RBSystem sys) {
		super(sys);
	}

	// Logging tag
	private static final String DEBUG_TAG = "RBnSCMCSystem";

	private double[] B_min;
	private double[] B_max;
	private int n_mubar;
	private int[] n_muhat;
	private Vector<double[]> mu_bar;
	private Vector<double[]>[] mu_hat;
	private double[] beta_bar;
	private double[][] beta_hat;
	private double[][][] zval;

	/**
	 * @param m
	 * @throws IOException
	 * @throws InconsistentStateException
	 */
	@SuppressWarnings("unchecked")
	public void loadOfflineData(AModelManager m) throws IOException, InconsistentStateException {

		BufferedReader reader = m.getBufReader("SCMdata.dat");

		String line = reader.readLine();
		String[] tokens = line.split(" ");

		int count = 0;

		B_min = new double[getQa()];
		B_max = new double[getQa()];
		for (int i = 0; i < B_min.length; i++) {
			B_max[i] = Double.parseDouble(tokens[count]);
			B_min[i] = -B_max[i];
			count++;
		}

		n_mubar = Integer.parseInt(tokens[count]);
		count++;

		int np = sys.getParams().getNumParams();
		mu_bar = new Vector<double[]>(0);
		for (int i = 0; i < n_mubar; i++) {
			mu_bar.add(new double[np]);
			for (int j = 0; j < np; j++) {
				mu_bar.get(i)[j] = Double.parseDouble(tokens[count]);
				count++;
			}
		}

		beta_bar = new double[n_mubar];
		for (int i = 0; i < n_mubar; i++) {
			beta_bar[i] = Double.parseDouble(tokens[count]);
			count++;
		}

		mu_hat = (Vector<double[]>[]) new Vector<?>[n_mubar];
		n_muhat = new int[n_mubar];
		beta_hat = new double[n_mubar][];
		zval = new double[n_mubar][][];
		for (int i = 0; i < n_mubar; i++) {
			n_muhat[i] = Integer.parseInt(tokens[count]);
			count++;
			beta_hat[i] = new double[n_muhat[i]];
			zval[i] = new double[n_muhat[i]][getQa() * 2];

			mu_hat[i] = new Vector<double[]>(0);
			for (int j = 0; j < n_muhat[i]; j++) {
				mu_hat[i].add(new double[np]);
				for (int k = 0; k < np; k++) {
					mu_hat[i].get(j)[k] = Double.parseDouble(tokens[count]);
					count++;
				}
			}

			for (int j = 0; j < n_muhat[i]; j++) {
				beta_hat[i][j] = Double.parseDouble(tokens[count]);
				count++;
			}

			for (int k = 0; k < getQa() * 2; k++)
				for (int j = 0; j < n_muhat[i]; j++) {
					zval[i][j][k] = Double.parseDouble(tokens[count]);
					count++;
				}
		}

		reader.close();

		Log.d(DEBUG_TAG, "Finished reading SCMdata.dat");

	}

	/**
	 * TODO create local double[] currentParam field and remove save_/restore_params fcns
	 */
	public double get_SCM_LB() {
		// return 0.01;

		double min_J_obj = 0.;
		double[] min_Jlocal_obj = new double[n_mubar];

		// Sort the indices of mu_bar based on distance from current_parameters
		List<Integer> sortedmubarIndices = getSorted_CJ_Indices(mu_bar);
		int icount = 0;
		// while ((min_J_obj<=0) && (icount < sortedmubarIndices.size())){
		while ((min_J_obj <= 0) && (icount < sortedmubarIndices.size())) {
			int imubar = sortedmubarIndices.get(icount);

			// First, declare the constraints
			Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();

			// Add bounding box constraints for the get_Q_a() variables
			for (int q = 0; q < getQa(); q++) {
				double[] index = new double[getQa() * 2];
				index[q] = 1.;

				constraints.add(new LinearConstraint(index, Relationship.GEQ, B_min[q] / beta_bar[imubar]));
				constraints.add(new LinearConstraint(index, Relationship.LEQ, B_max[q] / beta_bar[imubar]));

				index[q] = 0.;
				index[q + getQa()] = 1.;

				constraints.add(new LinearConstraint(index, Relationship.GEQ, B_min[q] / beta_bar[imubar]));
				constraints.add(new LinearConstraint(index, Relationship.LEQ, B_max[q] / beta_bar[imubar]));
			}

			// Save the current_parameters since we'll change them in the loop
			// below
			save_current_parameters();

			// Add the constraint rows
			if (n_muhat[imubar] > 0) {
				for (int imuhat = 0; imuhat < n_muhat[imubar]; imuhat++) {
					sys.getParams().setCurrent(mu_hat[imubar].get(imuhat));

					double[] constraint_row = new double[getQa() * 2];
					for (int q = 0; q < getQa(); q++) {
						Complex theta_q_a = sys.complex_eval_theta_q_a(q);
						constraint_row[q] = theta_q_a.getReal() * beta_bar[imubar];
						constraint_row[q + getQa()] = theta_q_a.getImaginary() * beta_bar[imubar];
					}

					constraints.add(new LinearConstraint(constraint_row, Relationship.GEQ, beta_hat[imubar][imuhat]));
				}
			}

			// Now load the original parameters back into current_parameters
			// in order to set the coefficients of the objective function
			reload_current_parameters();

			// Create objective function object
			double[] objectiveFn = new double[getQa() * 2];
			for (int q = 0; q < getQa(); q++) {
				Complex theta_q_a = sys.complex_eval_theta_q_a(q);
				objectiveFn[q] = theta_q_a.getReal() * beta_bar[imubar];
				objectiveFn[q + getQa()] = theta_q_a.getImaginary() * beta_bar[imubar];
			}
			LinearObjectiveFunction f = new LinearObjectiveFunction(objectiveFn, 0.);

			try {
				SimplexSolver solver = new SimplexSolver(1e-6);
				RealPointValuePair opt_pair = solver.optimize(f, constraints, GoalType.MINIMIZE, false);
				min_Jlocal_obj[icount] = opt_pair.getValue();
			} catch (OptimizationException e) {
				Log.e("RBSCMSYSTEM_TAG", "Optimal solution not found");
				e.printStackTrace();
			} catch (Exception e) {
				Log.e("RBSCMSYSTEM_TAG", "Exception occurred during SCM_LB calculation");
				e.printStackTrace();
			}

			min_J_obj = min_J_obj > min_Jlocal_obj[icount] ? min_J_obj : min_Jlocal_obj[icount];
			icount++;
		}
		return min_J_obj;
	}

	public double get_SCM_UB() {

		// cheating, since betaUB is rarely good
		double maxbetabar = beta_bar[0];
		for (int i = 0; i < n_mubar; i++)
			maxbetabar = maxbetabar < beta_bar[i] ? beta_bar[i] : maxbetabar;
		return maxbetabar;

		/*
		 * double[] theta_a = new double [get_Q_a()*2]; for(int q=0;
		 * q<get_Q_a(); q++){ Complex theta_q_a = complex_eval_theta_q_a(q);
		 * theta_a[q] = theta_q_a.getReal(); theta_a[q+get_Q_a()] =
		 * theta_q_a.getImaginary(); } double betaUB = -1e8; for (int i = 0; i <
		 * n_mubar; i++){ double localbetaUB = 1e8; for (int j = 0; j <
		 * n_muhat[i]; j++){ double calbetaUB = 0.; for (int k = 0; k <
		 * get_Q_a()*2; k++) calbetaUB += zval[i][j][k]*theta_a[k]; localbetaUB
		 * = localbetaUB > calbetaUB ? calbetaUB : localbetaUB; } betaUB =
		 * betaUB < localbetaUB ? localbetaUB : betaUB; } return betaUB;
		 */

	}

	private List<Integer> getSorted_CJ_Indices(Vector<double[]> C_J) {

		int J = C_J.size();

		LinkedHashMap<Double, Integer> dist_from_mu = new LinkedHashMap<Double, Integer>(J);

		for (int j = 0; j < J; j++) {
			double dist = param_dist(get_current_parameters(), C_J.get(j));
			dist_from_mu.put(dist, j);
		}

		List<Map.Entry<Double, Integer>> list = new LinkedList<Map.Entry<Double, Integer>>(dist_from_mu.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<Double, Integer>>() {
			public int compare(Map.Entry<Double, Integer> o1, Map.Entry<Double, Integer> o2) {
				return o1.getKey().compareTo(o2.getKey());
				/*
				 * return ((Comparable<?>) ((Map.Entry<Double,Integer>)
				 * (o1)).getKey()) .compareTo(((Map.Entry<Double,Integer>)
				 * (o2)).getKey());
				 */
			}
		});

		// Create a sorted list of values to return
		List<Integer> result = new LinkedList<Integer>();
		for (Map.Entry<Double, Integer> entry : list) {
			result.add(entry.getValue());
		}

		return result;
	}

}
