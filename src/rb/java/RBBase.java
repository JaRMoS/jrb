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

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.linear.ArrayFieldVector;
import org.apache.commons.math.linear.FieldVector;

import rb.java.affinefcn.IAffineFunctions;
import rmcommon.Log;
import rmcommon.ModelType;
import rmcommon.Parameters;
import rmcommon.io.AModelManager;

// 
// This class is modeled on the RBBase class in rbOOmit

/**
 * Base class for RB and SCM systems, stores the current parameter value,
 * parameter ranges, as well as other data common to RBSystems and SCMSystems.
 * 
 * Changes made by
 * 
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 * 
 */
public abstract class RBBase {

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
	
	/**
	 * @deprecated Here for rbappmit compatibility.
	 */
	public Class<?> oldAffFcnCl;
	/**
	 * @deprecated Here for rbappmit compatibility.
	 */
	public Object oldAffFcnObj;

	/**
	 * The number of terms in the affine expansion of the bilinear form
	 */
	private int fQa;

	/**
	 * 
	 */
	public boolean isReal = true;

	/**
	 * The system's parameters object containing the parameter values and
	 * descriptions
	 */
	private Parameters params;

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
	 * Evaluate theta_q_a (for the q^th bilinear form) at the current parameter.
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
	 * @return
	 */
	public boolean is_custom_mesh_transform() {
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

	/**
	 * @param mu
	 * @param x
	 * @return
	 */
	public float[] mesh_transform(double[] mu, float[] x) {
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
	}

	protected void readConfigurationRBAppMIT(GetPot infile) {
		int n_parameters = infile.call("n_parameters", 1);
		Log.d("RBBase", "n_parameters = " + n_parameters);
		if (n_parameters > 0) {
			params = new Parameters();
			for (int i = 0; i < n_parameters; i++) {
				String label = infile.call("param" + Integer.toString(i)
						+ "_label", "mu_" + i);
				// Read in the min/max for the i^th parameter
				String min_string = new String("mu" + i + "_min");
				double mu_i_min = infile.call(min_string, 0.);
				String max_string = new String("mu" + i + "_max");
				double mu_i_max = infile.call(max_string, 0.);
				params.addParam(label, mu_i_min, mu_i_max);
				Log.d("RBBase", "Parameter " + i + ": Min = " + mu_i_min
						+ ", Max = " + mu_i_max);
			}
		}
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

}
