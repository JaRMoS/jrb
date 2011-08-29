/**
 * Created on Aug 29, 2011 in Project JRB
 * Location: rb.java.affinefcn.rbappmitAffineFunctions.java
 */
package rb.java.affinefcn;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * @author Daniel Wirtz
 * @date Aug 29, 2011
 *
 */
public class rbappmitAffineFunctions implements IAffineFunctions, ITransient,
		IWithuL {
	
	/*
	 * The wrapped classes.
	 */
	private Class<?> affineFunctionsClass;
	private Object affineFunctionsInstance;
	
	/**
	 * Creates a wrapper for the old rbappmit affine functions classes.
	 * @param theclass
	 * @param inst
	 */
	public rbappmitAffineFunctions(Class<?> theclass, Object inst) {
		affineFunctionsClass = theclass;
		affineFunctionsInstance = inst;
	}

	/**
	 * @see rb.java.affinefcn.IWithuL#getQuL()
	 */
	@Override
	public int getQuL() {
		/*
		 * Extra for QuL, as it checks via exception if function is there
		 */
		Method meth;
		boolean noQ_uLdefined = false;
		try {
			// Get a reference to get_n_A_functions, which does not
			// take any arguments
			meth = affineFunctionsClass.getMethod("get_n_uL_functions", (Class<?>[]) null);
		} catch (NoSuchMethodException nsme) {
			// throw new
			// RuntimeException("getMethod for get_n_uL_functions failed",
			// nsme);
			noQ_uLdefined = true;
			meth = null;
		}

		if (noQ_uLdefined)
			return 0;
		else {
			Integer Q_uL;
			try {
				Object Q_uL_obj = meth.invoke(affineFunctionsInstance, (Object[]) null);
				Q_uL = (Integer) Q_uL_obj;
			} catch (IllegalAccessException iae) {
				throw new RuntimeException(iae);
			} catch (InvocationTargetException ite) {
				throw new RuntimeException(ite.getCause());
			}
			return Q_uL.intValue();
		}
	}

	/**
	 * 
	 * @see rb.java.affinefcn.ITransient#thetaQm(int, double[], double)
	 */
	@Override
	public double thetaQm(int i, double[] p, double t) {
		Method meth;

		try {
			// Get a reference to get_n_M_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[2];
			partypes[0] = Integer.TYPE;
			partypes[1] = double[].class;

			meth = affineFunctionsClass.getMethod("evaluateM", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateM failed", nsme);
		}

		Double theta_val;
		try {
			Object arglist[] = new Object[2];
			arglist[0] = new Integer(i);
			arglist[1] = p;

			Object theta_obj = meth.invoke(affineFunctionsInstance, arglist);
			theta_val = (Double) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return theta_val.doubleValue();
	}

	/**
	 * @see rb.java.affinefcn.ITransient#getQm()
	 */
	@Override
	public int getQm() {
		Method meth;
		try {
			// Get a reference to get_n_A_functions, which does not
			// take any arguments
			meth = affineFunctionsClass.getMethod("get_n_M_functions", (Class<?>[]) null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_n_M_functions failed", nsme);
		}
		try {
			Object Q_m_obj = meth.invoke(affineFunctionsInstance, (Object[]) null);
			return (Integer)Q_m_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#isTimeDependentL()
	 */
	@Override
	public boolean isTimeDependentL() {
		// TODO Auto-generated method stub
		return false;
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#isTimeDependentA()
	 */
	@Override
	public boolean isTimeDependentA() {
		// TODO Auto-generated method stub
		return false;
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#getQl()
	 */
	@Override
	public int[] getQl() {
		Method meth;
		try {
			Class<?> partypes[] = new Class[1];
			partypes[0] = Integer.TYPE;
			meth = affineFunctionsClass.getMethod("get_Q_l", partypes);
			int no = getNumOutputs();
			int[] Ql_values = new int[no];
			Object arglist[] = new Object[1];
			for (int i = 0; i < no; i++) {
				arglist[0] = i;
				Ql_values[i] = (Integer) meth.invoke(affineFunctionsInstance, arglist);
			}
			return Ql_values;
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod failed: " + nsme.getMessage(), nsme);
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#thetaQl(int, int, double[], double)
	 */
	@Override
	public double thetaQl(int i, int q_l, double[] p, double t) {
		Method meth;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[3];
			partypes[0] = Integer.TYPE;
			partypes[1] = Integer.TYPE;
			partypes[2] = double[].class;

			meth = affineFunctionsClass.getMethod("evaluateL", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateL failed", nsme);
		}

		Double theta_val;
		try {
			Object arglist[] = new Object[3];
			arglist[0] = new Integer(i);
			arglist[1] = new Integer(q_l);
			arglist[2] = p;

			Object theta_obj = meth.invoke(affineFunctionsInstance, arglist);
			theta_val = (Double) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return theta_val.doubleValue();
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#getQf()
	 */
	@Override
	public int getQf() {
		try {
			Method meth = affineFunctionsClass.getMethod("get_n_F_functions", (Class<?>[]) null);
			return (Integer) meth.invoke(affineFunctionsInstance, (Object[]) null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod failed: " + nsme.getMessage(), nsme);
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#thetaQf(int, double[], double)
	 */
	@Override
	public double thetaQf(int i, double[] p, double t) {
		Method meth;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[2];
			partypes[0] = Integer.TYPE;
			partypes[1] = double[].class;

			meth = affineFunctionsClass.getMethod("evaluateF", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateF failed", nsme);
		}

		Double theta_val;
		try {
			Object arglist[] = new Object[2];
			arglist[0] = i;
			arglist[1] = p;

			Object theta_obj = meth.invoke(affineFunctionsInstance, arglist);
			theta_val = (Double) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return theta_val.doubleValue();
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#getQa()
	 */
	@Override
	public int getQa() {
		Method meth;

		try {
			// Get a reference to get_n_A_functions, which does not
			// take any arguments
			meth = affineFunctionsClass.getMethod("get_n_A_functions", (Class<?>[]) null);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for get_n_A_functions failed", nsme);
		}

		Integer Q_a;
		try {
			Object Q_a_obj = meth.invoke(affineFunctionsInstance, (Object[]) null);
			Q_a = (Integer) Q_a_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return Q_a.intValue();
	}

	/**
	 * @see rb.java.affinefcn.IAffineFunctions#thetaQa(int, double[], double)
	 */
	@Override
	public double thetaQa(int i, double[] p, double t) {
		Method meth;

		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments

			Class<?> partypes[] = new Class[2];
			partypes[0] = Integer.TYPE;
			partypes[1] = double[].class;

			meth = affineFunctionsClass.getMethod("evaluateA", partypes);
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod for evaluateA failed", nsme);
		}

		Double theta_val;
		try {
			Object arglist[] = new Object[2];
			arglist[0] = i;
			arglist[1] = p;

			Object theta_obj = meth.invoke(affineFunctionsInstance, arglist);
			theta_val = (Double) theta_obj;
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}

		return theta_val.doubleValue();
	}

	
	/**
	 * @see rb.java.affinefcn.IAffineFunctions#getNumOutputs()
	 */
	@Override
	public int getNumOutputs() {
		try {
			// Get a reference to get_n_L_functions, which does not
			// take any arguments
			Method meth = affineFunctionsClass.getMethod("get_n_outputs", (Class<?>[]) null);
			Object n_outputs_obj = meth.invoke(affineFunctionsInstance, (Object[]) null);
			Integer n_outputs = (Integer) n_outputs_obj;
			return n_outputs.intValue();
		} catch (NoSuchMethodException nsme) {
			throw new RuntimeException("getMethod failed: " + nsme.getMessage(), nsme);
		} catch (IllegalAccessException iae) {
			throw new RuntimeException(iae);
		} catch (InvocationTargetException ite) {
			throw new RuntimeException(ite.getCause());
		}
	}

	/**
	 * @see rb.java.affinefcn.ITransient#isTimeDependentM()
	 */
	@Override
	public boolean isTimeDependentM() {
		return false;
	}

}
