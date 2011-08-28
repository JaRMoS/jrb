package rb.java;

import java.io.IOException;

import rmcommon.Log;
import rmcommon.Parameters;
import rmcommon.io.AModelManager;

/**
 * Base class for RB models and systems, extracted from the old RBActivity class
 * of rbappmit in order to separate data from the activity.
 * 
 * @author Daniel Wirtz
 * @date Aug 23, 2011
 * 
 */
public class RBContainer {

	// String for log printing
	static final String DEBUG_TAG = "RBContainer";

	/**
	 * Some url :-)
	 */
	public String descriptionURL;
	private SCMType fSCMType = SCMType.NONE;

	private SystemType fSystemType = SystemType.NONE;
	/**
	 * The main RBSCMSystem object
	 */
	public RBSCMSystem mRbScmSystem = null;
	/**
	 * The RBSystem object
	 */
	public RBSystem mRbSystem = null;

	/**
	 * The second RBSCMSystem object, needed in some time-dependent problems
	 */
	public RBSCMSystem mSecondRbScmSystem = null;

	/**
	 * Labels for the parameters
	 */
	public String[] paramLabels;
	/**
	 * A short problem description. Though not noticed where it ist used, yet.
	 */
	public String problemDescription;
	/**
	 * Descriptive member variables for the problem title, variable names, and
	 * general info.
	 */
	public String problemTitle;

	/**
	 * @return the SCMType
	 */
	public SCMType getSCMType() {
		return fSCMType;
	}

	/**
	 * @return the SystemType
	 */
	public SystemType getSystemType() {
		return fSystemType;
	}

	private void loadAffineFunctions(AModelManager m) throws Exception {

		ClassLoader cl = m.getClassLoader();
		Class<?> af = cl.loadClass("AffineFunctions");

		Log.d(DEBUG_TAG, "Loaded AffineFunctions class");

		if (mRbSystem != null) {
			mRbSystem.mAffineFnsClass = af;
			mRbSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();

			// Set Q_a, Q_f and n_outputs from the loaded class
			mRbSystem.read_in_Q_a();
			Log.d(DEBUG_TAG, "Q_a = " + mRbSystem.get_Q_a());

			mRbSystem.read_in_Q_f();
			Log.d(DEBUG_TAG, "Q_f = " + mRbSystem.get_Q_f());

			mRbSystem.read_in_n_outputs();
			Log.d(DEBUG_TAG, "n_outputs = " + mRbSystem.get_n_outputs());

			mRbSystem.read_in_Q_uL();

			if (fSystemType == SystemType.LINEAR_UNSTEADY
					|| fSystemType == SystemType.QN_UNSTEADY) {
				TransientRBSystem trans_rb = (TransientRBSystem) mRbSystem;
				trans_rb.read_in_Q_m();
				Log.d(DEBUG_TAG, "Q_m = " + trans_rb.get_Q_m());
			}
		}

		if (mRbScmSystem != null) {
			mRbScmSystem.mAffineFnsClass = af;
			mRbScmSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();

			// set Q_a
			mRbScmSystem.read_in_Q_a();
		}

	}

	/**
	 * Loads an rb/rbappmit type model using a provided ModelManager.
	 * 
	 * @param m
	 *            The ModelManager instance to use loading the model
	 * @return true if succeeded, false otherwise
	 */
	public boolean loadModel(AModelManager m) {

		// First, clear all the systems in case we're doing a new (different)
		// problem
		mRbScmSystem = null;
		mRbSystem = null;

		// Read system types and misc data into RBContainer
		if ("rbappmit".equals(m.getModelType())) {
			if (!readSystemDescriptionsRBAppMit(m)) return false;
		} else {
			readSystemDescriptionsJRB(m);
		}

		// Init the main systems
		mRbScmSystem = fSCMType.getNewRBSCMSystem();
		mRbSystem = fSystemType.getNewRBSystem();
		
		// Assign SCM system
		mRbSystem.setPrimarySCM(mRbScmSystem);

		try {
			loadAffineFunctions(m);
		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception occurred while attaching affine functions: "
					+ e.getMessage(), e);
			return false;
		}

		// Finally, initialize the RB and SCM systems
		try {
			if (mRbSystem != null) {
				// Read parameters into RB systems
				if (!mRbSystem.readConfiguration(m)) return false;
				// Load offline data
				mRbSystem.loadOfflineData(m);
				Log.d(DEBUG_TAG, "Finished reading offline data for RBSystem.");
			}

			if (mRbScmSystem != null) {
				// Read parameters into SCM systems
				if (!mRbScmSystem.readConfiguration(m)) return false;
				mRbScmSystem.loadOfflineData(m);
				Log.d(DEBUG_TAG, "Finished reading offline data for RBSCMSystem.");
			}

			/*
			 * A second SCM system seems not to be used within any of the current demos.
			 * If needed, uncomment this and one line in the loadSecondSCMSystem method to enable
			 * its use. Even so, it seems the second SCM system does not get loaded any offline data
			 * at the current state of the code.
			 */
//			loadSecondSCMSystem(m);
			
		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception occurred while reading offline data: "
					+ e.getMessage(), e);
			return false;
		}
		return true;
	}

	@SuppressWarnings("unused")
	private boolean loadSecondSCMSystem(AModelManager m) throws InstantiationException, IllegalAccessException {
		mSecondRbScmSystem = null;
		if (fSCMType == SCMType.COERCIVE_ALPHASIGMA) {
			mSecondRbScmSystem = fSCMType.getNewRBSCMSystem();

			if (mSecondRbScmSystem != null) {
				if (fSystemType == SystemType.LINEAR_UNSTEADY) {
					/*
					 * Uncomment this line to restore previous status using an optional second SCM system.
					 */
//					((TransientRBSystem) mRbSystem).setSecondarySCM(mSecondRbScmSystem);
				}

				if (!mSecondRbScmSystem.readConfiguration(m)) return false;

				// Attach AffineFunctions class also to this system
				mSecondRbScmSystem.mAffineFnsClass = mRbScmSystem.mAffineFnsClass;
				mSecondRbScmSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();

				// set Q_a
				mSecondRbScmSystem.read_in_Q_a();
			}
		}
		return true;
	}

	private void readSystemDescriptionsJRB(AModelManager m) {
		problemTitle = m.getModelXMLTagValue("description.name");

		Parameters p = m.getParameters();
		int params = p.getParamNumber();
		// TODO: remove paramLabels field and access Parameters object for names
		paramLabels = new String[params];
		for (int i = 0; i < params; i++) {
			paramLabels[i] = p.getParams().get(i).label;
		}

		// TODO: Dont know where this value is used
		descriptionURL = m.getModelXMLTagValue("description.infohtml");

		fSystemType = SystemType.parse(m.getModelXMLTagValue("rb_model.systype"));
		fSCMType = SCMType.parse(m.getModelXMLTagValue("rb_model.scmtype"));

		Log.d(DEBUG_TAG, "RB system type = " + fSystemType);
		Log.d(DEBUG_TAG, "SCM type = " + fSCMType);
	}

	private boolean readSystemDescriptionsRBAppMit(AModelManager m) {
		// GetPot infile = m.getParamFileGetPot();
		GetPot infile = null;
		try {
			infile = new GetPot(m.getInStream(Const.parameters_filename), Const.parameters_filename);
		} catch (IOException e) {
			Log.e("RBContainer", "Exception loading infile.in", e);
			e.printStackTrace();
			return false;
		}

		problemTitle = infile.call("title", "RB Online");
		int parameter_number = infile.call("n_parameters", 0);
		paramLabels = new String[parameter_number];
		for (int n = 0; n < parameter_number; n++) {
			paramLabels[n] = infile.call("param" + Integer.toString(n)
					+ "_label", "DEFAULT");

			// If DEFAULT was read in, then replace with a default
			// mu label
			if (paramLabels[n] == "DEFAULT") {
				paramLabels[n] = "\u00B5" + (n + 1);
			}
		}

		descriptionURL = infile.call("descriptionURL", "");

		String SystemTypeEnum_in = infile.call("system_type", "NONE");
		fSystemType = SystemType.parse(SystemTypeEnum_in);

		String SCMTypeEnum_in = infile.call("scm_type", "NONE");
		fSCMType = SCMType.parse(SCMTypeEnum_in);

		Log.d(DEBUG_TAG, "RB system type = " + fSystemType);
		Log.d(DEBUG_TAG, "SCM type = " + fSCMType);
		return true;
	}

}