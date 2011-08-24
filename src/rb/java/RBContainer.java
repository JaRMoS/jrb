package rb.java;

import java.io.IOException;

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
	 * The System type
	 */
	public RBEnums.SystemTypeEnum mSystemType;

	/**
	 * The SCM type
	 */
	public RBEnums.SCMTypeEnum mSCMType;

	/**
	 * Descriptive member variables for the problem title, variable names, and
	 * general info.
	 */
	public String problemTitle;
	/**
	 * Labels for the parameters
	 */
	public String[] paramLabels;
	/**
	 * A short problem description. Though not noticed where it ist used, yet.
	 */
	public String problemDescription;
	/**
	 * Some url :-)
	 */
	public String descriptionURL;

	/**
	 * The RBSystem object
	 */
	public RBSystem mRbSystem;
	/**
	 * The main RBSCMSystem object
	 */
	public RBSCMSystem mRbScmSystem;
	/**
	 * The second RBSCMSystem object, needed in some time-dependent problems
	 */
	public RBSCMSystem mSecondRbScmSystem;

	/**
	 * Creates a new RB Container with null values and NONE enum types.
	 */
	public RBContainer() {
		// Initialize the RB system and SCM types to NONE
		mSystemType = RBEnums.SystemTypeEnum.NONE;
		mSCMType = RBEnums.SCMTypeEnum.NONE;

		// Set rb.mRbScmSystem and rb.mRbSystem to null initially
		mRbSystem = null;
		mRbScmSystem = null;

		// Set the secondary SCM system to null also
		mSecondRbScmSystem = null;
	}

	private static RBEnums.SystemTypeEnum getSystemEnumFromString(String s) {

		for (RBEnums.SystemTypeEnum type : RBEnums.SystemTypeEnum.values()) {
			if (type.toString().equals(s)) {
				return type;
			}
		}

		return RBEnums.SystemTypeEnum.NONE;
	}

	private void read_system_types_from_input_file(AModelManager m)
			throws IOException {

		// GetPot infile = m.getParamFileGetPot();
		GetPot infile = new GetPot(m.getInStream(Const.parameters_filename), Const.parameters_filename);

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
		mSystemType = getSystemEnumFromString(SystemTypeEnum_in);

		String SCMTypeEnum_in = infile.call("scm_type", "NONE");
		mSCMType = getSCMEnumFromString(SCMTypeEnum_in);

		Log.d(DEBUG_TAG, "RB system type = " + mSystemType);
		Log.d(DEBUG_TAG, "SCM type = " + mSCMType);
	}

	/**
	 * Loads an rb/rbappmit type model using a provided ModelManager.
	 * 
	 * @param m
	 *            The ModelManager instance to use loading the model
	 * @return true if succeeded, false otherwise
	 */
	public boolean loadModel(AModelManager m) {

		boolean success = true;

		// Note: the local copy of the parameters_filename that has been
		// created here
		// if the source was a remote website has moved to
		// ModelManager.getParamFileGetPot()

		try {
			initialize_systems(m);
		} catch (InconsistentStateException e) {
			Log.e(DEBUG_TAG, "Inconsistent state exception occurred when parsing input file: "
					+ e.getMessage());
			success = false;
		} catch (IOException e) {
			Log.e(DEBUG_TAG, "I/O Exception thrown when accessing "
					+ Const.parameters_filename + ": " + e.getMessage());
			success = false;
		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception thrown when accessing "
					+ Const.parameters_filename, e);
			success = false;
		}

		try {
			attach_affine_functions(m);
		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception occurred while attaching affine functions: "
					+ e.getMessage(), e);
			success = false;
		}

		// Finally, initialize the RB and SCM systems
		try {
			if (mRbSystem != null) {
				mRbSystem.read_offline_data(m);
				Log.d(DEBUG_TAG, "Finished reading offline data for RBSystem.");
			}

			if (mRbScmSystem != null) {
				mRbScmSystem.read_offline_data(m);
				Log.d(DEBUG_TAG, "Finished reading offline data for RBSCMSystem.");
			}

		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception occurred while reading offline data: "
					+ e.getMessage(), e);
			success = false;
		}
		return success;
	}

	private void initialize_systems(AModelManager m) throws Exception {

		// First, clear all the systems in case we're doing a new (different)
		// problem
		mRbScmSystem = null;
		mSecondRbScmSystem = null;
		mRbSystem = null;

		// Find out which type of systems we'll be initializing
		read_system_types_from_input_file(m);

		// Initialize the SCM systems
		if (mSCMType == RBEnums.SCMTypeEnum.COERCIVE_ALPHASIGMA) {
			mRbScmSystem = RBSCMSystem.buildSCMSystem(RBEnums.SCMTypeEnum.COERCIVE);
			mSecondRbScmSystem = RBSCMSystem.buildSCMSystem(RBEnums.SCMTypeEnum.COERCIVE);
		} else if (mSCMType == RBEnums.SCMTypeEnum.COERCIVE
				|| mSCMType == RBEnums.SCMTypeEnum.QN_TRANSIENT_SCM) {
			mRbScmSystem = RBSCMSystem.buildSCMSystem(mSCMType);
		} else if (mSCMType == RBEnums.SCMTypeEnum.COMPLEX_NONCOERCIVE) {
			mRbScmSystem = RBSCMSystem.buildSCMSystem(mSCMType);
		} else {
			mRbScmSystem = null;
		}

		// Read parameters into SCM systems
		if (mRbScmSystem != null) {
			mRbScmSystem.parse_parameters_file(m);
		}
		if (mSecondRbScmSystem != null) {
			mSecondRbScmSystem.parse_parameters_file(m);
		}

		mRbSystem = RBSystem.buildRBSystem(mSystemType);
		mRbSystem.setPrimarySCM(mRbScmSystem);
		if (mSystemType == RBEnums.SystemTypeEnum.LINEAR_UNSTEADY) {
			TransientRBSystem trans_rb = (TransientRBSystem) mRbSystem;
			trans_rb.setSecondarySCM(mSecondRbScmSystem);
		}

		// Read parameters into RB systems
		if (mRbSystem != null) {
			mRbSystem.parse_parameters_file(m);
			/*
			 * if(mRbSystem.get_mfield() > 0) mRbModel = new
			 * GLObject(RBActivity.this);
			 */
		}
	}

	private void attach_affine_functions(AModelManager m) throws Exception {

		ClassLoader	cl = m.getClassLoader();	
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

			if (mSystemType == RBEnums.SystemTypeEnum.LINEAR_UNSTEADY
					|| mSystemType == RBEnums.SystemTypeEnum.QN_UNSTEADY) {
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
		if (mSecondRbScmSystem != null) {
			mSecondRbScmSystem.mAffineFnsClass = af;
			mSecondRbScmSystem.mTheta = mRbSystem.mAffineFnsClass.newInstance();

			// set Q_a
			mSecondRbScmSystem.read_in_Q_a();
		}

	}

	private static RBEnums.SCMTypeEnum getSCMEnumFromString(String s) {

		for (RBEnums.SCMTypeEnum type : RBEnums.SCMTypeEnum.values()) {
			if (type.toString().equals(s)) return type;
		}

		return RBEnums.SCMTypeEnum.NONE;
	}

}