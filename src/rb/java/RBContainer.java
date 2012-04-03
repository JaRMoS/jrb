package rb.java;

import jarmos.Log;
import jarmos.ModelType;
import jarmos.io.AModelManager;

import java.io.IOException;

import rb.java.affinefcn.IAffineFunctions;
import rb.java.affinefcn.rbappmitAffineFunctions;

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

	@SuppressWarnings("unchecked")
	private void loadAffineFunctions(AModelManager m) throws Exception {

		ClassLoader cl = m.getClassLoader();
		Class<?> af = cl.loadClass(m.getModelPackageStr() + "AffineFunctions");

		Log.d(DEBUG_TAG, "Loaded AffineFunctions class");

		if (mRbSystem != null) {
			mRbSystem.oldAffFcnCl = af;
			mRbSystem.oldAffFcnObj = af.newInstance();
			if (m.getModelType() == ModelType.rbappmit) {
				mRbSystem.affineFunctionsClass = IAffineFunctions.class;
				mRbSystem.affineFunctionsInstance = new rbappmitAffineFunctions(
						af, af.newInstance());
			} else {
				mRbSystem.affineFunctionsClass = (Class<IAffineFunctions>) af;
				mRbSystem.affineFunctionsInstance = mRbSystem.affineFunctionsClass
						.newInstance();
			}
		}

//		if (mRbScmSystem != null) {
//			mRbScmSystem.oldAffFcnCl = af;
//			mRbScmSystem.oldAffFcnObj = af.newInstance();
//			if (m.getModelType() == ModelType.rbappmit) {
//				mRbScmSystem.affineFunctionsClass = IAffineFunctions.class;
//				mRbScmSystem.affineFunctionsInstance = new rbappmitAffineFunctions(
//						af, af.newInstance());
//			} else {
//				mRbScmSystem.affineFunctionsClass = (Class<IAffineFunctions>) af;
//				mRbScmSystem.affineFunctionsInstance = mRbSystem.affineFunctionsClass
//						.newInstance();
//			}
//		}
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
		if (m.getModelType() == ModelType.rbappmit) {
			if (!readSystemDescriptionsRBAppMit(m))
				return false;
		} else {
			readSystemDescriptionsJRB(m);
		}

		// Init the main systems
		mRbSystem = fSystemType.getNewRBSystem();
		mRbScmSystem = fSCMType.getNewRBSCMSystem(mRbSystem);

		// Assign SCM system
		mRbSystem.setPrimarySCM(mRbScmSystem);

		try {
			/*
			 * Requires the system configurations to be set already for
			 * backwards compatibility! See rbappmitAffineFunctions class.
			 */
			loadAffineFunctions(m);
		} catch (Exception e) {
			Log.e(DEBUG_TAG,
					"Exception occurred while loading affine functions: "
							+ e.getMessage(), e);
			return false;
		}

		// Finally, load the RB and SCM system offline data
		try {
			if (mRbSystem != null) {
				// Read parameters into RB systems
				if (!mRbSystem.readConfiguration(m))
					return false;
				mRbSystem.loadOfflineData(m);
				Log.d(DEBUG_TAG, "Finished reading offline data for RBSystem.");
			}

			if (mRbScmSystem != null) {
				// Read parameters into SCM systems
				mRbScmSystem.readConfiguration(m);
				mRbScmSystem.loadOfflineData(m);
				Log.d(DEBUG_TAG,
						"Finished reading offline data for RBSCMSystem.");
			}

			/*
			 * A second SCM system seems not to be used within any of the
			 * current demos. If needed, uncomment this and one line in the
			 * loadSecondSCMSystem method to enable its use. Even so, it seems
			 * the second SCM system does not get loaded any offline data at the
			 * current state of the code.
			 */
			// loadSecondSCMSystem(m);

		} catch (Exception e) {
			Log.e(DEBUG_TAG, "Exception occurred while reading offline data: "
					+ e.getMessage(), e);
			return false;
		}
		return true;
	}

//	@SuppressWarnings("unused")
//	private boolean loadSecondSCMSystem(AModelManager m)
//			throws InstantiationException, IllegalAccessException {
//		mSecondRbScmSystem = null;
//		if (fSCMType == SCMType.COERCIVE_ALPHASIGMA) {
//			mSecondRbScmSystem = fSCMType.getNewRBSCMSystem();
//
//			if (mSecondRbScmSystem != null) {
//				if (fSystemType == SystemType.LINEAR_UNSTEADY) {
//					/*
//					 * Uncomment this line to restore previous status using an
//					 * optional second SCM system.
//					 */
//					// ((TransientRBSystem)
//					// mRbSystem).setSecondarySCM(mSecondRbScmSystem);
//				}
//
//				if (!mSecondRbScmSystem.readConfiguration(m))
//					return false;
//
//				// Attach AffineFunctions class also to this system
//				mSecondRbScmSystem.affineFunctionsClass = mRbScmSystem.affineFunctionsClass;
//				mSecondRbScmSystem.oldAffFcnCl = mRbSystem.oldAffFcnCl;
//				mSecondRbScmSystem.affineFunctionsInstance = mRbSystem.affineFunctionsClass
//						.newInstance();
//				mSecondRbScmSystem.oldAffFcnObj = mRbSystem.oldAffFcnObj;
//			}
//		}
//		return true;
//	}

	private void readSystemDescriptionsJRB(AModelManager m) {
		problemTitle = m.getModelXMLTagValue("description.name");

		// TODO: Dont know where this value is used
		descriptionURL = m.getModelXMLTagValue("description.infohtml");

		fSystemType = SystemType.parse(m
				.getModelXMLTagValue("rb_model.systype"));
		fSCMType = SCMType.parse(m.getModelXMLTagValue("rb_model.scmtype"));

		Log.d(DEBUG_TAG, "RB system type = " + fSystemType);
		Log.d(DEBUG_TAG, "SCM type = " + fSCMType);
	}

	private boolean readSystemDescriptionsRBAppMit(AModelManager m) {
		// GetPot infile = m.getParamFileGetPot();
		GetPot infile = null;
		try {
			infile = new GetPot(m.getInStream(Const.parameters_filename),
					Const.parameters_filename);
		} catch (IOException e) {
			Log.e("RBContainer", "Exception loading infile.in", e);
			e.printStackTrace();
			return false;
		}

		problemTitle = infile.call("title", "RB Online");
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