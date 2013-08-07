package rb;

/**
 * Enum for known RB system types in JRB
 * 
 * This class has been taken from the original @ref rbappmit package and modified to fit into the current JaRMoS
 * framework, formerly RBEnums.SystemTypeEnum.
 * 
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 * 
 */
public enum SystemType {
	/**
	 * No system type. ("Default" value, used for unknown system types)
	 */
	NONE,

	/**
	 * Linear, time-invariant/steady state rb system
	 */
	LINEAR_STEADY,

	/**
	 * Linear, time-dependent rb system
	 */
	LINEAR_UNSTEADY,

	/**
	 * 
	 */
	QN_UNSTEADY,

	/**
	 * Linear, complex valued time-independent/steady state rb system
	 */
	LINEAR_COMPLEX_STEADY;

	/**
	 * Parses a string into the corresponding enum.
	 * 
	 * @param s
	 * The string to parse
	 * @return Returns the matching system type or SystemType.NONE if no matching type is found or s is null.
	 */
	public static SystemType parse(String s) {
		for (SystemType type : SystemType.values()) {
			if (type.toString().equals(s))
				return type;
		}
		return SystemType.NONE;
	}

	/**
	 * Instantiates an RBSystem subclass corresponding to the current type.
	 * 
	 * @return An RBSystem instance
	 */
	public RBSystem getNewRBSystem() {
		switch (this) {
		case NONE:
			return null;
		case LINEAR_STEADY:
			return new RBSystem();
		case LINEAR_COMPLEX_STEADY:
			return new ComplexRBSystem();
		case LINEAR_UNSTEADY:
			return new TransientRBSystem();
		case QN_UNSTEADY:
			return new QNTransientRBSystem();
		default:
			return null;
		}
	}
}
