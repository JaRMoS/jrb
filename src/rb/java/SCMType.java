package rb.java;

/**
 * Enum containing the known SCM (successive constraint method) model types in
 * JRB.
 * 
 * Values taken from the previous rbappmit RBEnums.SCMTypeEnum class.
 * 
 * @author Daniel Wirtz
 * @date Aug 28, 2011
 * 
 */
public enum SCMType {
	/**
	 * Not an SCM system
	 */
	NONE,

	/**
	 * A coervice system
	 */
	COERCIVE,

	/**
	 * A coercive system with alpha/sigma
	 */
	COERCIVE_ALPHASIGMA,

	/**
	 * Dunno :-)
	 */
	QN_TRANSIENT_SCM,

	/**
	 * A complex-valued non-coercive system.
	 */
	COMPLEX_NONCOERCIVE;

	/**
	 * @param s
	 *            The string to parse
	 * @return The matching SCM type or SCMType.NONE if no match was found or s
	 *         is null.
	 */
	public static SCMType parse(String s) {
		for (SCMType t : SCMType.values()) {
			if (t.toString().equals(s)) return t;
		}
		return SCMType.NONE;
	}

	/**
	 * Creates a new instance of an RBSCMSystem subclass according to the current SCMType.
	 * @return A new RBSCMSystem instance
	 */
	public RBSCMSystem getNewRBSCMSystem(RBSystem sys) {
		switch (this) {
		case NONE:
			return null;
		case COERCIVE:
			return new RBSCMSystem(sys);
		case COERCIVE_ALPHASIGMA:
			return new RBSCMSystem(sys);
		case QN_TRANSIENT_SCM:
			return new QNTransientSCMSystem(sys);
		case COMPLEX_NONCOERCIVE:
			return new RBnSCMCSystem(sys);
		default:
			return null;
		}
	}
}
