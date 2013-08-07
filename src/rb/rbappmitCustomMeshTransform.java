package rb;

import jarmos.geometry.MeshTransform;

/**
 * Compatibility class for old @ref rbappmit models with custom mesh/geometry transformations.
 * 
 * @author Daniel Wirtz
 * 
 */
public class rbappmitCustomMeshTransform implements MeshTransform {

	private double[] mu;
	private RBSystem sys;

	public rbappmitCustomMeshTransform(double[] mu, RBSystem sys) {
		this.mu = mu.clone();
		this.sys = sys;
	}

	/**
	 * (non-Javadoc)
	 * 
	 * @see jarmos.geometry.MeshTransform#transformMesh(float[])
	 */
	@Override
	public float[] transformMesh(float[] vertices) {
		return sys.rbappmitCustomMeshTransform(mu, vertices);
	}

}
