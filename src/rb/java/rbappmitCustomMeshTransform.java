/**
 * 
 */
package rb.java;

import jarmos.geometry.MeshTransform;

/**
 * @author CreaByte
 *
 */
public class rbappmitCustomMeshTransform implements MeshTransform {
	
	private double[] mu;
	private RBSystem sys;
	
	public rbappmitCustomMeshTransform(double[] mu, RBSystem sys) {
		this.mu = mu.clone();
		this.sys = sys;
	}

	/** (non-Javadoc)
	 * @see jarmos.geometry.MeshTransform#transformMesh(float[])
	 */
	@Override
	public float[] transformMesh(float[] vertices) {
		return sys.rbappmitCustomMeshTransform(mu, vertices);
	}

}
