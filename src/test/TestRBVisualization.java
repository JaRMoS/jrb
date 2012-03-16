/**
 * 
 */
package test;

import jarmos.SimulationResult;
import jarmos.geometry.GeometryData;
import jarmos.io.AModelManager.ModelManagerException;
import jarmos.io.FileModelManager;
import jarmos.visual.ColorGenerator;
import jarmos.visual.JOGLRenderer;
import jarmos.visual.VisualizationData;
import rb.java.RBContainer;
import rb.java.RBSystem;

/**
 * @author CreaByte
 *
 */
public class TestRBVisualization {

	/**
	 * @param args
	 * @throws ModelManagerException 
	 */
	public static void main(String[] args) throws ModelManagerException {
		FileModelManager f = new FileModelManager("models");
		f.useModel("demo8");
		
		RBContainer rb = new RBContainer();
		rb.loadModel(f);
		
		// Perform the solve
		RBSystem s=rb.mRbSystem;
		double[] par = s.getParams().getRandomParam();
		s.getParams().setCurrent(par);
		s.computeRBSolution(s.getNBF());
		SimulationResult res = s.getSimulationResults();
		
//		s.performSweep(0, 4);
//		SimulationResult res = s.getSweepSimResults();
		
		GeometryData g = rb.mRbSystem.getGeometry();
		VisualizationData v = new VisualizationData(g);
		
		v.useResult(res);
		v.computeVisualFeatures(new ColorGenerator());
				
		JOGLRenderer.render(v);
	}

}
