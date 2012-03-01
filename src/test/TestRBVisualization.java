/**
 * 
 */
package test;

import rb.java.RBContainer;
import rb.java.RBSystem;
import rb.java.TransientRBSystem;
import rmcommon.SimulationResult;
import rmcommon.geometry.GeometryData;
import rmcommon.io.AModelManager.ModelManagerException;
import rmcommon.io.FileModelManager;
import rmcommon.visual.ColorGenerator;
import rmcommon.visual.JOGLRenderer;
import rmcommon.visual.VisualizationData;

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
		s.solveRB(s.getNBF());
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
