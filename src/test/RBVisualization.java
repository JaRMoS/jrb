/**
 * 
 */
package test;

import rb.java.RBContainer;
import rb.java.RBSystem;
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
public class RBVisualization {

	/**
	 * @param args
	 * @throws ModelManagerException 
	 */
	public static void main(String[] args) throws ModelManagerException {
		FileModelManager f = new FileModelManager("models");
		f.useModel("demo3");
		
		RBContainer rb = new RBContainer();
		rb.loadModel(f);
		
		// Perform the solve
		RBSystem s=rb.mRbSystem;
		double[] par = s.getParams().getRandomParam();
//		double[] par = new double[]{.5, .5};
		s.getParams().setCurrent(par);
		s.solveRB(s.getNBF()/2);
		
		SimulationResult res = s.getSimulationResults();
		GeometryData g = rb.mRbSystem.getGeometry();
		VisualizationData v = new VisualizationData(g);
		v.useResult(res);
		
		v.computeVisualFeatures(new ColorGenerator());
		
		//JOGLRenderer.render(v);
		
		s.performSweep(0, 4);
		res = s.getSweepSimResults();
		v.useResult(res);
		v.computeVisualFeatures(new ColorGenerator());
		
		JOGLRenderer.render(v);

	}

}
