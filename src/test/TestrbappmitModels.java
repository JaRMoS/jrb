/**
 * 
 */
package test;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

import rb.java.RBContainer;
import rb.java.RBSystem;
import rmcommon.ModelDescriptor;
import rmcommon.ModelType;
import rmcommon.SimulationResult;
import rmcommon.geometry.GeometryData;
import rmcommon.io.AModelManager.ModelManagerException;
import rmcommon.io.FileModelManager;
import rmcommon.visual.ColorGenerator;
import rmcommon.visual.VisualizationData;

/**
 * @author CreaByte
 * 
 */
public class TestrbappmitModels {
	/**
	 * Tests loading & solving of rbappmit models (type "rbappmit")
	 * 
	 * @throws ModelManagerException
	 */
	@Test
	public void testRBSolve_rbappmit_Models() throws ModelManagerException {
		FileModelManager f = new FileModelManager("models");
		for (ModelDescriptor md : f.getModelDescriptors()) {
			if (md.type == ModelType.rbappmit) {
				f.useModel(md.modeldir);

				RBContainer rb = new RBContainer();
				assertTrue(rb.loadModel(f));

				// Perform the solve
				RBSystem s = rb.mRbSystem;
				double[] par = s.getParams().getRandomParam();
				s.getParams().setCurrent(par);
				s.computeRBSolution(2);

				SimulationResult sol = s.getSimulationResults();

				GeometryData g = new GeometryData();
				g.loadModelGeometry(f);

				VisualizationData d = new VisualizationData(g);
				d.useResult(sol);

				d.computeVisualFeatures(new ColorGenerator());
			}
		}
	}
}
