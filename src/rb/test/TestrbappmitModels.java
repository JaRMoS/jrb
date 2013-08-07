package rb.test;

import static org.junit.Assert.assertTrue;
import jarmos.ModelDescriptor;
import jarmos.ModelType;
import jarmos.SimulationResult;
import jarmos.geometry.GeometryData;
import jarmos.io.AModelManager.ModelManagerException;
import jarmos.io.FileModelManager;
import jarmos.visual.ColorGenerator;
import jarmos.visual.VisualizationData;

import org.junit.Test;

import rb.RBContainer;
import rb.RBSystem;

/**
 * @author Daniel Wirtz
 * 
 */
public class TestrbappmitModels {
	/**
	 * Tests loading & solving of @ref rbappmit models (type "rbappmit")
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
