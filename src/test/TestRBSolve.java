/**
 * Created on Aug 29, 2011 in Project JRB
 * Location: test.TestRBSolve.java
 */
package test;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;

import org.junit.Test;

import rb.java.RBContainer;
import rb.java.RBSystem;
import rmcommon.SimulationResult;
import rmcommon.geometry.GeometryData;
import rmcommon.io.AModelManager.ModelManagerException;
import rmcommon.io.FileModelManager;
import rmcommon.visual.ColorGenerator;
import rmcommon.visual.VisualizationData;

/**
 * @author Daniel Wirtz
 * @date Aug 29, 2011
 *
 */
public class TestRBSolve {

	/**
	 * Tests loading & solving of JRB models (type "rb")
	 */
	//@Test
	public void testRBSolveJRB() {
		FileModelManager f = new FileModelManager("models");
		try {
			f.useModel("rbm_advec");
		} catch (ModelManagerException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		RBContainer rb = new RBContainer();
		assertTrue(rb.loadModel(f));
		
		// Perform the solve
		RBSystem s=rb.mRbSystem;
//		double[] par = s.getParams().getRandomParam();
		double[] par = new double[]{.5, .5};
		s.getParams().setCurrent(par);
		s.solveRB(4);
	}
	
	/**
	 * Tests loading & solving of JRB models (type "rb")
	 */
	//@Test
	public void testRBSolveJRB_TimeConst() {
		FileModelManager f = new FileModelManager("models");
		try {
			f.useModel("rbm_advec");
		} catch (ModelManagerException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		RBContainer rb = new RBContainer();
		assertTrue(rb.loadModel(f));
		
		// Perform the solve
		RBSystem s=rb.mRbSystem;
//		double[] par = s.getParams().getRandomParam();
		double[] par = new double[]{.5, .5, .5};
		s.getParams().setCurrent(par);
		s.solveRB(4);
		
		SimulationResult sol = s.getSimulationResults();
		
		GeometryData g = new GeometryData();
		g.loadModelGeometry(f);
		VisualizationData d = new VisualizationData(g);
		d.useResult(sol);
		d.computeVisualFeatures(new ColorGenerator());
	}
	
	
	
	/**
	 * Tests loading & solving of JRB models (type "rb")
	 */
	@Test
	public void testRBSolverbappmit() {
		FileModelManager f = new FileModelManager("models");
		try {
			f.useModel("demo3");
		} catch (ModelManagerException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		RBContainer rb = new RBContainer();
		assertTrue(rb.loadModel(f));
		
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
		
		s.performSweep(0, 4);
		res = s.getSweepSimResults();
		v.useResult(res);
		v.computeVisualFeatures(new ColorGenerator());
	}

}
