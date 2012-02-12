/**
 * Created on Aug 29, 2011 in Project JRB
 * Location: test.TestRBSolve.java
 */
package test;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import org.junit.Test;

import rb.java.RBContainer;
import rb.java.RBSystem;
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

//	/**
//	 * Tests loading & solving of JRB models (type "rb")
//	 */
//	@Test
//	public void testRBSolveJRB() {
//		FileModelManager f = new FileModelManager("models");
//		try {
//			f.setModelDir("adv_diff_rb");
//		} catch (ModelManagerException e) {
//			e.printStackTrace();
//			fail(e.getMessage());
//		}
//		
//		RBContainer rb = new RBContainer();
//		assertTrue(rb.loadModel(f));
//		
//		// Perform the solve
//		RBSystem s=rb.mRbSystem;
////		double[] par = s.getParams().getRandomParam();
//		double[] par = new double[]{.5, .5};
//		s.getParams().setCurrent(par);
//		s.RB_solve(4);
//	}
	
	/**
	 * Tests loading & solving of JRB models (type "rb")
	 */
	@Test
	public void testRBSolveJRB_TimeConst() {
		FileModelManager f = new FileModelManager("models");
		try {
			f.setModelDir("rbm_advec");
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
		s.RB_solve(4);
		
		float[][][] sol = s.getFullSolution();
		
		GeometryData g = new GeometryData();
		g.loadModelGeometry(f);
		VisualizationData d = new VisualizationData(g);
		d.set1FieldData(sol[0][0]);
		d.computeColorData(new ColorGenerator());
	}
	
//	/**
//	 * Tests loading & solving of JRB models (type "rb")
//	 */
//	@Test
//	public void testRBSolverbappmit() {
//		FileModelManager f = new FileModelManager("models");
//		try {
//			f.setModelDir("demo8");
//		} catch (ModelManagerException e) {
//			e.printStackTrace();
//			fail(e.getMessage());
//		}
//		
//		RBContainer rb = new RBContainer();
//		assertTrue(rb.loadModel(f));
//		
//		// Perform the solve
//		RBSystem s=rb.mRbSystem;
//		double[] par = s.getParams().getRandomParam();
////		double[] par = new double[]{.5, .5};
//		s.getParams().setCurrent(par);
//		s.RB_solve(s.getNBF()/2);
//	}

}
