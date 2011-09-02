/**
 * Created on Aug 29, 2011 in Project JRB
 * Location: test.TestGeometryLoading.java
 */
package test;

import static org.junit.Assert.fail;

import org.junit.Test;

import rmcommon.geometry.GeometryData;
import rmcommon.io.AModelManager.ModelManagerException;
import rmcommon.io.FileModelManager;

/**
 * @author Daniel Wirtz
 * @date Aug 29, 2011
 *
 */
public class TestGeometryLoading {

	/**
	 * Performs simple geometry loading for an JRB model
	 */
	@Test
	public void testGeometryLoading() {
		FileModelManager f = new FileModelManager("models");
		try {
			f.setModelDir("adv_diff_rb_tc");
		} catch (ModelManagerException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		GeometryData g = new GeometryData();
		g.loadModelGeometry(f);
	}

}
