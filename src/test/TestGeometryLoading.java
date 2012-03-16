/**
 * Created on Aug 29, 2011 in Project JRB
 * Location: test.TestGeometryLoading.java
 */
package test;

import static org.junit.Assert.fail;
import jarmos.geometry.GeometryData;
import jarmos.io.AModelManager.ModelManagerException;
import jarmos.io.FileModelManager;

import org.junit.Test;


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
			f.useModel("adv_diff_rb_tc");
		} catch (ModelManagerException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		
		GeometryData g = new GeometryData();
		g.loadModelGeometry(f);
	}

}
