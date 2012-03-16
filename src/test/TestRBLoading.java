/**
 * Created on Aug 23, 2011 in Project JRB
 * Location: test.TestRBLoading.java
 */
package test;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;
import jarmos.ModelDescriptor;
import jarmos.io.AModelManager;
import jarmos.io.AModelManager.ModelManagerException;
import jarmos.io.FileModelManager;

import org.junit.Test;

import rb.java.RBContainer;

/**
 * @author Daniel Wirtz
 * @date Aug 23, 2011
 * 
 */
public class TestRBLoading {

	/**
	 * Tests loading all the rb-models contained in the models folder (if there are any)
	 */
	@Test
	public void testModelLoading() {
		FileModelManager f = new FileModelManager("models");
		try {
			for (ModelDescriptor m : f.getModelDescriptors()) {
				if (("rb".equals(m.type) || "rbappmit".equals(m.type))) {
					System.out.println("\n----------------- Trying to load model "+m.title+" from folder "+m.modeldir+" -----------------\n");
					try {
						f.useModel(m.modeldir);
					} catch (ModelManagerException e) {
						e.printStackTrace();
						fail(e.getMessage());
					}

					assertTrue(f.modelFileExists(AModelManager.CLASSES_JARFILE));

					RBContainer rb = new RBContainer();
					assertTrue(rb.loadModel(f));
				}
			}
		} catch (ModelManagerException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

}
