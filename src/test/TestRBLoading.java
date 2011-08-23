/**
 * Created on Aug 23, 2011 in Project JRB
 * Location: test.TestRBLoading.java
 */
package test;

import static org.junit.Assert.*;

import org.junit.Test;

import rb.java.RBContainer;
import rmcommon.io.AModelManager.ModelManagerException;
import rmcommon.io.FileModelManager;

/**
 * @author Daniel Wirtz
 * @date Aug 23, 2011
 * 
 */
public class TestRBLoading {

	/**
	 * 
	 */
	@Test
	public void test() {

		FileModelManager f = new FileModelManager("models");
		try {
			f.setModelDir("aghdemo");
		}
		catch (ModelManagerException e) {
			fail(e.getMessage());
		}

		RBContainer rb = new RBContainer();
		rb.loadModel(f);
	}

}
