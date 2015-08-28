package supervised_bgp;

import weka.core.Instances;

public abstract class TransductiveClassifier {
	
	protected double[][] fUnlabeledDocs; 
	
	public abstract void buildClassifier(Instances dataTrain, Instances dataTest);

	
}
