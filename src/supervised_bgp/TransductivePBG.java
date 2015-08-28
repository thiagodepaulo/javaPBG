package supervised_bgp;

import java.util.Arrays;
import java.util.HashMap;

import weka.core.Instances;
import bgp.BipartiteGraph;
import bgp.Matrix;
import bgp.PBG;

public class TransductivePBG extends TransductiveClassifier {

	private int numTrain; // Number of labeled documents
	private int numTest; // Number of unlabeled documents
	private int numClasses; // Number of classes
	private int numTerms; // Number of terms
	
	private double alpha;
	private double beta;

	private PBG pbg;	

	public void buildClassifier(Instances dataTrain, Instances dataTest) {
		this.numTrain = dataTrain.numInstances();
		this.numTest = dataTest.numInstances();
		this.numClasses = dataTrain.numClasses();
		this.numTerms = dataTrain.numAttributes() - 1;

		BipartiteGraph g = createBipartiteGraph(dataTrain, dataTest);
		
		int op_inicializacao = 0;
		
		int local_max_itr = 100;
		int global_max_itr = 100;
		this.pbg = new PBG(g, alpha, beta, local_max_itr, global_max_itr);
		
		Matrix A = Matrix.createRand(this.numTest + this.numTrain,
				this.numClasses);
		Matrix B = Matrix.createRand(this.numClasses, this.numTerms)
				.transpose();
		
		if (op_inicializacao == 1) {
			initialize1(B, dataTrain);
		}
		
		// atribui classe para a matrix A para exemplos de treino
		HashMap<Integer, Integer> labeled = new HashMap<>();
		for(int inst = this.numTest ; inst < this.numTest + this.numTrain; inst++) {
			int k = (int)dataTrain.instance(this.numTest - inst).classValue();
			Arrays.fill(A.mat[inst], 0);
			A.mat[inst][k] = 1;
			labeled.put(inst, k); 
		}
		
		this.pbg.bgp(A, B, labeled);	// propagacao
		
		super.fUnlabeledDocs = new double[numTest][numClasses];
		
		for(int inst = 0; inst < numTest; inst++)
			fUnlabeledDocs[inst][A.maxIdLine(inst)] = 1;		
	}
	
	public void initialize1(Matrix B, Instances dataTrain) {		
		for (int term = 0; term < numTerms; term++) {
			Arrays.fill(B.mat[term], 0);
			for (int inst=0; inst < numTrain; inst++) {						
				double v = dataTrain.instance(inst).value(term);
				if (v > 0) {					
					int k = (int)dataTrain.instance(inst).classValue();						
					B.mat[term][k] += v;
				}
			}
		}
		B.normalizebycolumn();
	}

	public BipartiteGraph createBipartiteGraph(Instances dataTrain,
			Instances dataTest) {
		BipartiteGraph g = new BipartiteGraph();
		
		for (int inst = 0; inst < numTest; inst++) {
			for (int term = 0; term < numTerms; term++) {
				double v = dataTest.instance(inst).value(term);
				if (v > 0)
					g.add_edge(inst, term, v);
			}
		}

		for (int inst = numTest; inst < numTest + numTrain; inst++) {
			for (int term = 0; term < numTerms; term++) {
				double v = dataTrain.instance(inst).value(term);
				if (v > 0)
					g.add_edge(inst, term, v);
			}
		}
		return g;
	}
}
