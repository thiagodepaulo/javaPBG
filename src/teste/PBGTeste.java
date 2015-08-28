package teste;

import static java.lang.System.out;

import java.io.IOException;

import bgp.BipartiteGraph;
import bgp.Matrix;
import bgp.PBG;

public class PBGTeste {

	public static void main(String args[]) throws IOException {
		BipartiteGraph g = new BipartiteGraph();

		g.load("/exp/thiagopf/bgp/jBGP/teste.graph");
		
		out.println(g);
		
		Matrix A = Matrix.loadtxt("/exp/thiagopf/bgp/jBGP/A");
		//A.print();		
		Matrix B = Matrix.loadtxt("/exp/thiagopf/bgp/jBGP/B");
		//B.print();
		
		PBG pbg = new PBG(g, 0, 0, 10, 10);
		
		double[] r = pbg.localPropag(0, A.mat[0], B);
		r = pbg.localPropag(0, r, B);
		print(r);
		
		A.print();
		out.println();
		B.print();
		out.println();
		B = pbg.globalPropag(A, B);
		B.print();
		
	}
	
	public static void print(double[] vet) {
		for(double d: vet) {
			out.print(d+" ");
		}
		out.println();
	}

}
