package teste;

import java.io.FileNotFoundException;
import java.io.IOException;

import bgp.Matrix;

public class MatrixTeste {
	
	public static void main(String args[]) throws IOException {		
		
		teste4();
		
	}
	
	public static void teste1() {
		Matrix m = new Matrix();
		m.mat = new double[][]{{1,2},{3,1}};
		double r = Matrix.mean(m.mat[0], m.mat[1]);
		
		System.out.println(r);
		System.out.println(m);
		
		Matrix.save("a", m);
		
		Matrix m2 = Matrix.load("a");
		System.out.println(m2);
		
		double[] a = new double[3];
		
		for(double d: a)
			System.out.println(d);
	}
	
	public static void teste2() {
		Matrix m = Matrix.createRand(2, 3);
		System.out.println(m);						
	}
	
	public static void teste3() {
		Matrix m = new Matrix();
		m.mat = new double[][]{{1,2,8,6,4,7,-1},{9,7,3,5,8,2,100}};
		
		int[][] r = m.topL(2);
		for(int i=0; i<r.length; i++) {
			for (int j=0; j<r[i].length; j++) {
				System.out.print(r[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	public static void teste4() throws IOException {
		Matrix m = new Matrix();
		m.mat = new double[][]{{1,2},{3,1}};
		double r = Matrix.mean(m.mat[0], m.mat[1]);
		
		System.out.println(r);
		System.out.println(m);
		
		Matrix.savetxt("a", m);
		
		Matrix m2 = Matrix.loadtxt("a");
		System.out.println(m2);
				
	}

}
