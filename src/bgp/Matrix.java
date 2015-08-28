package bgp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import util.RandomSamplers;

public class Matrix {

	public double[][] mat;

	public Matrix() {

	}

	public Matrix(double[][] mat) {
		this.mat = mat;
	}

	public Matrix transpose() {
		int n = this.mat.length;
		int m = this.mat[0].length;
		double[][] mat_t = new double[m][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				mat_t[j][i] = this.mat[i][j];
			}
		}
		return new Matrix(mat_t);
	}

	public int[] shape() {
		int[] dim = new int[2];
		dim[0] = this.mat.length; // num. rows
		dim[1] = this.mat[0].length; // num columns
		return dim;
	}

	public int maxIdLine(int i) {
		int max_id = 0; 		
		for(int j = 1; j < this.mat[i].length; j++)
			if(this.mat[i][j] > this.mat[i][max_id])
				max_id = j;
		return max_id;
	}

	public int[][] topL(int l) {
		int n = mat.length;
		int tops[][] = new int[mat.length][l];

		for (int i = 0; i < n; i++) {
			int min = 0;
			int last = 0;
			for (int j = 0; j < mat[i].length; j++) {
				double v = mat[i][j];
				if (last < l)
					tops[i][last++] = j;
				else if (v > mat[i][min])
					tops[i][min] = j;
				for (int r = 0; r < l && r < last; r++)
					min = mat[i][tops[i][min]] > mat[i][tops[i][r]] ? r : min;
			}
		}
		return tops;
	}

	public static Matrix createRand(int n, int m) {
		RandomSamplers rand = new RandomSamplers();
		double[] ones = new double[m];
		java.util.Arrays.fill(ones, 1);
		double[][] mat = rand.randDir(ones, n);
		return new Matrix(mat);
	}

	public double[] getCopyLine(int d_j) {
		double[] A_j = mat[d_j];
		double[] copyA_j = new double[A_j.length];
		for (int k = 0; k < A_j.length; k++) {
			copyA_j[k] = A_j[k];
		}
		return copyA_j;
	}

	public static double mean(double[] v1, double[] v2) {
		double l = v1.length;
		double sum = 0;
		for (int i = 0; i < l; i++) {
			sum += Math.abs(v1[i] - v2[i]);
		}
		return sum / l;
	}

	public static double sum(double[] v1) {
		double sum = 0;
		for (double d : v1) {
			sum += d;
		}
		return sum;
	}

	public static double[] hadamard(double[] v1, double[] v2) {
		int size = v1.length;
		double[] r = new double[size];
		for (int i = 0; i < size; i++) {
			r[i] = v1[i] * v2[i];
		}
		return r;
	}

	public static double[] normalize(double[] v) {
		int size = v.length;
		double[] norm = new double[size];
		double sum = 0;
		for (double d : v) {
			sum += d;
		}
		for (int i = 0; i < size; i++) {
			norm[i] = v[i] / sum;
		}
		return norm;
	}

	public static double[] normalize(double[] v, double fac) {
		int size = v.length;
		double[] norm = new double[size];
		double sum = 0;
		for (double d : v) {
			sum += d;
		}
		for (int i = 0; i < size; i++) {
			norm[i] = (fac * v[i]) / sum;
		}
		return norm;
	}

	public void normalizebycolumn() {
		for (int k = 0; k < this.mat[0].length; k++) {
			double sum_k = 0;
			for (int i = 0; i < this.mat.length; i++) {
				sum_k += this.mat[i][k];
			}
			for (int i = 0; i < this.mat.length; i++) {
				this.mat[i][k] /= sum_k;
			}
		}
	}

	public static double prod(double[] v1, double[] v2) {
		double sum = 0;
		for (int i = 0; i < v1.length; i++) {
			sum += v1[i] * v2[i];
		}
		return sum;
	}

	public static double[] prod(double c, double[] v) {
		int size = v.length;
		double[] r = new double[size];
		for (int i = 0; i < size; i++) {
			r[i] = c * v[i];
		}
		return r;
	}

	public static void save(String arq, Matrix A) {
		try {
			FileOutputStream fos = new FileOutputStream(arq);
			GZIPOutputStream gzos = new GZIPOutputStream(fos);
			ObjectOutputStream out = new ObjectOutputStream(gzos);
			out.writeObject(A.mat);
			out.flush();
			out.close();
		} catch (IOException e) {
			System.out.println(e);
		}
	}

	public static Matrix loadtxt(String arq) throws IOException {
		// Open the file
		FileInputStream fstream = new FileInputStream(arq);
		BufferedReader br = new BufferedReader(new InputStreamReader(fstream));

		String strLine;
		HashMap<Integer, String[]> map = new HashMap<>();
		int i = 0;
		// Read File Line By Line
		while ((strLine = br.readLine()) != null) {
			String[] subs = strLine.split("\\s+");
			map.put(i++, subs);
		}
		// Close the input stream
		br.close();

		int n = map.size();
		double[][] mat = new double[n][];
		for (int l = 0; l < n; l++) {
			String[] subs = map.get(l);
			int m = subs.length;
			mat[l] = new double[m];
			for (int j = 0; j < m; j++) {
				mat[l][j] = Double.parseDouble(subs[j]);
			}
		}
		return new Matrix(mat);
	}

	public static void savetxt(String arq, Matrix A)
			throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(new File(arq));
		pw.print(A.toString());
		pw.close();
	}

	public static Matrix load(String arq) {
		try {
			FileInputStream fis = new FileInputStream(arq);
			GZIPInputStream gzis = new GZIPInputStream(fis);
			ObjectInputStream in = new ObjectInputStream(gzis);
			double[][] A = (double[][]) in.readObject();
			in.close();
			return new Matrix(A);
		} catch (Exception e) {
			System.out.println(e);
		}
		return null;
	}

	public String toString() {
		StringBuilder str = new StringBuilder();
		String nl = "\n";
		for (int i = 0; i < this.mat.length; i++) {
			for (int j = 0; j < this.mat[i].length; j++) {
				str.append(this.mat[i][j]);
				str.append(" ");
			}
			str.append(nl);
		}
		return str.toString();
	}

	public void print() {
		System.out.println(toString());
	}

	public void printLine(int i) {
		for (int j = 0; j < mat[i].length; j++)
			System.out.print(mat[i][j] + " ");
		System.out.println();
	}

}
