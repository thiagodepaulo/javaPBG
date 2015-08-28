package bgp;

import static bgp.Matrix.hadamard;
import static bgp.Matrix.mean;
import static bgp.Matrix.normalize;
import static bgp.Matrix.prod;
import static bgp.Matrix.save;
import static java.lang.Math.log;
import static java.lang.System.out;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

public class PBG {

	protected final BipartiteGraph G;
	protected final int LOCAL_MAX_ITR;
	protected final int GLOBAL_MAX_ITR;
	protected final double LOCAL_CONV_THRESHOLD;
	protected final double GLOBAL_CONV_THRESHOLD;
	protected final double alpha;
	protected final double beta;

	// private final double AUX = 1000;

	public boolean save = false;
	public String out_dir;

	public PBG(BipartiteGraph G, double alpha, double beta, int local_max_itr,
			int global_max_itr) {
		this(G, alpha, beta, local_max_itr, global_max_itr, 1e-6, 1e-6);
	}

	public PBG(BipartiteGraph G, double alpha, double beta, int local_max_itr,
			int global_max_itr, double local_threshold, double global_thresold) {
		this.G = G;
		this.alpha = alpha;
		this.beta = beta;
		this.LOCAL_MAX_ITR = local_max_itr;
		this.GLOBAL_MAX_ITR = global_max_itr;
		this.LOCAL_CONV_THRESHOLD = local_threshold;
		this.GLOBAL_CONV_THRESHOLD = global_thresold;
	}

	public void bgp(Matrix A, Matrix B) {
		double oldq = Double.POSITIVE_INFINITY;
		int global_niter = 0;
		int local_niter;
		double[] oldA_j = null;
		double mean_change;
		double q;

		while (global_niter < this.GLOBAL_MAX_ITR) {
			global_niter++;
			for (int d_j : G.a_vertices()) {
				if (d_j % 1000 == 0)
					System.out.println(d_j);
				local_niter = 0;
				while (local_niter < this.LOCAL_MAX_ITR) {
					local_niter++;
					oldA_j = A.mat[d_j]; // A.getCopyLine(d_j);
					A.mat[d_j] = this.localPropag(d_j, A.mat[d_j], B);
					mean_change = mean(A.mat[d_j], oldA_j);
					if (mean_change <= this.LOCAL_CONV_THRESHOLD) {
						break;
					}
				}
			}

			B = this.globalPropag(A, B);

			q = this.Q(A, B);
			out.println("itr=" + global_niter + " Q=" + q);
			if (Math.abs(q - oldq) <= this.GLOBAL_CONV_THRESHOLD) {
				out.println("Global convergence in " + global_niter
						+ " iterations");
				break;
			}
			oldq = q;
			if (this.save) {
				save(this.out_dir + "/A_" + global_niter + ".npj", A);
				save(this.out_dir + "/B_" + global_niter + ".npj", B);
			}
		}
	} 
	
	public void suppress(double[] vet, int k) {
		for (int i = 0; i < vet.length; i++)
			vet[i] = (i != k) ? 0 : vet[i];
	}

	public void bgp(Matrix A, Matrix B, HashMap<Integer, Integer> labeled) {
		double oldq = Double.POSITIVE_INFINITY;
		int global_niter = 0;
		int local_niter;
		double[] oldA_j = null;
		double mean_change;
		double q;

		while (global_niter < this.GLOBAL_MAX_ITR) {
			global_niter++;
			for (int d_j : G.a_vertices()) {
				if (d_j % 1000 == 0)
					System.out.println(d_j);
				local_niter = 0;
				while (local_niter < this.LOCAL_MAX_ITR) {
					local_niter++;
					oldA_j = A.mat[d_j]; // A.getCopyLine(d_j);
					A.mat[d_j] = this.localPropag(d_j, A.mat[d_j], B);
					mean_change = mean(A.mat[d_j], oldA_j);
					if (mean_change <= this.LOCAL_CONV_THRESHOLD) {
						break;
					}
				}
				if (labeled.containsKey(d_j))
					suppress(A.mat[d_j], labeled.get(d_j));
			}

			B = this.globalPropag(A, B);

			q = this.Q(A, B);
			out.println("itr=" + global_niter + " Q=" + q);
			if (Math.abs(q - oldq) <= this.GLOBAL_CONV_THRESHOLD) {
				out.println("Global convergence in " + global_niter
						+ " iterations");
				break;
			}
			oldq = q;
			if (this.save) {
				save(this.out_dir + "/A_" + global_niter + ".npj", A);
				save(this.out_dir + "/B_" + global_niter + ".npj", B);
			}
		}
	}	

	public double Q(Matrix A, Matrix B) {
		int w_i;
		double f_ij;
		double sumAjBi;
		double _sum = 0;
		for (int d_j : G.a_vertices()) {
			for (Entry<Integer, Double> e : G.w_a_neig(d_j)) {
				w_i = e.getKey();
				f_ij = e.getValue();
				sumAjBi = prod(A.mat[d_j], B.mat[w_i]);
				_sum += f_ij * log(f_ij / sumAjBi) - f_ij + sumAjBi;
			}
		}
		return _sum;
	}

	public Matrix globalPropag(Matrix A, Matrix B) {
		List<Entry<Integer, Double>> neig;
		int d_j;
		double f_ji;
		double[] C;
		Matrix rB = new Matrix(new double[B.shape()[0]][B.shape()[1]]);
		for (int w_i : G.b_vertices()) {
			neig = G.w_b_neig(w_i);
			for (Entry<Integer, Double> e : neig) {
				d_j = e.getKey();
				f_ji = e.getValue();

				C = normalize(hadamard(A.mat[d_j], B.mat[w_i]));
				for (int k = 0; k < B.mat[w_i].length; k++) {
					rB.mat[w_i][k] += this.beta + f_ji * C[k];
				}
			}
		}
		rB.normalizebycolumn();
		return rB;
	}

	public double[] localPropag(int d_j, double[] A_j, Matrix B) {
		List<Entry<Integer, Double>> neig = G.w_a_neig(d_j);
		double[] C;
		int w_i;
		double f_ji;
		double[] rA_j = new double[A_j.length];
		for (Entry<Integer, Double> e : neig) {
			w_i = e.getKey();
			f_ji = e.getValue();

			C = normalize(hadamard(A_j, B.mat[w_i]));
			for (int k = 0; k < A_j.length; k++) {
				rA_j[k] += this.alpha + f_ji * C[k];
			}
		}
		return rA_j;
	}

	public static void main(String args[]) throws IOException {

		if (args.length < 7)
			System.out
					.println("usage: <arq_graph> <k> <alpha> <beta> <local_niter> <global_niter> <out>");
		String arq_graph = args[0];
		int k = Integer.parseInt(args[1]);
		double alpha = Double.parseDouble(args[2]);
		double beta = Double.parseDouble(args[3]);
		int local_niter = Integer.parseInt(args[4]);
		int global_niter = Integer.parseInt(args[5]);
		String out = args[6];

		BipartiteGraph g = new BipartiteGraph();
		System.out.println("Carregando Grafo");
		g.load(arq_graph);
		System.out.println("pronto");

		PBG pbg = new PBG(g, alpha, beta, local_niter, global_niter);

		Matrix A = Matrix.createRand(g.a_size(), k);
		Matrix B = Matrix.createRand(k, g.b_size()).transpose();

		pbg.bgp(A, B);

		Matrix.savetxt(out + "/A.jnp", A);
		Matrix.savetxt(out + "/B.jnp", B);
	}

}
