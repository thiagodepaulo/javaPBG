package bgp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

public class BipartiteGraph {

	private HashMap<Integer, List<Integer>> typea;
	private HashMap<Integer, List<Integer>> typeb;
	private HashMap<Edge, Double> weight;
	
	public BipartiteGraph() {
		typea = new HashMap<>();
		typeb = new HashMap<>();
		weight = new HashMap<>();
	}
	
	public void add_edge(int a, int b, double w_ab) {
		List<Integer> lista = null;
		List<Integer> listb = null;
		if (typea.containsKey(a)) {
			lista = typea.get(a);
		} else {
			lista = new LinkedList<>();
			typea.put(a, lista);
		}
		if (typeb.containsKey(b)) {
			listb = typeb.get(b);
		} else {
			listb = new LinkedList<>();
			typeb.put(b, listb);
		}		
		lista.add(b);		
		listb.add(a);
		Edge e_ab = new Edge(a, b);
		weight.put(e_ab, w_ab);
	}

	public void add_weight(Edge e, double w) {
		weight.put(e, w);
	}

	public int a_len() {
		return typea.size();
	}

	public int b_len() {
		return typeb.size();
	}

	public int edges_len() {
		return weight.size();
	}

	public List<Integer> a_neig(int a) {
		return typea.get(a);
	}

	public List<Integer> b_neig(int b) {
		return typeb.get(b);
	}

	public double w(int a, int b) {
		return this.weight.get(new Edge(a, b));
	}

	public List<Entry<Integer, Double>> w_a_neig(int a) {
		List<Entry<Integer, Double>> r = new ArrayList<>(typea.size());
		for (int b : a_neig(a)) {
			Edge e = new Edge(a, b);
			PairEntry entry = new PairEntry(b, this.weight.get(e));
			r.add(entry);
		}
		return r;
	}

	public List<Entry<Integer, Double>> w_b_neig(int b) {
		List<Entry<Integer, Double>> r = new ArrayList<>(typeb.size());
		for (int a : b_neig(b)) {
			Edge e = new Edge(a, b);
			PairEntry entry = new PairEntry(a, this.weight.get(e));
			r.add(entry);
		}
		return r;
	}

	public List<Edge> edges() {
		List<Edge> l = new LinkedList<>();
		for (int a : this.a_vertices()) {
			for (Entry<Integer, Double> bw : this.w_a_neig(a)) {
				int b = (int) bw.getKey();
				double w = (double) bw.getValue();
				Edge e = new Edge(a, b, w);
				l.add(e);
			}
		}
		return l;
	}

	public Set<Integer> a_vertices() {
		return typea.keySet();
	}

	public Set<Integer> b_vertices() {
		return typeb.keySet();
	}

	public void load(String arq) throws IOException {
		// Open the file
		FileInputStream fstream = new FileInputStream(arq);
		BufferedReader br = new BufferedReader(new InputStreamReader(fstream));

		String strLine;

		// Read File Line By Line
		while ((strLine = br.readLine()) != null) {
			String[] subs = strLine.split("\\s+");
			int a = Integer.parseInt(subs[0]);
			int b = Integer.parseInt(subs[1]);
			double w = Integer.parseInt(subs[2]);
			this.add_edge(a, b, w);
		}

		// Close the input stream
		br.close();
	}

	public void save(String arq) throws IOException {
		FileOutputStream fstream = new FileOutputStream(arq);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fstream));

		for (Edge e : this.edges()) {
			bw.write(e.a + " " + e.b + " " + e.w + "\n");
		}
		bw.close();
	}

	public String toString() {
		StringBuilder strb = new StringBuilder();
		for(Edge e: this.edges()) {
			strb.append(e);
			strb.append(" ");
		}
		return strb.toString();
	}
	
	public int a_size() {
		return this.a_vertices().size();
	}
	
	public int b_size() {
		return this.b_vertices().size();
	}
}

class Edge {

	public int a;
	public int b;
	public double w;

	public Edge(int a, int b) {
		this(a, b, -1);
	}

	public Edge(int a, int b, double w) {
		this.a = a;
		this.b = b;
		this.w = w;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + this.a;
		result = prime * result + this.b;
		return result;
	}

	public boolean equals(Object obj) {
		if (obj instanceof Edge) {
			Edge edge2 = (Edge) obj;
			return edge2.a == this.a && edge2.b == this.b;
		}
		return false;
	}

	public String toString() {
		return "(" + this.a + ", " + this.b + ", " + this.w + ")";
	}
}

class PairEntry implements Entry<Integer, Double> {

	int key;
	double value;

	public PairEntry(int key, double value) {
		this.key = key;
		this.value = value;
	}

	@Override
	public Integer getKey() {
		return key;
	}

	@Override
	public Double getValue() {
		return value;
	}

	@Override
	public Double setValue(Double value) {
		this.value = value;
		return value;
	}
}