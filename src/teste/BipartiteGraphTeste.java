package teste;

import java.util.Map.Entry;

import bgp.BipartiteGraph;

public class BipartiteGraphTeste {
	
	public static void main(String[] args) {
		BipartiteGraph g = new BipartiteGraph();
		
		g.add_edge(1, 2, 4.5);
		g.add_edge(3,4,6.90);
		
		System.out.println("oi");
		System.out.println(g);
		
		for (int b: g.a_neig(1)) {
			System.out.println(b);
		}
		
		for (Entry<Integer, Double> e: g.w_a_neig(1)) {
			System.out.println(e.getKey());
			System.out.println(e.getValue());
		}
	}

}
