package phylonet.coalescent;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.lang.reflect.Array;
import java.util.AbstractMap;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.LinkedList;
import java.util.Queue;
import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;
import phylonet.coalescent.MGDInference_DP.TaxonNameMap;
import phylonet.coalescent.TreeAllCluster.FrequencyData;

public class DuplicationWeightCounter {

	HashMap<STBipartition, Integer> weights;
	String[] gtTaxa;
	String[] stTaxa;


	// private List<Set<STBipartition>> X;

	private Map<STBipartition, Integer> geneTreeSTBCount;
	private Map<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer> geneTreeInvalidSTBCont;

	private boolean rooted;

	private TaxonNameMap taxonNameMap;

	private ClusterCollection clusters;

	// Used for dup loss calculations
	List<STITreeCluster> treeAlls = new ArrayList<STITreeCluster>();

	HashMap<STBipartition, Set<STBipartition>> alreadyWeigthProcessed = new HashMap<STBipartition, Set<STBipartition>>();

	public DuplicationWeightCounter(String[] gtTaxa, String[] stTaxa,
			boolean rooted, TaxonNameMap taxonNameMap,
			ClusterCollection clusters2) {
		this.gtTaxa = gtTaxa;
		this.stTaxa = stTaxa;
		this.rooted = rooted;
		this.taxonNameMap = taxonNameMap;
		this.clusters = clusters2;
	}

	private String getSpeciesName(String geneName) {
		String stName = geneName;
		if (taxonNameMap != null) {
			stName = taxonNameMap.getTaxonName(geneName);
		}
		return stName;
	}

	private boolean addToClusters(STITreeCluster c, int size,
			boolean geneTreeCluster) {
		Vertex nv = c.new Vertex(); 
		//System.out.println("nv = " + nv + "size : "+  size + " " + c);
		return clusters.addCluster(nv, size);
	}


	private class NodeData{
		Double mainF , alt1F , alt2F ; 
		Integer effn ; 
	}

	
	private NodeData getNodeData(Double m, Double a1, Double a2, Integer en) {
		NodeData nd;
		nd = new NodeData();
		nd.mainF = m;
		nd.alt1F=a1;
		nd.alt2F=a2;
		nd.effn = en;
		return nd;
	}


	Integer calculateHomomorphicCost(List<Integer> El, STITreeCluster cluster,
			Vertex smallV, Vertex bigv, List<Tree> trees) {
		Integer e = 0;
		for (int k = 0; k < trees.size(); k++) {
			Tree tr = trees.get(k);
			STITreeCluster treeAll = treeAlls.get(k);
			if (smallV.getCluster().isDisjoint(treeAll)
					|| bigv.getCluster().isDisjoint(treeAll)) {
				continue;
			}
			if (El.get(k) == null) {
				if (taxonNameMap == null) {
					El.set(k, DeepCoalescencesCounter.getClusterCoalNum_rooted(
							tr, cluster));
				} else {
					El.set(k, DeepCoalescencesCounter.getClusterCoalNum_rooted(
							tr, cluster, taxonNameMap));
				}
			}
			e += El.get(k);
		}
		return e;
	}

	/*
	 * Calculates the cost of a cluster based on the standard definition
	 * */
	int calculateDLstdClusterCost(STITreeCluster cluster, List<Tree> trees) {
		/*
		 * if (XLweights.containsKey(cluster)) { return XLweights.get(cluster);
		 * }
		 */
		int weight = 0;
		for (Entry<STBipartition, Integer> entry : geneTreeSTBCount.entrySet()) {
			STBipartition otherSTB = entry.getKey();
			/*if (cluster.containsCluster(otherSTB.c))
				continue;*/
			boolean c1 = cluster.containsCluster(otherSTB.cluster1);
			boolean c2 = cluster.containsCluster(otherSTB.cluster2);
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (Entry<SimpleEntry<STITreeCluster, STITreeCluster>, Integer> entry : geneTreeInvalidSTBCont.entrySet()) {
			SimpleEntry<STITreeCluster, STITreeCluster> otherSTB = entry.getKey();
			boolean c1 = cluster.containsCluster(otherSTB.getKey());
			boolean c2 = cluster.containsCluster(otherSTB.getValue());
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (STITreeCluster treeAll : treeAlls) {
			if (cluster.containsCluster(treeAll)) {
				weight++;
			}
		}
		int ret = weight;
		// XLweights.put(cluster, ret);
		return ret;
	}





	double computeTreeSTBipartitions(MGDInference_DP inference) {
		System.out.println("Inside computeTreeSTBipartitions in DuplicationWeightCounter");


		double unweigthedConstant = 0;
		double weightedConstant = 0;
		int k = inference.trees.size(); 
		System.out.println("k / inference trees size = " + k); // number of gene trees
		String[] leaves = stTaxa;
		int n = leaves.length;
		System.out.println("n / leaves length / StTaxa size = " + n);
		boolean duploss = (inference.optimizeDuploss == 3);		

		geneTreeSTBCount = new HashMap<STBipartition, Integer>(k * n); // rooted tree with n leaves has n bipartitions for each node in each gene tree 

		geneTreeInvalidSTBCont = new HashMap<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer>();
		// geneTreeRootSTBs = new HashMap<STBipartition, Integer>(k*n);
		// needed for fast version
		// clusterToSTBs = new HashMap<STITreeCluster, Set<STBipartition>>(k*n);

		STITreeCluster all = new STITreeCluster(stTaxa);
		String as[];
		int j = (as = stTaxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i]; 
			System.out.println(t);
			all.addLeaf(t);
		}
		addToClusters(all, leaves.length, false);


		if(inference.STBforST == false)
		{

	
			for (int t = 0; t < inference.trees.size(); t++) {
				Tree tr = inference.trees.get(t);


				STITreeCluster allInducedByGT = new STITreeCluster(stTaxa);

				String[] gtLeaves = tr.getLeaves();
				for (int i = 0; i < gtLeaves.length; i++) {
					String l = gtLeaves[i];
					allInducedByGT.addLeaf(getSpeciesName(l));
				}
				for (String l : allInducedByGT.getTaxa()){
					System.out.println("allInducedByGT taxa: "+l);
				}

				treeAlls.add(allInducedByGT);
				int allInducedByGTSize = allInducedByGT.getClusterSize();
				//System.out.println("allInducedByGTSize = " + allInducedByGTSize);

				weightedConstant += duploss ? 
						2 * (allInducedByGTSize - 1) : 0;
						
				unweigthedConstant += (tr.getLeafCount() - 1);

				Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
				// Map<TNode,STITreeCluster> nodeToGTCluster = new HashMap<TNode,
				// STITreeCluster>(n);

				for (TNode node : tr.postTraverse()) {				
					//System.err.println("Node is:" + node);
					if (node.isLeaf()) {
						String nodeName = node.getName();

						// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
						// gtCluster.addLeaf(nodeName);

						nodeName = getSpeciesName(nodeName);
						//System.out.println("nodename" + nodeName);
						STITreeCluster cluster = new STITreeCluster(leaves);
						cluster.addLeaf(nodeName);

						addToClusters(cluster, 1, true);

						nodeToSTCluster.put(node, cluster);

						// nodeToGTCluster.put(node, gtCluster);

						if (!rooted) {
							throw new RuntimeException("Unrooted not implemented.");
							// tryAddingSTB(cluster, treeComplementary(null
							// /*gtCluster*/,leaves), null, node, true);
						}
					} 
					else {
						//System.out.println("This is not a leaf node");
						int childCount = node.getChildCount();
						//System.out.println(childCount);
						STITreeCluster childbslist[] = new STITreeCluster[childCount];
						BitSet bs = new BitSet(leaves.length); // This BitSet will be updated by performing OR operations with the BitSet of each child.

						// BitSet gbs = new BitSet(leaves.length);
						int index = 0;
						for (TNode child: node.getChildren()) {
							childbslist[index++] = nodeToSTCluster.get(child);
							bs.or(nodeToSTCluster.get(child).getBitSet());
							// gbs.or(nodeToGTCluster.get(child).getBitSet());
						}

						// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
						// gtCluster.setCluster(gbs);

						STITreeCluster cluster = new STITreeCluster(leaves);
						cluster.setCluster((BitSet) bs.clone());

						int size = cluster.getClusterSize();

						addToClusters(cluster, size, true);
						nodeToSTCluster.put(node, cluster);
						// nodeToGTCluster.put(node, gtCluster);

						if (rooted) {

							if (index > 2) {
								throw new RuntimeException(
										"None bifurcating tree: " + tr + "\n"
												+ node);
							}

							STITreeCluster l_cluster = childbslist[0];
							STITreeCluster r_cluster = childbslist[1];
							// System.out.println("l_cluster: " + l_cluster);
							// System.out.println("r_cluster: " + r_cluster);
							// System.out.println("cluster: " + cluster);
							// System.out.println("node: " + node);

							tryAddingSTB(l_cluster, r_cluster, cluster, node, true);

						} else {
							throw new RuntimeException("Unrooted not implemented.");
							/*
							* if (childCount == 2) { STITreeCluster l_cluster =
							* childbslist[0];
							* 
							* STITreeCluster r_cluster = childbslist[1];
							* 
							* STITreeCluster allMinuslAndr_cluster =
							* treeComplementary(null this should be
							* gtCluster?,leaves);
							* 
							* STITreeCluster lAndr_cluster = cluster;
							* 
							* if (allMinuslAndr_cluster.getClusterSize() != 0) { //
							* add Vertex STBs tryAddingSTB(l_cluster, r_cluster,
							* cluster, node, true); tryAddingSTB( r_cluster,
							* allMinuslAndr_cluster, null, node, true);
							* tryAddingSTB(l_cluster, allMinuslAndr_cluster, null,
							* node, true);
							* 
							* // Add the Edge STB tryAddingSTB(lAndr_cluster,
							* allMinuslAndr_cluster, null, node, true); }
							* 
							* } else if (childCount == 3 && node.isRoot()) {
							* STITreeCluster l_cluster = childbslist[0];
							* 
							* STITreeCluster m_cluster = childbslist[1];
							* 
							* STITreeCluster r_cluster = childbslist[2];
							* 
							* tryAddingSTB(l_cluster, r_cluster, null, node, true);
							* tryAddingSTB(r_cluster, m_cluster, null, node, true);
							* tryAddingSTB(l_cluster, m_cluster, null, node, true);
							* } else { throw new
							* RuntimeException("None bifurcating tree: "+ tr+ "\n"
							* + node); }
							*/}
					}
				}

			}

			int s = 0;
			for (Integer c : geneTreeSTBCount.values()) {
				s += c;
			}
			// System.err.println("STBs in gene trees (count): "
			// 		+ geneTreeSTBCount.size());
			// System.err.println("STBs in gene trees (sum): " + s);

			s = clusters.getClusterCount();

			//System.err.println("Number of Clusters: " + s);

			weights = new HashMap<STBipartition, Integer>(
					geneTreeSTBCount.size() * 2);
			// System.err.println("sigma n is "+sigmaN);

			if (inference.DLbdWeigth == -1) {			
				inference.DLbdWeigth = (weightedConstant + 2*k + 0.0D) / 2*(k*n);
				System.out.println("Estimated bd weight = " + inference.DLbdWeigth);
			}
				
			return (unweigthedConstant + (1 - inference.DLbdWeigth) * weightedConstant);
	}
	else {


		System.out.println("************For the Species Tree************");
		Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);

		for (TNode node : inference.SpeciesTree.postTraverse()) {				
			//System.err.println("Node is:" + node);
			if (node.isLeaf()) {
				String nodeName = node.getName();

				// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
				// gtCluster.addLeaf(nodeName);

				nodeName = getSpeciesName(nodeName);
				//System.out.println("nodename" + nodeName);
				STITreeCluster cluster = new STITreeCluster(leaves);
				cluster.addLeaf(nodeName);

				addToClusters(cluster, 1, true);

				nodeToSTCluster.put(node, cluster);

				// nodeToGTCluster.put(node, gtCluster);

				if (!rooted) {
					throw new RuntimeException("Unrooted not implemented.");
					// tryAddingSTB(cluster, treeComplementary(null
					// /*gtCluster*/,leaves), null, node, true);
				}
			} 
			else {
				//System.out.println("This is not a leaf node");
				int childCount = node.getChildCount();
				//System.out.println(childCount);
				STITreeCluster childbslist[] = new STITreeCluster[childCount];
				BitSet bs = new BitSet(leaves.length); // This BitSet will be updated by performing OR operations with the BitSet of each child.

				// BitSet gbs = new BitSet(leaves.length);
				int index = 0;
				for (TNode child: node.getChildren()) {
					childbslist[index++] = nodeToSTCluster.get(child);
					bs.or(nodeToSTCluster.get(child).getBitSet());
					// gbs.or(nodeToGTCluster.get(child).getBitSet());
				}

				// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
				// gtCluster.setCluster(gbs);

				STITreeCluster cluster = new STITreeCluster(leaves);
				cluster.setCluster((BitSet) bs.clone());

				int size = cluster.getClusterSize();

				addToClusters(cluster, size, true);
				nodeToSTCluster.put(node, cluster);
				// nodeToGTCluster.put(node, gtCluster);

				if (rooted) {

					if (index > 2) {
						throw new RuntimeException(
								"None bifurcating tree: " + "\n"
										+ node);
					}

					STITreeCluster l_cluster = childbslist[0];
					STITreeCluster r_cluster = childbslist[1];
					// System.out.println("l_cluster: " + l_cluster);
					// System.out.println("r_cluster: " + r_cluster);
					// System.out.println("cluster: " + cluster);
					// System.out.println("node: " + node);

					tryAddingSTB(l_cluster, r_cluster, cluster, node, true);

				}
			}
		}

		System.out.println("Finished adding STBs for the species tree");

	}
		return 0.0;
	}


	double computeTripartitionsforST(MGDInference_DP inference)
	{
		System.out.println("Inside computeTripartitionsforST in DuplicationWeightCounter");
		DecimalFormat df ;
		df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		DecimalFormatSymbols dfs = DecimalFormatSymbols.getInstance();
		dfs.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(dfs);

		Queue<NodeData> NodeDataList = new LinkedList<NodeData>(); 
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		for (TNode node : inference.SpeciesTree.postTraverse()) 
		{
			if(node.isLeaf())
			{
				// System.out.println("Leaf Node: " + node.getName());
				STITreeCluster cluster = new STITreeCluster(stTaxa);
				cluster.addLeaf(getSpeciesName(node.getName()));
				stack.add(cluster);
			}
			else
			{
				NodeData nd = null;


				if(node.isRoot()){
					NodeDataList.add(null) ; 
					
				}
				// System.out.println("Internal Node: ");
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(stTaxa.length);
				for (TNode child: node.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				STITreeCluster cluster = new STITreeCluster(stTaxa);
				cluster.setCluster((BitSet) bs.clone());
				stack.add(cluster);
				STITreeCluster remaining = cluster.complementaryCluster();
				if(remaining.getClusterSize() != 0)
				{
					childbslist.add(remaining);
				}
				// System.out.println("Childbslist size: " + childbslist.size());
				for (int i = 0; i < childbslist.size(); i++) {
					for (int j = i+1; j < childbslist.size(); j++) {
						for (int k = j+1; k < childbslist.size(); k++) {
							
							Tripartition[] threeTripartitions = new Tripartition[] {
								new Tripartition(childbslist.get(i), childbslist.get(j), childbslist.get(k)) ,
								new Tripartition(childbslist.get(j), childbslist.get(k), childbslist.get(i)) , 
								new Tripartition(childbslist.get(k), childbslist.get(i), childbslist.get(j))
							};
							// Tripartition tri = new Tripartition(childbslist.get(i), childbslist.get(j), childbslist.get(k));
							// System.out.println(inference.BLUE+ "Tripartition: " + tri + inference.RESET);

							TreeAllCluster tas = inference.tas ; 
							FrequencyData f = tas.WeightCalculator(threeTripartitions) ; 
							System.out.println(inference.YELLOW+ "Tripartition:  " + f.freq_s[0] + " " + f.freq_s[1] + " " + f.freq_s[2] + " " + f.effn + inference.RESET);
							nd = getNodeData(f.freq_s[0], f.freq_s[1], f.freq_s[2], f.effn) ; 
							NodeDataList.add(nd);





							
						}
					}					       
				}
			}
		}

		NodeData nd = null ; 
		for (TNode n : inference.SpeciesTree.postTraverse())
		{
			STINode node = (STINode) n ;
			if(node.isLeaf()){
				node.setData(null);
				continue;
			}
			nd = NodeDataList.poll();
			if(nd==null){
				node.setData(null);
				continue;
			}
			Double f1 = nd.mainF ;
			Double f2 = nd.alt1F ;
			Double f3 = nd.alt2F ;
			Double effni = nd.effn + 0.0 ;
			//System.out.println("Before going to final post " + f1 + " " + f2 + " " + f3) ;
			Posterior post = new Posterior(f1 ,f2,f3,effni, 0.5) ; 
			double postT1 = post.getPost() ; 
			post = new Posterior(f2 ,f1,f3,effni, 0.5) ;
			double postT2 = post.getPost() ;
			post = new Posterior(f3 ,f1,f2,effni, 0.5) ;
			double postT3 = post.getPost() ;
			node.setData(":'[pp1="+df.format(postT1)+";pp2="+df.format(postT2)+";pp3="+df.format(postT3)+"]'");
		}


		for(TNode node: inference.SpeciesTree.postTraverse()) {
			STINode stnode = (STINode) node;
			if(stnode.getData() == null){
				continue;
			}
			System.out.println(inference.GREEN+stnode.getData()+inference.RESET);
		}
		return 2.0 ; 
	}

	void addAllPossibleSubClusters(STITreeCluster cluster, int size) {
		BitSet bs = (BitSet) cluster.getBitSet().clone();
		bs.clear(0, size);
		while (true) {
			int tsb = bs.nextClearBit(0);
			if (tsb >= size) {
				break;
			}
			bs.set(tsb);
			bs.clear(0, tsb);
			STITreeCluster c = new STITreeCluster(cluster.getTaxa());
			c.setCluster((BitSet) bs.clone());
			addToClusters(c, c.getClusterSize(), false);
		}
		System.err
				.println("Number of Clusters After Adding All possible clusters: "
						+ clusters.getClusterCount());
	}

	void addExtraBipartitionsByInput(ClusterCollection extraClusters,
			List<Tree> trees, boolean extraTreeRooted) {

		String[] leaves = stTaxa;
		int n = leaves.length;

		// STITreeCluster all = extraClusters.getTopVertex().getCluster();

		for (Tree tr : trees) {

			/*
			 * String[] treeLeaves = tr.getLeaves(); STITreeCluster treeAll =
			 * new STITreeCluster(treeLeaves); for (int i = 0; i <
			 * treeLeaves.length; i++) { String l = treeLeaves[i];
			 * treeAll.addLeaf(l); }
			 */
			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(
					n);

			for (Iterator<TNode> nodeIt = tr.postTraverse().iterator(); nodeIt
					.hasNext();) {
				TNode node = nodeIt.next();
				if (node.isLeaf()) {
					String treeName = node.getName();
					String nodeName = getSpeciesName(treeName);

					STITreeCluster tb = new STITreeCluster(leaves);
					tb.addLeaf(nodeName);

					nodeToSTCluster.put(node, tb);

					if (!extraTreeRooted) {
						// TODO: fix the following (first null)
						throw new RuntimeException("Unrooted not implemented.");
						// tryAddingSTB(tb, tb.complementaryCluster(),
						// all,node,false);
					}
				} else {
					int childCount = node.getChildCount();
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(leaves.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster(leaves);
					cluster.setCluster((BitSet) bs.clone());

					addToClusters(cluster, cluster.getClusterSize(), false);
					nodeToSTCluster.put(node, cluster);

					if (extraTreeRooted) {

						if (index > 2) {
							throw new RuntimeException(
									"None bifurcating tree: " + tr + "\n"
											+ node);
						}

						STITreeCluster l_cluster = childbslist[0];

						STITreeCluster r_cluster = childbslist[1];

						tryAddingSTB(l_cluster, r_cluster, cluster, node, false);
					} else {
						throw new RuntimeException("Unrooted not implemented.");
						/*
						 * if (childCount == 2) { STITreeCluster l_cluster =
						 * childbslist[0];
						 * 
						 * STITreeCluster r_cluster = childbslist[1]; // Fix the
						 * following (first null) STITreeCluster
						 * allMinuslAndr_cluster = treeComplementary(null,
						 * leaves);
						 * 
						 * STITreeCluster lAndr_cluster = cluster;
						 * 
						 * // add Vertex STBs tryAddingSTB( l_cluster,
						 * r_cluster, cluster,node,false); if
						 * (allMinuslAndr_cluster.getClusterSize() != 0) {
						 * tryAddingSTB( r_cluster, allMinuslAndr_cluster,
						 * null,node,false); tryAddingSTB( l_cluster,
						 * allMinuslAndr_cluster, null,node,false);
						 * tryAddingSTB( lAndr_cluster, allMinuslAndr_cluster,
						 * all,node,false); }
						 * 
						 * } else if (childCount == 3 && node.isRoot()) {
						 * STITreeCluster l_cluster = childbslist[0];
						 * 
						 * STITreeCluster m_cluster = childbslist[1];
						 * 
						 * STITreeCluster r_cluster = childbslist[2];
						 * 
						 * tryAddingSTB( l_cluster, r_cluster, null,node,false);
						 * tryAddingSTB( r_cluster, m_cluster, null,node,false);
						 * tryAddingSTB( l_cluster, m_cluster, null,node,false);
						 * } else { throw new
						 * RuntimeException("None bifurcating tree: "+ tr+ "\n"
						 * + node); }
						 */
					}
				}
			}

		}
		int s = extraClusters.getClusterCount();
		/*
		 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
		 */
		System.err
				.println("Number of Clusters After additions from extra Trees: "
						+ s);
	}

	private void tryAddingSTB(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster cluster, TNode node,
			boolean fromGeneTrees) {
		// System.err.println("before adding: " + STBCountInGeneTrees);
		// System.err.println("Trying: " + l_cluster + "|" + r_cluster);
		int size = cluster.getClusterSize();
		if (l_cluster.isDisjoint(r_cluster)) {

			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);
			((STINode) node).setData(stb);
			if (fromGeneTrees) {
				clusters.addGeneTreeSTB(stb, size);
				// gtNodeToSTBs.put(node,stb);
				// addSTBToX(stb,size);
				// System.out.println(stb + " hashes to " + stb.hashCode());
				// if (! hash.containsKey(stb.hashCode()))
				// hash.put(stb.hashCode(), new HashSet<STBipartition>());
				// hash.get(stb.hashCode()).add(stb);
				geneTreeSTBCount.put(
						stb,
						geneTreeSTBCount.containsKey(stb) ? geneTreeSTBCount
								.get(stb) + 1 : 1);
				//*It also updates the geneTreeSTBCount, which keeps track of how many times a particular bipartition has been observed.
			}

			/*
			 * if (size == allInducedByGTSize){ if (!
			 * geneTreeRootSTBs.containsKey(stb)) { geneTreeRootSTBs.put(stb,
			 * 1); } else { geneTreeRootSTBs.put(stb,
			 * geneTreeRootSTBs.get(stb)+1); } }
			 */
		} else {
			System.out.println("Not disjoint");
			AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster> stb = 
					l_cluster.getBitSet().cardinality() > r_cluster.getBitSet().cardinality() ? 
					new AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>(l_cluster, r_cluster) : 
					new AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>(r_cluster, l_cluster);
			geneTreeInvalidSTBCont.put(stb,	geneTreeInvalidSTBCont.containsKey(stb) ? 
					geneTreeInvalidSTBCont.get(stb) + 1 : 
					1);
			// System.err.println("Adding only to extra");
			// This case could happen for multiple-copy
			BitSet and = (BitSet) l_cluster.getBitSet().clone();
			and.and(r_cluster.getBitSet());

			BitSet l_Minus_r = (BitSet) and.clone();
			l_Minus_r.xor(l_cluster.getBitSet());
			STITreeCluster lmr = new STITreeCluster(stTaxa);
			lmr.setCluster(l_Minus_r);

			BitSet r_Minus_l = (BitSet) and.clone();
			r_Minus_l.xor(r_cluster.getBitSet());
			STITreeCluster rml = new STITreeCluster(stTaxa);
			rml.setCluster(r_Minus_l);

			if (!rml.getBitSet().isEmpty()) {
				addToClusters(rml, rml.getClusterSize(), false);
				// addSTBToX( new STBipartition(l_cluster, rml, cluster),size);
			}
			if (!lmr.getBitSet().isEmpty()) {
				addToClusters(lmr, lmr.getClusterSize(), false);
				// addSTBToX(new STBipartition(lmr, r_cluster, cluster), size);
			}
		}
	}

	public Integer getCalculatedBiPartitionDPWeight(STBipartition bi) {
		if (!weights.containsKey(bi)) {
			// weights.put(bi,calculateMissingWeight(bi));
			return null;
		}
		return weights.get(bi);
	}

	// static public int cnt = 0;

	void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {

		if (rooted && taxonNameMap == null && stTaxa.length > trees.size()) {
			//calculateWeightsByLCA(trees, trees);
			if (extraTrees != null) {
				//calculateWeightsByLCA(extraTrees, trees);
			}
		}

	}

	void calculateWeightsByLCA(List<Tree> stTrees, List<Tree> gtTrees) {

		for (Tree stTree : stTrees) {
			SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(stTree);
			for (Tree gtTree : gtTrees) {
				Stack<TNode> stack = new Stack<TNode>();
				for (TNode gtNode : gtTree.postTraverse()) {
					if (gtNode.isLeaf()) {
						stack.push(stTree.getNode(gtNode.getName()));
					} else {
						TNode rightLCA = stack.pop();
						TNode leftLCA = stack.pop();
						// If gene trees are incomplete, we can have this case
						if (rightLCA == null || leftLCA == null) {
							stack.push(null);
							continue;
						}
						TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
						stack.push(lca);
						if (lca != leftLCA && lca != rightLCA) {
							// LCA in stTree dominates gtNode in gene tree
							// gtTree
							STBipartition stSTB = (STBipartition) ((STINode) lca)
									.getData();
							STBipartition gtSTB = (STBipartition) ((STINode) gtNode)
									.getData();
							Set<STBipartition> alreadyProcessedSTBs = alreadyWeigthProcessed
									.get(gtSTB);

							if (alreadyProcessedSTBs == null) {
								alreadyProcessedSTBs = new HashSet<STBipartition>(
										gtTrees.size() / 4);
								alreadyWeigthProcessed.put(gtSTB,
										alreadyProcessedSTBs);
							}

							if (alreadyProcessedSTBs.contains(stSTB)) {
								continue;
							}

							weights.put(
									stSTB,
									(weights.containsKey(stSTB) ? weights
											.get(stSTB) : 0)
											+ geneTreeSTBCount.get(gtSTB));
							alreadyProcessedSTBs.add(stSTB);
						}
					}
				}
			}
		}
	}

	class CalculateWeightTask {

		/**
		 * 
		 */
		private static final long serialVersionUID = -2614161117603289345L;
		private STBipartition stb;
		private ClusterCollection containedClusterCollection;

		public CalculateWeightTask(STBipartition stb,
				ClusterCollection collection) {
			this.stb = stb;
			this.containedClusterCollection = collection;
		}

		int calculateMissingWeight() {
			// System.err.print("Calculating weight for: " + biggerSTB);
			int weight = 0;
			// System.err.print("Calculating weight for: " + biggerSTB);                              
			BitSet X =  stb.cluster1.getBitSet() ;
			BitSet Y =  stb.cluster2.getBitSet() ;
			//System.out.println("Calculating missing Weight");
			//System.out.println(stb.toString());
			for (STBipartition smallerSTB : clusters.getContainedGeneTreeSTBs()) {
					int temp = 0;
					// possible place of implementation.
				//	System.out.println("Loop "+ smallerSTB.toString() +"count =  "+ geneTreeSTBCount.get(smallerSTB)); 
					BitSet A = smallerSTB.cluster1.getBitSet();
					BitSet B = smallerSTB.cluster2.getBitSet();
					
					BitSet X1 = and(X,A);
					BitSet Y1 = and(Y,B);
					temp = apply(X1, Y1);

					BitSet X2 = and(X,B);
					BitSet Y2 = and(Y,A);
					temp += apply(X2, Y2);
								
					//System.out.println(geneTreeSTBCount.get(smallerSTB));
					temp *= geneTreeSTBCount.get(smallerSTB);
					weight += temp;
					//System.out.println(smallerSTB.toString() + " :: "+ temp);
			}
						
		//	System.out.println("STB score of + " + stb.toString() + " is =  "+weight);
			// System.err.print(" ... " + weight);
			
			if (!rooted) {
				throw new RuntimeException("Unrooted not implemented.");
				/*
				 * for (STBipartition rootSTB : geneTreeRootSTBs.keySet()) { int
				 * c = geneTreeRootSTBs.get(rootSTB); STBipartition inducedSTB =
				 * biggerSTB.getInducedSTB(rootSTB.c); if
				 * (inducedSTB.equals(rootSTB)){ weight -= 2 * c;
				 * //System.err.print(" .. (" + rootSTB +" )" +c+" "+ weight);
				 * if (inducedSTB.cluster1.getClusterSize() != 1 &&
				 * inducedSTB.cluster2.getClusterSize() != 1) { weight -= 2 * c;
				 * //System.err.print(" . " + weight); } } }
				 */
			}
			weights.put(stb, weight);
			// System.err.println("Weight of " + biggerSTB + " is " + weight);
			return weight;
		}

		protected Integer compute() {
			return calculateMissingWeight();
		}

	}
	int apply(BitSet x, BitSet y) {
		int res = 0;
		int c1 = x.cardinality();
		int c2 = y.cardinality();

		res = c1*(c1-1)*c2 + c2*(c2-1)*c1;
		res /=2;
		return res;
	}
	BitSet and(BitSet x, BitSet y) {
		BitSet _x = (BitSet) x.clone();
		BitSet _y = (BitSet) y.clone();
		_x.and(_y);

		return _x;
	}

	/*
	 * public void addGoodSTB (STBipartition good, int size) {
	 * goodSTBs.get(size).add(good); }
	 */
	/*
	 * public Set<STBipartition> getClusterBiPartitions(STITreeCluster cluster)
	 * {
	 * 
	 * return clusterToSTBs.get(cluster); }
	 */

	/*
	 * public boolean addCompleteryVertx(Vertex x, STITreeCluster refCluster) {
	 * STITreeCluster c = x._cluster; Vertex reverse = new Vertex();
	 * reverse._cluster = new STITreeCluster(refCluster);
	 * reverse._cluster.getCluster().xor(c.getCluster()); int size =
	 * reverse._cluster.getClusterSize(); if
	 * (!clusters.get(size).contains(reverse)){ clusters.get(size).add(reverse);
	 * return true; //System.err.println("Clusters: "+clusters); } return false;
	 * }
	 */
	/*
	 * private void addSTBToX(STBipartition stb, int size) {
	 * //System.err.println("Adding to X: "+stb+" "+stb.c
	 * +" "+clusterToSTBs.containsKey(stb.c)); if
	 * (clusterToSTBs.containsKey(stb.c) &&
	 * clusterToSTBs.get(stb.c).contains(stb)){ return; } //int size =
	 * stb.c.getClusterSize(); // TODO: following line is algorithmically
	 * harmless, // but inefficient. is it necessary?
	 * //geneTreeSTBCount.put(stb, 0); addToClusters(stb.c, size, false); //
	 * Following needed for Fast //Set<STBipartition> stbs =
	 * clusterToSTBs.get(c); //stbs = (stbs== null)? new
	 * HashSet<STBipartition>() : stbs; //stbs.add(stb); //clusterToSTBs.put(c,
	 * stbs); //System.err.println("X updated: "+STBCountInGeneTrees);
	 * //System.err.println("X updated: "+clusterToSTBs); }
	 */

	/*
	 * private STITreeCluster treeComplementary(STITreeCluster gtCluster,
	 * String[] leaves){ //System.err.print("Tree complementary of "+gtCluster);
	 * STITreeCluster newGTCluster = gtCluster.complementaryCluster();
	 * //System.err.println(" is: "+newGTCluster.getCluster()); STITreeCluster
	 * newSTCluster = new STITreeCluster(leaves); for (String s :
	 * newGTCluster.getClusterLeaves()) {
	 * newSTCluster.addLeaf(getSpeciesName(s)); }
	 * //System.err.println("Tree complementary of "
	 * +gtCluster+" is: "+newSTCluster); return newSTCluster; }
	 */

	/*
	 * private STITreeCluster treeComplementary(List<String> treeNames, Cluster
	 * c , TaxonNameMap taxonMap){ HashSet<String> set = new HashSet<String> ();
	 * set.add(cluster); return treeComplementary(treeNames, set, taxonMap); }
	 */

	/*
	 * void addExtraBipartitionsByHeuristics(ClusterCollection clusters2) {
	 * //goodSTBs = X; //if (true) return; int added = 0; for (int i=1;
	 * i<goodSTBs.size(); i++) { Set<STBipartition> curr_set = goodSTBs.get(i);
	 * for (STBipartition stb1:curr_set) { //if (Math.random() < 0.70) continue;
	 * for (int j=i; j<goodSTBs.size(); j++) { Set<STBipartition> other_set =
	 * goodSTBs.get(j); //if (Math.random() < 0.70) continue; for (STBipartition
	 * stb2:other_set) { //System.out.println(stb1 +" **AND** " + stb2); if
	 * (stb1.cluster1.getClusterSize() < 3 || stb1.cluster2.getClusterSize() < 3
	 * || stb2.cluster1.getClusterSize() < 3 || stb2.cluster2.getClusterSize() <
	 * 3) { if (tryToAdd(stb1,stb2,bipToAddToX) != null) added++; } }
	 * System.err.println(bipToAddToX.size() + " " + i); } } }
	 * 
	 * for (STBipartition stb: bipToAddToX) { //System.err.println( "Adding: " +
	 * stb); addSTBToX(clusters, stb); } System.out.println("\n\nAdded " +
	 * added+ " bipartitions:\n");
	 * 
	 * int s = 0; for (Integer c: clusters.keySet()){ s +=
	 * clusters.get(c).size(); }
	 * System.out.println("Number of Clusters After Addition: " +s);
	 * 
	 * }
	 * 
	 * 
	 * private STBipartition tryAddingExtraSTB_AndreRule(STBipartition stb1,
	 * STBipartition stb2, Set<STBipartition> bipToAddToX) { if
	 * (stb1.equals(stb2)) return null; if ( stb1.isDominatedBy(stb2) ||
	 * stb2.isDominatedBy(stb1) ) return null;
	 * 
	 * if ( stb1.c.isDisjoint(stb2.c) ) return null;
	 * 
	 * if ( stb1.cluster1.isDisjoint(stb2.cluster2) &&
	 * stb1.cluster2.isDisjoint(stb2.cluster1)) { STITreeCluster cl1 = new
	 * STITreeCluster(stb1.cluster1); cl1 = cl1.merge(stb2.cluster1);
	 * STITreeCluster cl2 = new STITreeCluster(stb1.cluster2); cl2 =
	 * cl2.merge(stb2.cluster2); STITreeCluster cl = new STITreeCluster(stb1.c);
	 * cl = cl.merge(stb2.c); STBipartition r = new STBipartition(cl1,cl2,cl);
	 * bipToAddToX.add(r); return r; } else if (
	 * stb1.cluster1.isDisjoint(stb2.cluster1) &&
	 * stb1.cluster2.isDisjoint(stb2.cluster2) ) { STITreeCluster cl1 = new
	 * STITreeCluster(stb1.cluster1); cl1 = cl1.merge(stb2.cluster2);
	 * STITreeCluster cl2 = new STITreeCluster(stb1.cluster2); cl2 =
	 * cl2.merge(stb2.cluster1); STITreeCluster cl = new STITreeCluster(stb1.c);
	 * cl = cl.merge(stb2.c); STBipartition r = new STBipartition(cl1,cl2,cl);
	 * bipToAddToX.add(r); return r; } return null; }
	 */
}
