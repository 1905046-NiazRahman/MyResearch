package phylonet.coalescent;


import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.util.Collapse;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;

public class MGDInference_DP {
	public static double _versinon =  0.0;
	static boolean _print = true;
	int optimizeDuploss = 1; //one means dup, 3 means duploss
	boolean rooted = true;
	boolean fast = false;
	boolean extrarooted = true;
	double DLbdWeigth;
	double CS;
	double CD;

	List<Tree> trees;
	private List<Tree> extraTrees = null;
	//Map<STITreeCluster, Vertex> clusterToVertex;
	double sigmaNs;
	DuplicationWeightCounter counter;
	TaxonNameMap taxonNameMap = null;
	private boolean exactSolution;
	private STITree st = null;
	public Tree SpeciesTree = null;
	public Boolean STBforST = false ; 
	public TreeAllCluster tas = null; 
	public String[] stTaxaList ; 
	boolean fakeST = true ; 
	
	class TaxonNameMap {
		Map<String, String> taxonMap;
		String pattern = null;
		String rep = null;
		public TaxonNameMap (Map<String, String> taxonMap) {
			this.taxonMap = taxonMap;
		}
		public TaxonNameMap (String pattern, String rep) {
			this.pattern = pattern;
			this.rep = rep;
		}
		public String getTaxonName(String geneName) {
			if (geneName == null || "".equals(geneName)) {
				throw new RuntimeException("Empty name?");
			}
			if (pattern != null) {
				String s = geneName.replaceAll(pattern,rep);
				//System.err.println("Mapped " + geneName + " to " + s);
				return s;
			} else {
				return taxonMap.get(geneName);
			}
		}		
	}

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees,
			Map<String, String> taxonMap) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		if (taxonMap != null) {
			this.taxonNameMap = new TaxonNameMap(taxonMap);
		}
	}

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees,
			String pattern, String rep) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		this.taxonNameMap = new TaxonNameMap (pattern, rep);
	}
	public static final String YELLOW = "\033[33m";
	public static final String BLUE = "\033[34m";
	public static final String MAGENTA = "\033[35m";
	public static final String CYAN = "\033[36m";
	public static final String GREEN = "\u001B[32m";
	public static final String RESET = "\033[0m";
	
	// a method that returns the trees of inference 
	public List<Tree> getTrees() {
		return trees;
	}
	public static void main(String[] args) {

		

        System.err.println("This is STELAR version " + _versinon+"\n");
		System.err.println("========================================================\n");
		
		

		if ((args == null) || args.length == 0 || (args[0].equals("-h"))
				|| (args.length < 1)) {
			printUsage();
			return;
		}
		
		int optimizeDuploss = 0;
		boolean rooted = true;
		boolean fast = false;
		boolean extrarooted = true;
		boolean exactSolution = false;
		
		Map<String, String> taxonMap = null;
		String rep = null;
		String pattern = null;
		List<Tree> trees = null;
		List<Tree> extraTrees = null;
		String output = null;
		String timeFile = null;
		STITree scorest = null;
		// boolean explore = false;
		// double proportion = 0.0D;
		// boolean exhaust = false;
		double bootstrap = 1.0D;
		double cs = 1.0D;
		double cd = 1.0D;
		double time = -1.0D;
		double wd = 1.0D;
		double wh = 1.0D;
		boolean unresolved = false;		
		long startTime = System.currentTimeMillis();
		String line;
		BufferedReader treeBufferReader = null;
		BufferedReader extraTreebuffer = null;
		String inputFilename = null;
		boolean fakeST = false ;
		
		try {
			List<String[]> options = getOptions(args);
			for (String[] option : options) {
				if (option[0].equals("-i")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					inputFilename = option[1];
					System.out.println(YELLOW+"Reading trees from file: "+inputFilename+RESET);
					
					treeBufferReader = new BufferedReader(new FileReader(option[1]));

					trees = new ArrayList();
					// String line;
				} else if (option[0].equals("-ex")) {
					if (option.length != 2) {
						printUsage();
						return;
					}					
					extraTreebuffer = new BufferedReader(new FileReader(
							option[1]));

					extraTrees = new ArrayList();					
				} else if (option[0].equals("-st")) {
					if (option.length != 2) {
						printUsage();
						return;
					}					
					BufferedReader tmp = new BufferedReader(new FileReader(
							option[1]));
					line = tmp.readLine();
					NewickReader nr = new NewickReader(new StringReader(line));
					scorest = new STITree(true);
					nr.readTree(scorest);				
				} else if (option[0].equals("-a")) {
					if ( (option.length != 2) && (option.length != 3)) {
						printUsage();
						return;
					}
					if (option.length == 2) {
						BufferedReader br = new BufferedReader(new FileReader(
								option[1]));
	
						taxonMap = new HashMap<String, String>();
						while ((line = br.readLine()) != null) {
							// String line;
							String[] mapString = line.trim().split(";");
							for (String s : mapString) {
								String species = s.substring(0, s.indexOf(":"))
										.trim();
								s = s.substring(s.indexOf(":") + 1);
								String[] alleles = s.split(",");
								for (String allele : alleles) {
									allele = allele.trim();
									if (taxonMap.containsKey(allele)) {
										System.err
												.println("The input file is not in correct format");
										System.err
												.println("Any gene name can only map to one species");
										System.exit(-1);
									} else {
										taxonMap.put(allele, species);
									}
								}
							}
						}
						br.close();
					} else {
						pattern = option[1];
						rep = option [2];
						if (rep.equals("/delete/")) {
							rep = "";
						}
					}
				} else if (option[0].equals("-o")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					output = option[1];
					setPrinting (false);
				} else if (option[0].equals("-x")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					// exhaust = true;
				} else if (option[0].equals("-e")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					// explore = true;
					if (option.length != 2)
						continue;
					try {
						// proportion = Double.parseDouble(option[1]);
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-b")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					if (option.length != 2)
						continue;
					try {
						bootstrap = Double.parseDouble(option[1]);
						if ((bootstrap <= 1.0D) && (bootstrap > 0.0D))
							continue;
						printUsage();
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-cs")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					try {
						cs = Double.parseDouble(option[1]);
						if ((cs <= 1.0D) && (cs >= 0.0D))
							continue;
						printUsage();
						return;
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-cd")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					try {
						cd = Double.parseDouble(option[1]);
						if ((cd <= 1.0D) && (cd >= 0.0D))
							continue;
						printUsage();
						return;
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-ur")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					unresolved = true;
				} else if (option[0].equals("-f")) {
					if ((option.length != 1)) {
						printUsage();
						return;
					}
					fast = true;
				} else if (option[0].equals("-u")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					rooted = false;
				} else if (option[0].equals("-xu")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					extrarooted = false;
				} else if (option[0].equals("-t")) {
					if ((option.length != 2) || (unresolved)) {
						printUsage();
						return;
					}
					if (option.length != 2)
						continue;
					try {
						time = Double.parseDouble(option[1]);
						if (time > 0.0D)
							continue;
						printUsage();
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}
				} else if (option[0].equals("-xt")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					exactSolution = true;
				} else if (option[0].equals("-dl")) {
					optimizeDuploss = 1;
					if (option.length != 2) {
						printUsage();
						return;
					}					
					try {
						if (option[1].equals("auto")) {
							wh = -1;
							continue;
						} else {
							wh = Double.parseDouble(option[1]);
							if (wh >= 0.0D)
								continue;
						}
						printUsage();
						return;
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter wd");
						printUsage();
						return;
					}
				} 
				else if(option[0].equals("-t")) {
					timeFile = option[1];
				}
				else {
					printUsage();
					return;
				}
			} //! finish reading options from command line arguments

			if (treeBufferReader == null) {
				System.err.println("The input file has not been specified.");
				printUsage();
				return;
			}

			System.err.println("Gene trees are treated as "
					+ (rooted ? "rooted" : "unrooted"));
			int l = 0;
			try {
				while ((line = treeBufferReader.readLine()) != null) {
					l++;
					// previousTreeTaxa keeps track of the taxa in the previous trees . 
					Set<String> previousTreeTaxa = new HashSet<String>();
					if (line.length() > 0) {
						NewickReader nr = new NewickReader(new StringReader(line));
						if (rooted) {
							STITree gt = new STITree(true); // true means rooted tree
							nr.readTree(gt); // read the tree and store it in gt 
							if (previousTreeTaxa.isEmpty()) 
							{
								// if the previous tree is empty, add all the leaves of the current tree to the previous tree taxa
								previousTreeTaxa.addAll(Arrays.asList(gt
										.getLeaves()));
							} else {
								// if the previous tree is not empty, check if the leaves of the current tree are the same as the previous tree
								if (!previousTreeTaxa.containsAll(Arrays.asList(gt
										.getLeaves()))) {
									throw new RuntimeException(
											"Not all trees are on the same set of taxa: "
													+ gt.getLeaves() + "\n"
													+ previousTreeTaxa);
								}
							}
							trees.add(gt);
						} else {						
							Tree tr = nr.readTree();
							trees.add(tr);
						}
					}
				}
				treeBufferReader.close();
			} catch (ParseException e) {
				treeBufferReader.close();
				throw new RuntimeException("Failed to Parse Tree number: " + l ,e);
			}			

			if (extraTreebuffer != null) {
				while ((line = extraTreebuffer.readLine()) != null) {
					if (line.length() > 0) {					
						NewickReader nr = new NewickReader(
								new StringReader(line));
						if (extrarooted) {
							STITree gt = new STITree(true);
							nr.readTree(gt);
							extraTrees.add(gt);
						} else {
							Tree tr = nr.readTree();
							extraTrees.add(tr);
						}
					}
				}				
				extraTreebuffer.close();
			}
		} catch (IOException e) {
			System.err.println("Error when reading trees. The function exits.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		} catch (ParseException e) {
			System.err
					.println("Error when parsing the Newick representation from input file.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		}

		if (trees.size() == 0) {
			System.err.println("Empty list of trees. The function exits.");
			return;
		}

		if (_print) {
			System.err.println("Reading trees in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		startTime = System.currentTimeMillis();
		System.out.println("List of trees : ");
		for (Tree tr : trees) {
			System.out.println(BLUE+tr+RESET);
		}

		
		
		MGDInference_DP inference;		
		
		if (rep != null) {			
			inference = new MGDInference_DP(trees, extraTrees,
				pattern, rep);
		} else if (taxonMap != null) {
			inference = new MGDInference_DP(trees, extraTrees,
				taxonMap);
		} else {
			inference = new MGDInference_DP(trees, extraTrees, null); // create a new instance of MGDInference_DP
		}
		
		inference.optimizeDuploss = optimizeDuploss > 0 ? 3 : 1;
		inference.DLbdWeigth = wh; 
		inference.rooted = rooted;
		inference.fast = fast;
		inference.extrarooted = extrarooted;
		inference.CS = cs;
		inference.CD = cd;
		inference.exactSolution = exactSolution;
		inference.st = scorest;

		if (scorest != null) {
			//System.out.println("in score gene tree");
			inference.scoreGeneTree();
			System.exit(0);
		}
		else {
			//System.out.println("not in score gene tree");
		}
		startTime = System.currentTimeMillis();
		List<Solution> solutions = inference.inferSpeciesTree(); //TODO : check this function which returns a list of solutions 
		long endTime = System.currentTimeMillis();


        System.err.println("Optimal tree inferred in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");

        if(inputFilename !=null){
            try {
                PrintWriter pw = new PrintWriter(new FileWriter(inputFilename+".log", true));
                pw.append("Optimal tree inferred in "
                        + (System.currentTimeMillis() - startTime) / 1000.0D
                        + " secs for "+inputFilename+"\n");
                pw.close();
            }catch(Exception ex){
                System.out.println("error while outputting time value");
            }

        }

		
		//! printing all possible node data of the species tree
		
		System.out.println(MAGENTA+"Showing all the subtree bipartitions"+RESET);
		
		for (Solution s : solutions) {
				inference.SpeciesTree = s._st; 
				inference.STBforST = true;
				List<Solution> x = inference.inferSpeciesTree();
				inference.STBforST  = false;
				for(TNode node: s._st.postTraverse()) {
					STINode stnode = (STINode) node;
					//System.out.println(GREEN+stnode.getData()+RESET);
				}
			
		}



		//! FOR TESTING PURPOSE , READ FROM A FILE HERE and update the inference.SpeciesTree 
		System.out.println(inference.SpeciesTree);
		if(fakeST) {
			try {
				System.out.println(CYAN+"Fake ST"+RESET); ;
				BufferedReader br = new BufferedReader(new FileReader("E:/Thesis/STELAR-master2/STELAR-master/main/Fake_Data/fakeST.tre"));
				line = br.readLine();
				System.out.print(line);
				NewickReader nr = new NewickReader(new StringReader(line));
				STITree st = new STITree(true);
				nr.readTree(st);
				inference.SpeciesTree = st;
				System.out.println("HEREEEEEEEEEEEEEEEEE") ; 
				System.out.println(YELLOW+inference.SpeciesTree+RESET);
				br.close();
			}catch(Exception e) {
				System.out.println("Error while reading species tree");
			}
		}
		 
		inference.computeSTtaxalist() ; 
		System.out.println(YELLOW+"Computing treeAllClusters for all the gene trees"+RESET);
		inference.tas = inference.computeallclusters() ; 
		System.out.println(BLUE+"Showing all the tripartitions of the ST"+RESET);
		for (Solution s : solutions)
		{
			//inference.SpeciesTree = s._st ;
			inference.computeTripartitions() ;  
		}

		System.out.println(YELLOW+"Inferred species tree: "+RESET);

		if ((_print)) {
			String metric;
			if (optimizeDuploss == 0) {
				metric = "duplications";
			} else if (optimizeDuploss == 1) {
				metric = "duplication+loss (homomorphic)";
			} else {
				metric = "duplication+loss (original)";
			}
			for (Solution s : solutions) {
				
				System.out.println(
					        s._st.toStringWD()
					//	+ " \n"
						//+ s._totalCoals
						//+ " " + metric + " in total");
					        );
			}
		} else
			try {
				FileWriter fw = new FileWriter(output);
			//	System.out.println(solutions.size());

				// for (Solution s : solutions) {
				// 	fw.write(s._st.toString()+ " \n");
				//}
				System.out.println(inference.SpeciesTree.toStringWD());
				fw.write(inference.SpeciesTree.toStringWD());
				fw.close();
			} catch (IOException e) {
				System.err.println("Error when writing the species tree");
				System.err.println(e.getMessage());
				e.printStackTrace();
			}
			

		

		
		
		//	*/
	}

	


	private void computeSTtaxalist() {
		String[] gtTaxa ; 
		String[] stTaxa ;
		Collapse.CollapseDescriptor cd = null;
		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}
		List<String> taxalist = new ArrayList<String>();
		for (Tree tr : trees) {
			for (TNode node : tr.postTraverse()) {
				if (node.isLeaf()) {
					taxalist.add(node.getName());
				}
			}
		}
		stTaxa = new String[taxalist.size()];
		int index = 0;
		for (String taxon : taxalist) {
			stTaxa[(index++)] = taxon;
		}
		gtTaxa = stTaxa;
		stTaxaList = stTaxa ;
		System.out.println("Taxa list : " + Arrays.toString(stTaxaList));
	}

	private TreeAllCluster computeallclusters() {
	
		TreeAllCluster treeallcluster = new TreeAllCluster() ; 
		int noOfGeneTrees = treeallcluster.computeAllClusters(this) ; 
		// print the tree all cluster 
		treeallcluster.printtreeallcluster();

		treeallcluster.setupgeneTreesAsInts(this) ;
		// PRINT THE WHOLE GENETREE AS INT 
		treeallcluster.printGeneTreesAsInts() ;
		treeallcluster.printGeneTreesAsString();
		return treeallcluster ;
	}

	private void computeTripartitions() {
		System.out.println("Computing tripartitions");

		ClusterCollection clusters = new BasicClusterCollection(stTaxaList.length);
		counter = new DuplicationWeightCounter(stTaxaList, stTaxaList, rooted,taxonNameMap, clusters);

		// set the data for the species tree
		


		double retval = counter.computeTripartitionsforST(this) ;


	
		
		


	}

	private int [] calc(Tree gtTree, SchieberVishkinLCA lcaLookup, Tree stTree) { // lcalookup is made from stTree
		int [] res = {0,0,0};
	//	Stack<TNode> stack = new Stack<TNode>();
		for(TNode t: gtTree.getNodes()) {
			if(t.isLeaf()) {
				continue;
			}
			
			List<TNode> gtLeaf = new ArrayList<TNode>();
			for(TNode child: t.getChildren()) {
		//		System.out.print(child.getLeaves());
				gtLeaf.add(child);
			}
	//		System.out.println("\n====");
			for(TNode t2: stTree.getNodes()) {
				List<TNode> stLeaf = new ArrayList<TNode>();
				if(t2.isLeaf()) {
					continue;
				}
				for(TNode child: t2.getChildren()) {
	//					System.out.print(child.getLeaves());
						stLeaf.add(child);
				}
	//			System.out.println("");
				int score=0;
				//System.out.println(stLeaf.size()+" "+gtLeaf.size());
				score += d (stLeaf.get(0),gtLeaf.get(0), stLeaf.get(1),gtLeaf.get(1));
				//	System.out.print(score+" ");
				score += d (stLeaf.get(1),gtLeaf.get(0), stLeaf.get(0),gtLeaf.get(1));
			//	System.out.println(score);
				res[0]+=score;
			}
	//		System.out.println("");
		//	TNode lca = lcaLookup.getLCA(gtLeaf.get(0),gtLeaf.get(1));
		//	if(lca == null || lca.isLeaf()) {
		//		continue;
		//	}
		//	for(TNode t1: lca.getChildren()) {
		//		stLeaf.add(t1);
		//	}
			
		}
		/*
		for (TNode gtNode : gtTree.postTraverse()) {
			if (gtNode.isLeaf()) {
			    	TNode node = stTree.getNode(this.taxonNameMap !=null ? 
					this.taxonNameMap.getTaxonName(gtNode.getName()):
						gtNode.getName());
			    	if (node == null) {
					throw new RuntimeException("Leaf " + gtNode.getName() +
						" was not found in species tree; mapped as: "+
						this.taxonNameMap.getTaxonName(gtNode.getName())); 
			    	}
			    	stack.push(node);
				//System.out.println("stack: " +this.taxonNameMap.getTaxonName(gtNode.getName()));
			} else {
				
				TNode leftLCA = stack.pop();
				TNode rightLCA = stack.pop();
				if (rightLCA == null || leftLCA == null) {
					System.out.println("Should not be printed");
					stack.push(null);
					continue;
				}
				TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
			//	System.out.println(rightLCA.getParent()+" is mapped to "+lca);
				List<TNode> stLeaf = new ArrayList<TNode>();
				List<TNode> gtLeaf = new ArrayList<TNode>();
				
				gtLeaf.add(leftLCA);
				gtLeaf.add(rightLCA);
				
				for(TNode node: lca.getChildren()) {
					stLeaf.add(node);
				}
				if(stLeaf.size()!=2) {
					System.out.println("Some thing is wrong!" + stLeaf.size());
					return res;
				}

				for(TNode t: stLeaf) {
					System.out.print(t.getLeaves());
				}
				System.out.print(" ");
				for(TNode t: gtLeaf) {
					System.out.print(t.getLeaves());
				}
				System.out.print(" ");
				
				int score=0;
				score += d (stLeaf.get(0),gtLeaf.get(0), stLeaf.get(1),gtLeaf.get(1));
				System.out.print(score+" ");
				score += d (stLeaf.get(1),gtLeaf.get(0), stLeaf.get(0),gtLeaf.get(1));
				System.out.println(score);
				res[0]+=score;
				/*
				for(TNode node: rightLCA.getChildren()){
					System.out.print(node.getID()+" ");
				}
				System.out.println("");
				
				// If gene trees are incomplete, we can have this case
				if (rightLCA == null || leftLCA == null) {
					stack.push(null);
					continue;
				}
				TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
				*/
		//		if(leftLCA.getParent()!=rightLCA.getParent()) {
		//			System.out.println("!OK");
		//		}
		//		stack.push(lca);
				//break;
				/*
				if (lca == leftLCA || lca == rightLCA) {
					// LCA in stTree dominates gtNode in gene tree
					res[0]++;
					if (lca == leftLCA && lca == rightLCA) {
						res[1] += 0;
					} else {
						res[1] += (lca == leftLCA) ?
									d(rightLCA,lca) + 1:
									d(leftLCA,lca) + 1;
					}
				} else {
					res[1] += (d(rightLCA,lca) + d(leftLCA,lca));
				}
			}
		}
	//	TNode rootLCA = stack.pop();
	//	res[2] = res[1];
		res[0] += d(rootLCA,stTree.getRoot()) + (rootLCA == stTree.getRoot()?0:1);
		*/
		return res;
	}
	
	public int d(TNode X, TNode A, TNode Y, TNode B) {
		int d1 = intersectionCount(X,A);
		int d2= intersectionCount(Y,B);
		return F(d1,d2) + F(d2,d1);
	}
	int F(int p,int q) {
		return (p*q*(q-1))/2;
	}

	private int intersectionCount(TNode n1, TNode n2) {
		int val = 0;
		for(TNode node1: n1.getLeaves()) {
			for(TNode node2: n2.getLeaves()) {
			//	System.out.println(node1.getName()+ "&"+node2.getName());
				if(node1.getName().equals(node2.getName())) {
					val++;
				}
			}
		}
	//	System.out.println(n1.getLeaves()+" & "+n2.getLeaves()+" & score is = "+val);
		return val;
	}

	private void scoreGeneTree() {
		
		// does this calculate st score when a st is given -st argument?
		// first calculated duplication cost by looking at gene trees. 
		
		SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(this.st); // lcaLookUp is a tree
	//	System.out.println(lcaLookup.getTree());
		Integer duplications = 0;
		Integer triplets = 0;
		Integer lossesstd = 0;
		
		for (Tree gtTree : this.trees) {
			int[] res = calc(gtTree,lcaLookup, this.st); // calc what? why we are sending lcaklookup and st since lcalookup is made from st?
			triplets += res[0];
		//	lossesstd += res[1];
			
		//	STITree hmst = new STITree(this.st);
			//hmst.constrainByLeaves(Arrays.asList(gtTree.getLeaves()));
		//	SchieberVishkinLCA hmlcaLookup = new SchieberVishkinLCA(hmst);
		//	int[] res2 = calc(gtTree,hmlcaLookup, hmst);
			
		//	lossesstd += res2[2];
		}
	//	System.out.println("Total number of duplications is: "+duplications);
	//	System.out.println("Total number of losses (bd) is: "+losses);
	//	System.out.println("Total number of losses (std) is: "+lossesstd);
	//	System.out.println("Total number of duploss (bd) is: " + (losses+duplications));
	//	System.out.println("Total number of duploss (st) is: " + (lossesstd+duplications));
	//	System.out.println("Total weighted (wd = "+this.DLbdWeigth+") loss is: " + (lossesstd + this.DLbdWeigth*(losses-lossesstd)));
		System.out.println("Total number of triplets satisfied in (st) is: " + triplets);
	}

	private int d (TNode down, TNode upp) {
		int ret = 0;
		TNode t = down;
		//System.err.println("Down: "+down+"\nUPP: "+upp);
		while (t != upp) {ret++; t=t.getParent();}
		return Math.max(ret-1,0);
	}
	protected static void printUsage() {
		System.out
				.println("This tool infers the species tree from rooted gene trees  by maximizing  the triplet distribution in the input gene trees.");
		System.out.println("Usage is:");
		System.out
				.println("\tMGDInference_DP -i input [-a mapping]  [-ex extra_trees] [-o output]  [-xt] [-s species tree] [-q score_tree ]");
		System.out
				.println("\t-i gene tree file: The file containing gene trees. (required)");
		System.out
				.println("\t-st species tree file: The file containing a species tree to be scored.\n" +
						 "\t                       If this option is provided the software only scores the species tree.");
		System.out
				.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled.\n" +
						 "\t                 Alternatively, two reqular expressions for automatic name conversion (optional)");
		System.out
				.println("\t-o species tree file: The file to store the species tree. (optional)");
		//System.out.println("\t-dl N: optimize duplications and losses. Use -dl 0 for standard (homomorphic) definition, and -dl 1 for ``bd'' definition. Any value in between weights the impact of missing taxa on the tree.");
		System.out.println("\t-xt find the exact solution by looking at all clusters.");
		//System.out.println("\t-u treat input gene trees as unrooted (Not implemented!)");
		System.out.println("\t-ex provide extra trees to add to set of STBs searched");
		//System.out.println("\t-xu treat extra trees input gene trees as unrooted (Not implemented!)");
		System.out.println("\t-cs and " +
						   "-cd thes two options set two parameters (cs and cd) to a value between 0 and 1. \n" +
						   "\t    For any cluster C if |C| >= cs*|taxa|, we add complementary clusters (with respect to C) of all subclusters of C\n" +
						   "\t    if size of the subcluster is >= cd*|C|.\n" +
						   "\t    By default cs = cd = 1; so no extra clusters are added. Lower cs and cd values could result in better scores\n" +
						   "\t    (especially when gene trees have low taxon occupancy) but can also increase the running time dramatically.");
		
		//System.out.println("\t-f perform fast and less-accurate subtree-bipartition based search (Not implemented!).");
		System.out.println();
	}

	public static void setPrinting(boolean print) {
		_print = print;
	}

	public List<Solution> inferSpeciesTree() {
		System.out.println(YELLOW+"Infering species tree"+RESET);
		long startTime = System.currentTimeMillis();

		String[] gtTaxa;
		String[] stTaxa;
		Collapse.CollapseDescriptor cd = null;

		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}
		if (taxonNameMap != null && taxonNameMap.taxonMap != null) {
			System.out.println(MAGENTA+"Using taxon mapping when both not null"+RESET);
			System.out.println(taxonNameMap.taxonMap);

			Map<String,String> taxonMap = taxonNameMap.taxonMap;
			String error = Trees.checkMapping(trees, taxonMap);
			if (error != null) {
				throw new RuntimeException("Gene trees have a leaf named "
						+ error
						+ " that hasn't been defined in the mapping file");
			}

			List temp1 = new LinkedList();
			List temp2 = new LinkedList();
			for (String s : taxonMap.keySet()) {
				temp1.add(s);
				if (!((List) temp2).contains(taxonMap.get(s))) {
					((List) temp2).add((String) taxonMap.get(s));
				}
			}
			gtTaxa = new String[temp1.size()];
			stTaxa = new String[temp2.size()];

			for (int i = 0; i < gtTaxa.length; i++) {
				gtTaxa[i] = ((String) temp1.get(i));
			}
			for (int i = 0; i < stTaxa.length; i++) {
				stTaxa[i] = ((String) ((List) temp2).get(i));
			}
		} else if (taxonNameMap != null && taxonNameMap.taxonMap == null) {
			System.out.println(MAGENTA+"Using taxon mapping when taxonMap is null"+RESET);

			Set<String> taxalist = new HashSet<String>();
			Set<String> genelist = new HashSet<String>();
			for (Tree tr : trees) {
				String[] leaves = tr.getLeaves();
				for (int i = 0; i < leaves.length; i++) {
					String leaf = leaves[i];				
					genelist.add(leaf);
					taxalist.add(taxonNameMap.getTaxonName(leaf));
				}
			}			

			stTaxa = new String[taxalist.size()];
			gtTaxa = new String[genelist.size()];

			int index = 0;
			for (String taxon : taxalist) {
				stTaxa[(index++)] = taxon;
			}
			index = 0;
			for (String gene : genelist) {
				gtTaxa[(index++)] = gene;
			}
		} else {
			System.out.println(MAGENTA+"Using taxon mapping when both null"+RESET);
			cd = null;
			if (rooted & extraTrees == null & taxonNameMap == null && false) {
				cd = doCollapse(trees);
			}

			List<String> taxalist = new ArrayList<String>();
			for (Tree tr : trees) {
				for (TNode node : tr.postTraverse()) {
					// System.out.println(node.getName());
					// System.out.println(node.getParent());
					// System.out.println("");
					// System.out.println("");
					if ((node.isLeaf()) && (!taxalist.contains(node.getName()))) {
						taxalist.add(node.getName());
					}
				}
			}
			System.out.println(YELLOW+"Taxa list: "+taxalist+RESET); //* taxalist is the list of taxa in the gene trees

			stTaxa = new String[taxalist.size()];

			int index = 0;
			for (String taxon : taxalist) {
				stTaxa[(index++)] = taxon;
			}
			gtTaxa = stTaxa;
		}

		System.err.println("Number of taxa: " + stTaxa.length);
		System.err.println("Taxa: " + Arrays.toString(stTaxa));

		ClusterCollection clusters = new BasicClusterCollection(stTaxa.length); // create a new instance of BasicClusterCollection , with the number of taxa in the species tree 
		

		List<Solution> solutions;
		// print all of the arguments of the function
		System.out.println("gtTaxa: "+Arrays.toString(gtTaxa)); //* the taxa in the gene trees
		System.out.println("stTaxa: "+Arrays.toString(stTaxa)); //* the taxa in the species tree
		System.out.println("rooted: "+rooted); //* true
		System.out.println("clusters: "+clusters); //* BasicClusterCollection
		counter = new DuplicationWeightCounter(gtTaxa, stTaxa, rooted,taxonNameMap, clusters);

		

		double sigmaN = counter.computeTreeSTBipartitions(this); 

		if(STBforST == true) {
			return null;
		}
		

		if (extraTrees != null) {		
			counter.addExtraBipartitionsByInput(clusters, extraTrees,extrarooted);					
		}
		
		if (exactSolution) {
			counter.addAllPossibleSubClusters(clusters.getTopVertex().getCluster(),  stTaxa.length);
		}

		//counter.addExtraBipartitionsByHeuristics(clusters);

		if (_print) {
			System.err.println("STBs formed in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}
		

		counter.preCalculateWeights(trees, extraTrees);
		
		if (_print) {
			System.err.println("DP starting after "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		sigmaNs = sigmaN ;
		
		//printing stTaxa , counter , trees , taxonNameMap , clusters if not null
		//System.out.println("trees: "+trees); //* the gene trees , same as before 
		//System.out.println("taxonNameMap: "+taxonNameMap); //* null
		ClusterPrinter.printClusters(clusters);


		System.out.println(MAGENTA+"Showing all the subtree bipartitions"+RESET);
		for(Tree t : trees) {
			for(TNode node: t.postTraverse()) {
				STINode stnode = (STINode) node;
				System.out.println(GREEN+stnode.getData()+RESET);
			}
		}


		solutions = findTreesByDP(stTaxa, counter, trees, taxonNameMap,clusters);
		System.out.println("got the solutions");

		if (taxonNameMap == null && rooted && extraTrees == null && false) {
			restoreCollapse(solutions, cd);
		}

		return (List<Solution>) solutions;
	}

	public class ClusterPrinter {
		public static void printClusters(ClusterCollection collection) {
			Iterable<Set<Vertex>> subClusters = collection.getSubClusters();
			for (Set<Vertex> cluster : subClusters) {
				System.out.println("Cluster:");
				for (Vertex vertex : cluster) {
					System.out.println("  " + vertex);
				}
				System.out.println(); // Blank line between clusters
			}
		}
	}

	


	private Collapse.CollapseDescriptor doCollapse(List<Tree> trees) {
		Collapse.CollapseDescriptor cd = Collapse.collapse(trees);
		return cd;
	}

	private void restoreCollapse(List<Solution> sols,
			Collapse.CollapseDescriptor cd) {
		for (Solution sol : sols) {
			Tree tr = sol._st;
			Collapse.expand(cd, (MutableTree) tr);
			for (TNode node : tr.postTraverse())
				if (((STINode) node).getData() == null)
					((STINode) node).setData(Integer.valueOf(0));
		}
	}

	  public static Tree buildTreeFromClusters(List<STITreeCluster> clusters)
	  {
	    if ((clusters == null) || (clusters.size() == 0)) {
	      System.err.println("Empty list of clusters. The function returns a null tree.");
	      return null;
	    }

	    MutableTree tree = new STITree();

	    String[] taxa = ((STITreeCluster)clusters.get(0)).getTaxa();
	    for (int i = 0; i < taxa.length; i++) {
	      tree.getRoot().createChild(taxa[i]);
	    }

	    for (STITreeCluster tc : clusters) {
	      if ((tc.getClusterSize() <= 1) || (tc.getClusterSize() == tc.getTaxa().length))
	      {
	        continue;
	      }

	      Set clusterLeaves = new HashSet();
	      TNode node;
	      for (String l : tc.getClusterLeaves()) {
	        node = tree.getNode(l);
	        clusterLeaves.add(node);
	      }

	      SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(tree);
	      TNode lca = lcaFinder.getLCA(clusterLeaves);

	      Object movedChildren = new LinkedList();
	      for (TNode child : lca.getChildren()) {
	        BitSet childCluster = new BitSet(taxa.length);
	        for (TNode cl : child.getLeaves()) {
	          for (int i = 0; i < taxa.length; i++) {
	            if (taxa[i].equals(cl.getName())) {
	              childCluster.set(i);
	              break;
	            }
	          }
	        }

	        BitSet temp = (BitSet)childCluster.clone();
	        temp.and(tc.getBitSet());
	        if (temp.equals(childCluster)) {
	          ((List)movedChildren).add(child);
	        }

	      }

	      STINode newChild = ((STINode)lca).createChild();

	      while (!((List)movedChildren).isEmpty()) {
	        newChild.adoptChild((TMutableNode)((List)movedChildren).get(0));
	        ((List)movedChildren).remove(0);
	      }
	    }

	    return (Tree)tree;
	  }
	  
	private List<Solution> findTreesByDP(String[] stTaxa,
			DuplicationWeightCounter counter, List<Tree> trees,
			TaxonNameMap taxonNameMap, ClusterCollection clusters) {
		List<Solution> solutions = new ArrayList<Solution>();

		

		Vertex all = (Vertex) clusters.getTopVertex();
		System.err.println("Sigma N: " + sigmaNs);

		System.err.println("Size of largest cluster: " +all.getCluster().getClusterSize());

		try {
			//vertexStack.push(all);
			ComputeMinCostTask allTask = new ComputeMinCostTask(this,all,clusters);
			//ForkJoinPool pool = new ForkJoinPool(1);
			allTask.compute();
			double v = all._max_score;
			if (v == Integer.MIN_VALUE) {
				throw new CannotResolveException(all.getCluster().toString());
			}
		} catch (CannotResolveException e) {
			System.err.println("Was not able to build a fully resolved tree. Not" +
					"enough STBs present in input gene trees ");
			e.printStackTrace();
			System.exit(1);
		}

		if (_print) {
			//System.err.println("Weights are: "
				//	+ counter.weights);
		}
		//System.out.println("domination calcs:" + counter.cnt);

		List minClusters = new LinkedList();
		List coals = new LinkedList();
		Stack minVertices = new Stack();
		if (all._min_rc != null) {
			minVertices.push(all._min_rc);
		}
		if (all._min_lc != null) {
			minVertices.push(all._min_lc);
		}
		if (all._subcl != null) {
			for (Vertex v : all._subcl) {
				minVertices.push(v);
			}
		}		
		while (!minVertices.isEmpty()) {
			Vertex pe = (Vertex) minVertices.pop();
			//System.out.println(pe._min_rc);
			//System.out.println(pe._min_lc);
			minClusters.add(pe.getCluster());
			// int k = sigmaNs/(stTaxa.length-1);

			if (pe._min_rc != null) {
				minVertices.push(pe._min_rc);
			}
			if (pe._min_lc != null) {
				minVertices.push(pe._min_lc);
			}
			if (pe._min_lc != null && pe._min_rc != null) {
				coals.add(pe._c);
			} else {
				coals.add(0D);
			}
			if (pe._subcl != null) {
				for (Vertex v : pe._subcl) {
					minVertices.push(v);
				}
			}
		}
		Solution sol = new Solution();
		if ((minClusters == null) || (minClusters.isEmpty())) {
			System.err.println("WARN: empty minClusters set.");
			Object tr = new STITree();
			for (String s : stTaxa) {
				((MutableTree) tr).getRoot().createChild(s);
			}
			sol._st = ((Tree) tr);
		} else {
			sol._st = buildTreeFromClusters(minClusters);
		}
		//System.err.println("SOL: " + sol._st);
		//System.err.println("coals: " + coals);
		//System.err.println("min cluster: " + minClusters);
		Object map = new HashMap();
		for (TNode node : sol._st.postTraverse()) {
			BitSet bs = new BitSet(stTaxa.length);
			if (node.isLeaf()) {
				for (int i = 0; i < stTaxa.length; i++) {
					if (node.getName().equals(stTaxa[i])) {
						bs.set(i);
						break;
					}
				}
				((Map) map).put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = (BitSet) ((Map) map).get(child);
					bs.or(childCluster);
				}
				((Map) map).put(node, bs);
			}
//            System.err.println("Node: "+node);
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.setCluster(bs);
//            System.err.println("m[0]: "+((STITreeCluster)minClusters.get(0)).toString2());
//            System.err.println("C: "+c.toString2());
//            System.err.println("Equals: "+((STITreeCluster)minClusters.get(0)).equals(c));
			if (c.getClusterSize() == stTaxa.length) {
				((STINode) node).setData(Double.valueOf(0));
			} else {
				int pos = minClusters.indexOf(c);                                
				((STINode) node).setData((Double) coals.get(pos));
			}
		}

		//sol._totalCoals = (int) (sigmaNs - all._max_score);
		//System.out.println("total cost: " + (sigmaNs - all._max_score));
	//	System.out.println("Score of constructed species tree = "+all._max_score);
		solutions.add(sol);

		return (List<Solution>) (List<Solution>) solutions;
	}

	protected static List<String[]> getOptions(String[] args) {
		LinkedList opts = new LinkedList();
		LinkedList arg_list = new LinkedList();

		int i = 0;
		while (i < args.length) {
			if (args[i].charAt(0) != '-') {
				printUsage();
				System.exit(-1);
			}

			arg_list.clear();
			arg_list.addFirst(args[i]);
			i++;

			while ((i < args.length) && (args[i].charAt(0) != '-')) {
				arg_list.addLast(args[i]);
				i++;
			}
			String[] arg_array = new String[arg_list.size()];
			arg_list.toArray(arg_array);

			opts.addLast(arg_array);
		}
		return opts;
	}

	private int getResolutionsNumber(int nodeNumber) {
		int total = 1;
		for (int i = 3; i <= nodeNumber; i++) {
			total *= (2 * i - 3);
		}
		return total;
	}

}
