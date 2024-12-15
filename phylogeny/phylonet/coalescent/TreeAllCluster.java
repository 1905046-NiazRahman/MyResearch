package phylonet.coalescent;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import java.util.List;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
public class TreeAllCluster {
    List<STITreeCluster> treeAllClusters = new ArrayList<STITreeCluster>();
    Integer[] geneTreesAsInts = null;
    String[] geneTreesAsString = null;


    public int computeAllClusters(MGDInference_DP mgdInference_DP)
    {
        System.out.println("No. of gene trees : " + mgdInference_DP.getTrees().size());


        for (Tree tr : mgdInference_DP.getTrees())
        {
            String[] leaves = tr.getLeaves();
            STITreeCluster cluster = new STITreeCluster(leaves);
            for (String leaf : leaves)
            {
                cluster.addLeaf(leaf);
            }
            treeAllClusters.add(cluster);

        }
        return mgdInference_DP.getTrees().size();
    }

    public void setupgeneTreesAsInts(MGDInference_DP mgdInference_DP)
    {
        List<Integer> temp = new ArrayList<Integer>();
        List<String> temp1 = new ArrayList<String>();
        for (Tree tr : mgdInference_DP.getTrees())
        {
            for(TNode node : tr.postTraverse())
            {
                if(node.isLeaf())
                {
                    temp.add(new STITreeCluster(mgdInference_DP.stTaxaList).addLeaf(node.getName())) ; 
                    temp1.add(node.getName());
                }
                else {
                    temp.add(-node.getChildCount());

                }
                if(node.isRoot()){
                    temp.add(Integer.MIN_VALUE);
                }

            }

        }
        geneTreesAsInts = temp.toArray(new Integer[] {});
        geneTreesAsString = temp1.toArray(new String[] {});
    }

    public void printGeneTreesAsInts() {
        for (Integer i : geneTreesAsInts) {
            System.out.print(i + " ");
        }
        System.out.println();
    }

    public void printGeneTreesAsString() {
        for (String i : geneTreesAsString) {
            System.out.print(i + " ");
        }
        System.out.println();
    }

    public void printtreeallcluster() {
        for (STITreeCluster cluster : treeAllClusters) {
            System.out.println(cluster.toString());
        }
    }


    class FrequencyData
    {
        double [] freq_s;
		int effn;
		
		FrequencyData (double[] weight, int n){
			freq_s = weight;
			effn = n;
		} 
    }


    private long F(long a,long b,long c) {
		if (a<0 || b<0 || c<0) {
			throw new RuntimeException("negative side not expected: "+a+" "+b+" "+c);
		}
		return a*b*c;
	}	

    class Intersects {
		long s0;
		long s1;
		long s2;



		public Intersects(long s0, long s1, long s2) {
			this.s0 = s0;
			this.s1 = s1;
			this.s2 = s2;
	
		}
		
        public Intersects(Intersects side1, Intersects side2) {
            this(side1.s0+side2.s0,
					side1.s1+side2.s1,
					side1.s2+side2.s2);               
        }	

		public Intersects(Intersects other) {
			this.s0 = other.s0;
			this.s1 = other.s1;
			this.s2 = other.s2;
		
		}

		public void addin(Intersects pop) {
			this.s0 += pop.s0;
			this.s1 += pop.s1;
			this.s2 += pop.s2;  
		
		}

		public void subtract(Intersects pop) {
			this.s0 -= pop.s0;
			this.s1 -= pop.s1;
			this.s2 -= pop.s2;         
			
		}

		public String toString() {
			return this.s0+","+this.s1+"|"+this.s2;
		}
		
        public boolean isNotEmpty() {
        	return (this.s0 + this.s1 + this.s2) != 0;
        }
        
        public boolean hasEmpty() {
        	return this.maxPossible() == 0;
        }
        
        public long maxPossible() {
        	return (this.s0 * this.s1 * this.s2);
        }
	}

    Intersects getSide(int i, Tripartition trip) {
      //System.out.println("getSide : " + i + " " + trip.cluster1.getBitSet() + " " + trip.cluster2.getBitSet() + " " + trip.cluster3.getBitSet());

		if (trip.cluster1.getBitSet().get(i)) {
			return new Intersects(1,0,0);
		} else if (trip.cluster2.getBitSet().get(i)) {
			return new Intersects(0,1,0);
		} else if (trip.cluster3.getBitSet().get(i)) {
			return  new Intersects(0,0,1);
		} 
		else {
			return  new Intersects(0,0,0);
		}
	}
	
	




    public FrequencyData WeightCalculator(Tripartition[] trips)
    {


        //  System.out.println("Triparitions : " + trips[0] + " " + trips[1] + " " + trips[2]);
        long[] frequencies = {0l , 0l , 0l};
        long mi = 0l; 
        double [] weight = {0l , 0l , 0l};
        int effn = 0;
        Intersects [] allsides = new Intersects[3];
        boolean newTree = true;
        boolean no_match = false; 

        Iterator<STITreeCluster> it = treeAllClusters.iterator();
        Deque<Intersects> [] stack = new Deque [] {new ArrayDeque<Intersects>(), new ArrayDeque<Intersects>(), new ArrayDeque<Intersects>()};

        for (Integer gtb :geneTreesAsInts)
        {
           
            if(newTree){
                // Bringing the full leaf cluster of that gene tree 
                STITreeCluster fullClusterOftheGT = it.next(); 
                
                // Given tripartition, we are checking if the full cluster of the gene tree is intersecting with the clusters of that tripartition
                for(int i = 0 ; i < 3 ; i++){
                    allsides[i] = new Intersects(
                        trips[i].cluster1.getBitSet().intersectionSize(fullClusterOftheGT.getBitSet()), trips[i].cluster2.getBitSet().intersectionSize(fullClusterOftheGT.getBitSet()), trips[i].cluster3.getBitSet().intersectionSize(fullClusterOftheGT.getBitSet())
                    );
                    // System.out.println("allsides : " + allsides[i]);
                }

                newTree = false;
                mi = allsides[0].maxPossible() ; 
                if(mi == 0){
                    no_match = true;
                }
                else {
                    effn++;
                }

            }
            // for new gene tree 
            if(gtb == Integer.MIN_VALUE)
            { 
                if(!no_match){
                    for(int i=0 ; i<3 ; i++)
                    {
                       
                        // double efffreq = (frequencies[i]+0.0)/(2.0*mi); // calculating the effective frequency by dividing the frequency by 2*mi , where mi is the maximum possible frequency [according to Astral]
                        double efffreq = (frequencies[i]+0.0)/(mi);
                        //System.out.println("efffreq : " + efffreq);
                        weight[i] += efffreq;
                       // System.out.println("weight : " + weight[i]);
                    }
                    for(int i=0 ; i<3 ; i++) stack[i].clear();
                    mi = 0;
                    newTree = true;
                    no_match = false;
                    frequencies = new long[]{0l,0l,0l};
                }
            }
            else{
                if(no_match){
                    continue;
                }
                if(gtb>=0) // for leaf node 
                {
                    
                    for(int i = 0 ; i<3 ; i++)
                    {
                        //  System.out.println("gtb for leaf :  "+ gtb + " in tripartition " + trips[i] + " with side "+ getSide(gtb, trips[i]));
                        stack[i].push(getSide(gtb, trips[i]));

                    }
                        
                    
                }
                else if(gtb == -2) // for internal nodes, find information about the bipartitions . 
                {
                    //System.out.println("Gtb for internal node " + gtb);
                    for(int i=0 ; i<3 ; i++)
                    {
                        Intersects side1 = stack[i].pop();
                        Intersects side2 = stack[i].pop();
                        //System.out.println("i = "+ i + "side1 : " + side1 + " side2 : " + side2);
                        Intersects side = new Intersects(side1, side2);
                        stack[i].push(side);
                        //System.out.println(MGDInference_DP.GREEN+ "allcases : " + allcases(side1, side2) + MGDInference_DP.RESET);

                        // System.out.println("i = "+ i + "side1 : " + side1 + " side2 : " + side2);
                        // System.out.println("allcases : " + allcases(side1, side2));


                        frequencies[i]+=allcases(side1,side2) ;
                    }
                }
            }




        }
        
        return new FrequencyData(weight, effn);


    }

    private long allcases(Intersects side1, Intersects side2) {
        // return F(side1.s0, side1.s1, side2.s2) + F(side1.s0, side1.s2, side2.s1) + F(side1.s1, side1.s2, side2.s0)+ F(side1.s0, side2.s1, side2.s2) + F(side1.s1, side2.s0, side2.s2) + F(side1.s2, side2.s0, side2.s1); //! works , accoding to Astral , but gives the same values for all the tripartitions
        //return F(side1.s0, side1.s1, side2.s2) + F(side1.s2, side2.s0, side2.s1) + F(side1.s0 , side2.s0 , side2.s1)+ F(side1.s1 , side2.s0 , side2.s1) + F(side1.s0, side1.s1, side2.s0) + F(side1.s0, side1.s1, side2.s1); //! seems the most logical one , but not working , gives error 
        return F(side1.s0, side1.s1, side2.s2) + F(side1.s2, side2.s0, side2.s1) ; //! Gives the best possible result 
        // return F(side1.s0,side1.s1,side2.s0 )+ F(side1.s0,side1.s1,side2.s1 )+ F(side1.s0,side1.s1,side2.s2 )+ F(side2.s0,side2.s1,side1.s0 )+ F(side2.s0,side2.s1,side1.s1 )+ F(side2.s0,side2.s1,side1.s2 ) ; //! Not working , gives error

       

    }

    
}
