package starbeast3.operators;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import beast.core.*;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.math.matrixalgebra.CholeskyDecomposition;
import beast.math.matrixalgebra.IllegalDimension;
import beast.util.Randomizer;

@Description("Adaptable variance multivariate normal operator on cube description of a tree")
public class AVMNCubeOperator extends TreeOperator {
	final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "start scaling factor, larger values give bolder moves (this is tuned during the run)", 1.0); 
	final public Input<Double> coefficientInput = new Input<>("coefficient", "determines diagonal correlation for variance matrix", 1.0);
	final public Input<Double> betaInput = new Input<>("beta", "fraction of proposal determined by non-covariance matrix");
	final public Input<Integer> initialInput = new Input<>("initial", "Number of proposals before covariance matrix is considered in proposal. "
			+ "Must be larger than burnin, if specified. "
			+ "If not specified (or < 0), the operator uses 200 * paramater dimension", -1); 
	final public Input<Integer> burninInput = new Input<>("burnin", "Number of proposals that are ignored before covariance matrix is being updated. "
			+ "If initial is not specified, uses half the default initial value (which equals 100 * paramater dimension)", 0); 
	final public Input<Integer> everyInput = new Input<>("every", "update interval for covarionce matrix, default 1 (that is, every step)", 1); 
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    
    TreeInterface tree;
	
    
    private double scaleFactor;
    private double beta;
    private int iterations, updates, initial, burnin, every;
    private int dim;

    private double[] oldMeans, newMeans;

    private double[][] matrix;
    private double[][] empirical;
    private double[][] cholesky;

    // temporary storage, allocated once.
    private double[] epsilon;
    private double[][] proposal;
    
	double [] heights;
	int [] order;
    
	public AVMNCubeOperator() {}
	
	public AVMNCubeOperator(Tree tree, double weight) {
		initByName("tree", tree, "weight", weight);
	}
	

	@Override
	public void initAndValidate() {
	}
	
	private void setUp() {
		tree = treeInput.get();
		leafNodeCount = tree.getLeafNodeCount();
		heights = new double[leafNodeCount-1];
		order = new int[2*leafNodeCount - 1];
		Node [] internalNodes = new Node[leafNodeCount - 1];
		getCube(tree.getRoot(), order, heights, new int[]{0}, internalNodes);

	
	
        this.scaleFactor = scaleFactorInput.get();
        this.beta = betaInput.get();
        this.iterations = 0;
        this.updates = 0;
        
        dim = heights.length;
        
        matrix = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
        	matrix[i][i] = Math.pow(coefficientInput.get(), 2) / (dim);            	
        }        
        
        
        // constantFactor = Math.pow(2.38, 2) / ((double) dim); // not necessary because scaleFactor is auto-tuned
        this.initial = initialInput.get();
        
        if (this.initial < 0) {
        	// options set according to recommendations in AVMVN paper
        	this.initial = 200 * dim;
        	this.burnin = this.initial / 2;
        } else {
            this.burnin = burninInput.get();
        }

        
        this.every = everyInput.get();
        this.empirical = new double[dim][dim];
        this.oldMeans = new double[dim];
        this.newMeans = new double[dim];

        this.epsilon = new double[dim];
        this.proposal = new double[dim][dim];

    	
    	
        if (burnin > initial || burnin < 0) {
            throw new IllegalArgumentException("Burn-in must be smaller than the initial period.");
        }


        if (every <= 0) {
            throw new IllegalArgumentException("Covariance matrix needs to be updated at least every single iteration.");
        }

        if (scaleFactor <= 0.0) {
            throw new IllegalArgumentException("ScaleFactor must be greater than zero.");
        }
        
        try {
            cholesky = (new CholeskyDecomposition(matrix)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Unable to decompose matrix in AdaptableVarianceMultivariateNormalOperator");
        }		
	}

	int leafNodeCount = 0;
	
	
    //act as if population mean is known
    private double calculateCovariance(int number, double currentMatrixEntry, double[] values, int firstIndex, int secondIndex) {

        // number will always be > 1 here
        /*double result = currentMatrixEntry * (number - 1);
        result += (values[firstIndex] * values[secondIndex]);
        result += ((number - 1) * oldMeans[firstIndex] * oldMeans[secondIndex] - number * newMeans[firstIndex] * newMeans[secondIndex]);
        result /= ((double) number);*/

        double result = currentMatrixEntry * (number - 2);
        result += (values[firstIndex] * values[secondIndex]);
        result += ((number - 1) * oldMeans[firstIndex] * oldMeans[secondIndex] - number * newMeans[firstIndex] * newMeans[secondIndex]);
        result /= (number - 1);
        return result;

    }


	@Override
	public double proposal() {
		if (order == null) {
			setUp();
		}

		tree2Cube(order, tree, heights);

        double[] x = new double[dim];
        System.arraycopy(heights, 0, x, 0, dim);

        //transform to the appropriate scale
        double[] transformedX = new double[dim];

        for (int i = 0; i < dim; i++) {
            transformedX[i] = transform(x[i]);
        }

        //store MH-ratio in logq
        double logJacobian = 0.0;

        //change this: make a rule for when iterations == burnin
        if (iterations > 1 && iterations > burnin) {

            if (iterations > (burnin+1)) {

                if (iterations % every == 0) {

                    updates++;

                    //first recalculate the means using recursion
                    for (int i = 0; i < dim; i++) {
                        newMeans[i] = ((oldMeans[i] * (updates - 1)) + transformedX[i]) / updates;
                    }

                    if (updates > 1) {
                        //here we can simply use the double[][] matrix
                        for (int i = 0; i < dim; i++) {
                            for (int j = i; j < dim; j++) {
                                empirical[i][j] = calculateCovariance(updates, empirical[i][j], transformedX, i, j);
                                empirical[j][i] = empirical[i][j];
                            }
                        }
                    }

                }

                /*if (iterations == 17) {
                    System.exit(0);
                }*/

            } else if (iterations == (burnin+1)) {

                //updates++;

                //i.e. iterations == burnin+1, i.e. first sample for C_t
                //this will not be reached when burnin is set to 0
                for (int i = 0; i < dim; i++) {
                    //oldMeans[i] = transformedX[i];
                    //newMeans[i] = transformedX[i];
                    oldMeans[i] = 0.0;
                    newMeans[i] = 0.0;
                }

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        empirical[i][j] = 0.0;
                    }
                }

            }

            // End of adaptable covariance -- move into separate class

        } else if (iterations == 1) {

            //iterations == 1
            for (int i = 0; i < dim; i++) {
                //oldMeans[i] = transformedX[i];
                //newMeans[i] = transformedX[i];
                oldMeans[i] = 0.0;
                newMeans[i] = 0.0;
            }

            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    empirical[i][j] = 0.0;
                    proposal[i][j] = matrix[i][j];
                }
            }

        }

        for (int i = 0; i < dim; i++) {
            epsilon[i] = scaleFactor * Randomizer.nextGaussian();
        }

        if (iterations > initial) {

            if (iterations % every == 0) {
                // For speed, it may not be necessary to update decomposition each and every iteration

                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) { // symmetric matrix
                        proposal[j][i] = (1 - beta) * // constantFactor *  /* auto-tuning using scaleFactor */
                                empirical[j][i] + beta * matrix[j][i];
                        proposal[i][j] = proposal[j][i] ;
                    }
                }

                // not necessary for first test phase, but will need to be performed when covariance matrix is being updated
                try {
                    cholesky = (new CholeskyDecomposition(proposal)).getL();
                } catch (IllegalDimension illegalDimension) {
                    throw new RuntimeException("Unable to decompose matrix in AdaptableVarianceMultivariateNormalOperator");
                }
            }

        }

        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++) {
                transformedX[i] += cholesky[j][i] * epsilon[j];
            }
        }

        for (int i = 0; i < dim; i++) {
             x[i] = inverse(transformedX[i]);
             logJacobian += getLogJacobian(x[i]) - getLogJacobian(heights[i]);
        }

        if (iterations % every == 0) {
            //copy new means to old means for next update iteration
            //System.arraycopy(newMeans, 0, oldMeans, 0, dim);
            double[] tmp = oldMeans;
            oldMeans = newMeans;
            newMeans = tmp; // faster to swap pointers
        }

		cube2Tree(order, x, tree);

        return logJacobian;

	}
 	
    public double transform(double value) {
        return Math.log(value);
    }

    public double inverse(double value) {
        return Math.exp(value);
    }

    public double gradientInverse(double value) { return Math.exp(value); }

    public double updateGradientLogDensity(double gradient, double value) {
        // gradient == gradient of inverse()
        // value == gradient of inverse() (value is untransformed)
        // 1.0 == gradient of log Jacobian of inverse()
        return gradient * value + 1.0;
    }

    protected double getGradientLogJacobianInverse(double value) {
        return 1.0;
    }

    public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
        // value == inverse()
        // diagonalHessian == hessian of inverse()
        // gradient == gradient of inverse()
        return value * (gradient + value * diagonalHessian);
    }

    public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transfomationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
        return offdiagonalHessian * valueI * valueJ + gradientJ * transfomationHessian;
    }

    public double gradient(double value) {
        return value;
    }

    public double getLogJacobian(double value) { return -Math.log(value); }

	
	
	
	/** for a given order and tree, determine heights of a cube **/
	private void tree2Cube(int[] order, TreeInterface tree, double[] heights) {
		boolean [] done = new boolean[tree.getNodeCount()];
		for (int j = 0; j < heights.length; j++) {
			Node node = tree.getNode(order[j]);
			node = node.getParent();
			while (done[node.getNr()]) {
				node = node.getParent();
			}
			heights[j] = node.getHeight();
			done[node.getNr()] = true;
		}
	}

	/** recursively determines order of tree compatible with (random) planar drawing of the tree 
	 * @param node current node in tree
	 * @param order is populated with ordering of leaf nodes
	 * @param heights of the cube
	 * @param index keeps track of next order index to add
	 * @param internalNodes nodes for which heights are recorded
	 */
	private void getCube(Node node, int [] order, double [] heights, int [] index, Node [] internalNodes) {
		if (node.isLeaf()) {
			order[index[0]] = node.getNr();
		} else {
			Node first = null, second = null;
			if (Randomizer.nextBoolean()) {
				first = node.getChild(0);
				second = node.getChild(1);
			} else {
				first = node.getChild(1);
				second = node.getChild(0);
			}
			getCube(first, order, heights, index, internalNodes);
			heights[index[0]] = node.getHeight();
			internalNodes[index[0]] = node;
			order[index[0]+leafNodeCount] = node.getNr();
			index[0]++;
			getCube(second, order, heights, index, internalNodes);
		}		
	}

	
	final static boolean debug = false;
	public static int attempts = 0;
	public static int attemptsMax = 0;

	
	private Node cube2Tree(int [] order, double [] values, TreeInterface tree) {
		int N = tree.getLeafNodeCount();
		
		Node [] nodes = tree.getNodesAsArray();
		
		double [] heights = new double[tree.getNodeCount() + 1];
		int [] parent_ = new int[tree.getNodeCount()];
		int [] left = new int[tree.getNodeCount() + 1];
		int [] right = new int[tree.getNodeCount() + 1];
		Arrays.fill(parent_, -1);
		Arrays.fill(left, -1);
		Arrays.fill(right, -1);
		
		int root_ = tree.getNodeCount();
		heights[root_] = Double.MAX_VALUE;
		left[root_] = order[0];
		parent_[order[0]] = root_;
		int current_;
		int next_ = N;
		for (int i = 0; i < N - 1; i++) {
			current_ = order[i];
			double target = values[i];
			while (target > heights[parent_[current_]]) {
				current_ = parent_[current_];
			}
			int parent = parent_[current_];
			int internal = order[next_];
			if (left[parent] == current_) {
				left[parent] = internal;
			} else {
				// note, then right[parent] == current_
				right[parent] = internal;
			}
			parent_[internal] = parent;
			
			left[internal] = current_;
			parent_[current_] = internal;
			right[internal] = order[i+1];
			parent_[order[i+1]] = internal;
			heights[internal] = target;
			next_++;
		}
		
		
		List<Node> dirtyNodes = new ArrayList<>();
		for (int i = N; i < nodes.length; i++) {
			Node node = nodes[i];
			if (Math.abs(node.getHeight() - heights[i]) > 1e-10) {
				node.setHeight(heights[i]);
			}
			if (node.getLeft().getNr() != left[i]) {
				node.setLeft(nodes[left[i]]);
				nodes[left[i]].setParent(node);
				dirtyNodes.add(node);
			} 
			if (node.getRight().getNr() != right[i]) {
				node.setRight(nodes[right[i]]);
				nodes[right[i]].setParent(node);
				dirtyNodes.add(node);
			}
			if (parent_[i] == nodes.length) {
				node.setParent(null);
			}
		}
		if (left[nodes.length] != nodes.length-1) {
			int treeroot = left[nodes.length];
			((Tree) tree).setRoot(nodes[treeroot]);
		}
		
		
		int max = 0;
		for (Node node : dirtyNodes) {
			int i = 0;
			while (node.getParent() != null) {
				node = node.getParent();
				i++;
			}
			max = Math.max(i,  max);
		}
		attempts++;
		attemptsMax += max;
		
		
		if (left[nodes.length] != nodes.length-1) {
			String newick2 = tree.toString();
			for (int i = N; i < nodes.length; i++) {
				nodes[i].removeAllChildren(false);
			}
			Node root = new Node();
			root.setHeight(Double.MAX_VALUE);
			root.addChild(nodes[order[0]]);
	
			Node current;
			int next = N;
			for (int i = 0; i < N - 1; i++) {
				current = nodes[order[i]];
				double target = values[i];
				while (target > current.getParent().getHeight()) {
					current = current.getParent();
				}
				Node parent = current.getParent();
				parent.removeChild(current);
				Node internal = nodes[order[next]];
				parent.addChild(internal);
				internal.addChild(current);
				internal.addChild(nodes[order[i+1]]);
				internal.setHeight(target);
				next++;
			}
			
			root = root.getChild(0);
			root.setParent(null);
			((Tree)tree).setRoot(root);
			//return root;

			String newick = tree.toString();
			if (!newick.equals(newick2)) {
				System.err.println("Something went wrong: \n" + newick + "\n" + newick2);
				cube2Tree(order, values, tree);
			}
		}
		return tree.getRoot();
	}
    
}
