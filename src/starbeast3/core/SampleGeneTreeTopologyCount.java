package starbeast3.core;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;

import beast.app.util.OutFile;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import starbeast3.util.TreeTopologySerialiser;

/**
 * 
 * @author Jordan Douglas
 * A class which: 
 * 1) Simulates a species tree from the prior
 * 2) Simulates G gene trees under the species tree (and the rest of the prior)
 * 3) Records the number of unique topologies among the G gene trees
 * 4) Repeats steps 1-3 'nSamples' times for each G in the range 'nGeneRange'
 * 5) Optionally logs the results using parameter 'numTopologies'
 */
public class SampleGeneTreeTopologyCount extends DirectSimulator{

	
    public Input<Integer> nGeneRangeInput = new Input<>("nGeneRange", "Range of number of gene tree samples.", 
    		Input.Validate.REQUIRED);

    
	final public Input<OutFile> svgOutputInput = new Input<>("svg", "svg output file. if not specified, no SVG output is produced.",
			new OutFile("[[none]]"));
	
	
	//final public Input<State> stateInput = new Input<>("state", "the state containing all parameters, a species tree, and a single gene tree.",
		//	Input.Validate.REQUIRED);
	
	final public Input<Tree> geneTreeInput = new Input<>("geneTree", "the single gene tree.",
			Input.Validate.REQUIRED);
	
	final public Input<GeneTreeForSpeciesTreeDistribution> geneTreePriorInput = new Input<>("geneTreePrior", "the prior distribution behind the gene tree).",
			Input.Validate.REQUIRED);
	
	final public Input<IntegerParameter> numTopologiesInput = new Input<>("numTopologies", "A parameter which will be populated with the number of topologies after each sample.",
			Input.Validate.OPTIONAL);
	
	
	
	GeneTreeForSpeciesTreeDistribution geneTreePrior;
	Tree geneTree;
	int nGeneRange;
	PrintStream svg;
	
	
	 @Override
	 public void initAndValidate() {
		super.initAndValidate();
		 
		
		nGeneRange = nGeneRangeInput.get();
		//state = stateInput.get();
		geneTree = geneTreeInput.get();
		geneTreePrior = geneTreePriorInput.get();
		
		try {
		
		
			svg = null;
			if (svgOutputInput.get() != null && !svgOutputInput.get().getName().equals("[[none]]")) {
				Log.warning("Writing to file " + svgOutputInput.get().getPath());
				svg = new PrintStream(svgOutputInput.get());
				// svg.println(header.replaceAll("file1", src1Input.get().getPath()).replaceAll("file2", src2Input.get().getPath()));
			}
		
		} catch(Exception e) {
			
			e.printStackTrace();
			return;
			
		}
		 
	 }
	

	@Override
    public void doASimulation() {
        clearSampledFlags(distribution);
        
        
        // Set the number of gene trees
        int G = this.nGeneRange;
        HashMap<String, Boolean> topologies = new HashMap<String, Boolean>();
        
        
        // Sample everything
        distribution.sample(state, random);
        
        
        // Cache the sampled gene tree topology
        String newick = TreeTopologySerialiser.serialiseNode(geneTree.getRoot());
        topologies.put(newick, true);
        
        
        // Sample the remaining G-1 gene trees under this species tree / population size
        for (int g = 1; g < G; g ++) {
        	clearSampledFlags(geneTreePrior);
        	geneTreePrior.speciesTreePriorInput.get().sampledFlag = true;
        	geneTreePrior.sample(state, random);
        	
        	//System.out.println(newick);
        	
        	 // Cache the newick
        	 newick = TreeTopologySerialiser.serialiseNode(geneTree.getRoot());
        	 topologies.put(newick, true);
        	
        }
        
        
        // Count the number of gene tree topologies (and possibly log it under numTopologiesInput)
        int numberOfTopologies = topologies.keySet().size();
        if (numTopologiesInput.get() != null) numTopologiesInput.get().setValue(numberOfTopologies);

        
    }
	
	
	
	
	
	
	

}
