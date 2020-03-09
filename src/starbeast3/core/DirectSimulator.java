package starbeast3.core;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.Runnable;
import beast.util.Randomizer;


import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

/**
 * Runnable for generating a fixed number of samples from a prior distribution using
 * direct simulation.
 *
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 16/06/17.
 * Modified for StarBeast3 by Jordan Douglas 
 */
@Description("Runnable for generating a fixed number of samples from a prior distribution" +
             "using direct simulation.")
public class DirectSimulator extends Runnable {

    public Input<Distribution> distributionInput = new Input<>("distribution",
            "Distribution to sample from via direct simulation.",
            Input.Validate.REQUIRED);
    
    // Init nodes
    final public Input<List<StateNodeInitialiser>> initialisersInput =
            new Input<>("init", "one or more state node initilisers used for determining " +
                    "the start state of the chain",
                    new ArrayList<>());

    public Input<List<Logger>> loggersInput = new Input<>("logger",
            "Log that simulated parameters and trees can be written to.",
            new ArrayList<>());

    public Input<Integer> nSamplesInput = new Input<>("nSamples",
            "Number of independent samples to generate.",
            Input.Validate.REQUIRED);
    


    State state;
    Distribution distribution;
    List<Logger> loggers;
    int nSamples;

    Random random;

    @Override
    public void initAndValidate() {
        distribution = distributionInput.get();
        loggers = loggersInput.get();
        nSamples = nSamplesInput.get();
        
        
        
        // StateNode initialisation, only required when the state is not read from file
        if (restoreFromFile) {
            final HashSet<StateNode> initialisedStateNodes = new HashSet<>();
            for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
                // make sure that the initialiser does not re-initialises a StateNode
                final List<StateNode> list = new ArrayList<>(1);
                initialiser.getInitialisedStateNodes(list);
                for (final StateNode stateNode : list) {
                    if (initialisedStateNodes.contains(stateNode)) {
                        throw new RuntimeException("Trying to initialise stateNode (id=" + stateNode.getID() + ") more than once. " +
                                "Remove an initialiser from MCMC to fix this.");
                    }
                }
                initialisedStateNodes.addAll(list);
                // do the initialisation
                //initialiser.initStateNodes();
            }
        }


        // Create new Random instance initialized with seed from Randomizer
        // (Necessary because sample() currently needs a Random instance and we
        // want to use the same seed specified on the command line.)
        random = new Random(Randomizer.getSeed());
    }

    public void clearSampledFlags(BEASTInterface obj) {
        if (obj instanceof Distribution)
            ((Distribution) obj).sampledFlag = false;

        for (String inputName : obj.getInputs().keySet()) {
            Input input = obj.getInput(inputName);

            if (input.get() == null)
                continue;

            if (input.get() instanceof List) {
                for (Object el : ((List)input.get())) {
                    if (el instanceof BEASTInterface)
                        clearSampledFlags((BEASTInterface)el);
                }
            } else if (input.get() instanceof BEASTInterface) {
                clearSampledFlags((BEASTInterface)(input.get()));
            }
        }
    }

    @Override
    public void run() throws Exception {

        
        
        //do {
        for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
            initialiser.initStateNodes();
        }
            //oldLogLikelihood = state.robustlyCalcPosterior(posterior);
            //initialisationAttempts += 1;
       // } while (Double.isInfinite(oldLogLikelihood) && initialisationAttempts < numInitializationAttempts.get());


        // Initialize loggers
        for (Logger logger : loggers)
            logger.init();

        // Perform simulations
        for (int i=0; i<nSamples; i++) {
        	
        	this.doASimulation();

            for (Logger logger : loggers) {
                logger.log(i);
            }
        }

        // Finalize loggers
        for (Logger logger: loggers)
            logger.close();

        System.out.println("Direct simulation of " + nSamples + " samples completed.");
    }
    
    
    public void doASimulation() {
        clearSampledFlags(distribution);
        distribution.sample(state, random);
    }
    
    
    
}





