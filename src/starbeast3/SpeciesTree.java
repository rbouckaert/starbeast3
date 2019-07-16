package starbeast3;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

/**
 * adapted from starbeast2
* @author Huw Ogilvie
 */

public class SpeciesTree extends Tree {
    Map<String, Integer> tipNumberMap;
    Map<Integer, Set<String>> numberTipMap;

    public void initAndValidate() {
        super.initAndValidate();
        tipNumberMap = new LinkedHashMap<>();
        numberTipMap = new LinkedHashMap<>();
        makeMaps(tipNumberMap, numberTipMap);
    }

    public Map<String, Integer> getTipNumberMap() {
        return tipNumberMap;
    }

    public Map<Integer, Set<String>> getNumberTipMap() {
        return numberTipMap;
    }

	public void adjustTreeNodeHeights() {
		adjustTreeNodeHeights(root);
	}
	
    // generate map of species tree tip node names to node numbers
    private void makeMaps(Map<String, Integer> tipNumberMapx, Map<Integer, Set<String>> numberTipMap) {
        final Map<String, Integer> speciesNumberMap = new LinkedHashMap<>();
        Node speciesTreeRoot = getRoot();
        for (Node leafNode: speciesTreeRoot.getAllLeafNodes()) {
            final String speciesName = leafNode.getID();
            final int speciesNumber = leafNode.getNr();
    
            speciesNumberMap.put(speciesName, speciesNumber);
        }

        // generate map of gene tree tip node names to species tree tip node numbers
        final TaxonSet taxonSuperSet = getTaxonset();
        final Set<Taxon> speciesSet = new LinkedHashSet<>(taxonSuperSet.taxonsetInput.get());

        for (Taxon species: speciesSet) {
            final String speciesName = species.getID();
            int speciesNumber = 0;
            if (speciesNumberMap.containsKey(speciesName)) { // skipped for BEAUTi
                speciesNumber = speciesNumberMap.get(speciesName);
            }
            final TaxonSet speciesTaxonSet = (TaxonSet) species;
            final Set<Taxon> tipSet = new LinkedHashSet<>(speciesTaxonSet.taxonsetInput.get());

            for (Taxon tip: tipSet) {
                final String tipName = tip.getID();
                tipNumberMap.put(tipName, speciesNumber);
                if (!numberTipMap.containsKey(speciesNumber)) {
                    numberTipMap.put(speciesNumber, new LinkedHashSet<>());                	
                }
                numberTipMap.get(speciesNumber).add(tipName);
            }
        }
    }
	
}
