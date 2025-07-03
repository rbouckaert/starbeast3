package starbeast3.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.evolution.alignment.Alignment;
import beastfx.app.inputeditor.BeautiDoc;
import morphmodels.app.beauti.BeautiMorphModelAlignmentProvider;
import starbeast3.tree.SpeciesTree;
import starbeast3.tree.StarBeast3TaxonSet;

import java.util.List;

@Description("Class for creating new partitions for morphological data to be edited by AlignmentListInputEditor")
public class StarBeast3MorphModelAlignmentProvider extends BeautiMorphModelAlignmentProvider {
    @Override
    public void processAlignment(Alignment alignment, List<BEASTInterface> filteredAlignments, boolean ascertained, BeautiDoc doc) throws Exception {
        StarBeast3TaxonSet ts = (StarBeast3TaxonSet) doc.pluginmap.get("taxonsuperset");
        ts.alignmentInput.set(alignment);
        ts.initAndValidate();

        SpeciesTree st = (SpeciesTree) doc.pluginmap.get("Tree.t:Species");
        st.m_taxonset.set(ts);
        st.makeCaterpillar(0, 1, false);

        super.processAlignment(alignment, filteredAlignments, ascertained, doc);
    }
}