/*
 * Copyright (c) 2012 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 1/4/12
 * Time: 11:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExonDeletionAnalysis extends RegionalAnalysisWalker {

    double normal_count = 0;
    double tumor_count = 0;

    PileupEvidence normal_evidence;
    PileupEvidence tumor_evidence;

    GeneContext geneContext;

    @Override
    public boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context) {
        normal_evidence = Evidence.get("rna_normal");
        tumor_evidence = Evidence.get("rna_tumor");
        geneContext = context;

        return (normal_evidence != null && tumor_evidence != null && context.context == GeneContext.GeneContextClass.Exon);
    }

    @Override
    public void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts) {
        normal_count += normal_evidence.Pileup.size();
        tumor_count += tumor_evidence.Pileup.size();
    }

    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true; //
    }

    @Override
    public boolean isTranscriptionalAnalysis() {
        return true;
    }

    @Override
    public List<Event> reduce() {
        List<Event> le = new ArrayList<Event>();

        if (normal_count + tumor_count != 0) {
            double freq = normal_count / (normal_count + tumor_count);
            if (freq > 0.9) {
                Event e = new Event("exon_silencing", freq, geneContext);
                le.add(e);
            }
        }
        return le;
    }
}
