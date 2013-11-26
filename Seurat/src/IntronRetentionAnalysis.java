/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
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
 * Date: 1/11/12
 * Time: 4:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class IntronRetentionAnalysis extends RegionalAnalysisWalker {

    SeuratArgumentCollection arguments;
    double min_event_p;
    List<Event> temp_events = new ArrayList<Event>();

    PileupEvidence rna_normal_evidence;
    PileupEvidence rna_tumor_evidence;


    @Override
    public boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context) {
        return true;
    }

    public void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts) {

    }

    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean isTranscriptionalAnalysis() {
        return true;
    }

    @Override
    public List<Event> reduce() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
