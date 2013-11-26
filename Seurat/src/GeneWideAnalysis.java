/*
 * Copyright (c) 2011 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.List;
import java.util.Map;


public abstract class GeneWideAnalysis {

    protected List<Event> eventList;
    protected String gene;

    public void beginGene(String geneName, List<Event> eventList) {
        this.gene = geneName;
        this.eventList = eventList;
    }

    public List<Event> endGene() {
        return eventList;
    }

    public abstract void map(Map<String, PileupEvidence> evidenceList, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);

}
