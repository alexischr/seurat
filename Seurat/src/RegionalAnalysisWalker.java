/*
 * Copyright (c) 2012 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/21/11
 * Time: 11:21 AM
 * To change this template use File | Settings | File Templates.
 */

public abstract class RegionalAnalysisWalker {

    public abstract boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context);

    public abstract void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts);

    public abstract boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);

    public abstract boolean isTranscriptionalAnalysis();

    public abstract List<Event> reduce();

}
