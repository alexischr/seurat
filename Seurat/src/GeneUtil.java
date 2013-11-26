/*
 * Copyright (c) 2011 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.refseq.Transcript;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 3/29/11
 * Time: 10:16 AM
 * To change this template use File | Settings | File Templates.
 */
enum GeneRegion {
    UNKNOWN,
    INTRON,
    CCDS,
    UTR,
    UPSTREAM_NEAR
}

public class GeneUtil {
    List<GeneRegion> getRegion(SeekableRODIterator rodIterator, GenomeLoc location) {

        RODRecordList annotationList = (rodIterator == null ? null : rodIterator.seekForward(location));
        List<GeneRegion> list = new ArrayList<GeneRegion>();

        for (GATKFeature c_feature : annotationList) {
            Transcript c_transcript = (Transcript) c_feature.getUnderlyingObject();

            if (c_transcript == null)
                list.add(GeneRegion.UNKNOWN);
            else {
                if (c_transcript.overlapsExonP(location)) {
                    if (c_transcript.overlapsCodingP(location))
                        list.add(GeneRegion.CCDS);
                    else
                        list.add(GeneRegion.UTR);
                } else {
                    list.add(GeneRegion.INTRON);
                }

            }
        }
        return list;
    }
}

