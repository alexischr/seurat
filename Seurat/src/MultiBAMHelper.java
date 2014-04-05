/*
 * Copyright (c) 2014 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.HashMap;
import java.util.Map;

public class MultiBAMHelper {

    //private
    private HashMap<String, String> sample_set;
    public HashMap<String, PileupEvidence> Evidence;

    public int getMinPerSampleCoverage() {
        return MinPerSampleCoverage;
    }

    public void setMinPerSampleCoverage(int minPerSampleCoverage)
    {
        MinPerSampleCoverage = minPerSampleCoverage;
        minimum_total_pileup = total_inputs * MinPerSampleCoverage;
    }

    private int MinPerSampleCoverage = 6;
    private int max_mismatches = -1;

    private int rna_tags = 0;
    private int total_inputs = 0;
    private int minimum_total_pileup = 0;

    public void setMaxMismatches(int maxMismatches) {
        max_mismatches = maxMismatches;
    }

    public int getMaxMismatches() {
        return max_mismatches;
    }


    public MultiBAMHelper(GenomeAnalysisEngine engine) {
        sample_set = new HashMap<String, String>();
        Evidence = new HashMap<String, PileupEvidence>();

        for (SAMReaderID rid : engine.getReadsDataSource().getReaderIDs()) {
            total_inputs++;
            Tags tags = rid.getTags();
            if (tags.getPositionalTags().isEmpty())
                throw new UserException.BadInput("This module requires that input BAMs are tagged (ie. -I:dna_normal normal.bam -I:dna_tumor tumor.bam): " +
                        engine.getSourceFileForReaderID(rid));

            for (String tag : tags.getPositionalTags()) {
                for (SAMReadGroupRecord rg : engine.getSAMFileHeader(rid).getReadGroups()) {

                    String existing_tag = sample_set.get(rg.getSample());

                    if (existing_tag != null && !existing_tag.equals(tag)) {
                        throw new UserException.BadInput(String.format("Sample '%s' is being assigned to tag '%s', but has already been assigned to tag '%s'. This is often a result of tumor/normal BAMs using the same sample IDs in their @RG tags. \n", rg.getSample(), tag, existing_tag));
                    }

                    sample_set.put(rg.getSample(), tag);
                    Evidence.put(tag, new PileupEvidence());

                    if (tag.startsWith("rna")) {
                        rna_tags++;
                    }
                }

            }

        }

    }

    public boolean SetBasePileup(ReadBackedExtendedEventPileup pileup, ReferenceContext ref_context, boolean contains_insertions) {

        if (rna_tags == 0 && pileup.getBases().length < minimum_total_pileup) //if there's not enough pileup to satisfy the minimum coverage requirement for all of the tags, quit early
            //note: this assumes non-overlapping tags!
            return false;

        byte ref_base;

        if (contains_insertions) { //pileup is actually of the *next* base
            ref_base = ref_context.getBases()[101];
        } else
            ref_base = ref_context.getBase();

        for (PileupEvidence bam_pileup : Evidence.values()) {
            bam_pileup.Reset(ref_base);
        }

        for (ExtendedEventPileupElement p : pileup.toExtendedIterable()) {

            GATKSAMRecord read = p.getRead();

            //filter
            if (isMismatchFiltered(ref_context, read)) continue;

            if (read.getDuplicateReadFlag() || read.getReadUnmappedFlag())
                continue;

            String sample = read.getReadGroup().getSample();

            Evidence.get(sample_set.get(sample)).Add(p);
        }

        for (Map.Entry<String, PileupEvidence> bam_pileup : Evidence.entrySet()) {
            if (!bam_pileup.getKey().startsWith("rna") && bam_pileup.getValue().Pileup.size() < MinPerSampleCoverage) {
                return false;
            }
        }

        return true;
    }

    public boolean SetBasePileup(ReadBackedPileup pileup, ReferenceContext ref_context) {

        if (rna_tags == 0 && pileup.getBases().length < minimum_total_pileup) //if there's not enough pileup to satisfy the minimum coverage requirement for all of the tags, quit early
            //note: this assumes non-overlapping tags!
            return false;

        for (PileupEvidence bam_pileup : Evidence.values()) {
            bam_pileup.Reset(ref_context.getBase());
        }

        for (PileupElement p : pileup) {

            GATKSAMRecord read = p.getRead();

            if (isMismatchFiltered(ref_context, read)) continue;

            //filter
            if (read.getDuplicateReadFlag() || read.getReadUnmappedFlag())
                continue;

            String sample = read.getReadGroup().getSample();

            Evidence.get(sample_set.get(sample)).Add(p);
        }

        for (Map.Entry<String, PileupEvidence> bam_pileup : Evidence.entrySet()) {
            if (!bam_pileup.getKey().startsWith("rna") && bam_pileup.getValue().Pileup.size() < MinPerSampleCoverage) {
                return false;
            }

        }

        return true;
    }

    private boolean isMismatchFiltered(ReferenceContext ref_context, GATKSAMRecord read) {
        if (max_mismatches > -1) {

            AlignmentUtils.MismatchCount mismatches = (AlignmentUtils.MismatchCount) read.getTemporaryAttribute("NMM");
            if (mismatches == null) {
                int window_offset = read.getAlignmentStart() - ref_context.getWindow().getStart();

                if (window_offset >= 0) {
                    mismatches = AlignmentUtils.getMismatchCount(read, ref_context.getBases(), window_offset);
                    read.setTemporaryAttribute("NMM", mismatches);
                }
            }
            if (mismatches != null && mismatches.numMismatches > max_mismatches)
                return true;
        }
        return false;
    }

}
