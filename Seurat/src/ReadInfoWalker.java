/*
 * Copyright (c) 2011 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.io.PrintStream;


/**
 * Paired Sample Mutation Caller
 */

@WalkerName("ReadInfoWalker")
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})
// This walker requires -I input.bam, it also requires -R reference.fasta
public class ReadInfoWalker extends LocusWalker<Integer, Long> {
    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return super.filter(tracker, ref, context);    //To change body of overridden methods use File | Settings | File Templates.
    }


    @Output
    private PrintStream out;

    // control the various parameters to be used
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALTY_SCORE = 10;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 10;

    @Override
    public void initialize() {
        //Separate read groups by BAM
    }


    /**
     * The map function runs once per single-base locus, and accepts a 'context', a
     * data structure consisting of the reads which overlap the locus, the sites over
     * which they fall, and the base from the reference that overlaps.
     *
     * @param tracker The accessor for reference metadata.
     * @param ref     The reference base that lines up with this locus.
     * @param context Information about reads aligning to this locus.
     * @return In this case, returns a count of how many loci were seen at thisb site (1).
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        char refBase = (char) ref.getBase();
        double ref_count = 0;
        double nonref_count = 0;
        double ar;


        for (PileupElement p : context.getBasePileup()) {
            SAMRecord read = p.getRead();

            //filter
            if (read.getMappingQuality() < MIN_MAPPING_QUALTY_SCORE || p.getQual() < MIN_BASE_QUALTY_SCORE)
                continue;

            char base = (char) p.getBase();

            if (base != refBase)
                nonref_count++;
            else
                ref_count++;

        }

        out.printf("%s\t%d\t%.0f\t%.0f\t%.3f\n", context.getLocation().getContig(), context.getLocation().getStart(), nonref_count, ref_count, nonref_count / (nonref_count + ref_count));

        return 0;


    }


    /**
     * Provides an initial value for the reduce function.  Hello walker counts loci,
     * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
     *
     * @return 0.
     */
    @Override
    public Long reduceInit() {
        return 0L;
    }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return The total count of loci processed so far.
     */
    @Override
    public Long reduce(Integer value, Long sum) {
        return sum + value;
    }

    /**
     * Retrieves the final result of the traversal.
     *
     * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
     *               by the reduce function.
     */
    @Override
    public void onTraversalDone(Long result) {
        //out.println("Number of mutations found: " + result);
    }
}
