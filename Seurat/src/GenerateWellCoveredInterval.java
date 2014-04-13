/*
 * Copyright (c) 2014 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.WalkerName;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * Walker for creating interval lists with a set minimum coverage (for use with the -L GATK parameter)
 */

@WalkerName("Generate ")
@Requires({DataSource.READS, DataSource.REFERENCE, DataSource.REFERENCE_BASES})
// This walker requires -I input.bam, it also requires -R reference.fasta
public class GenerateWellCoveredInterval extends LocusWalker<Integer, Long> {
    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return multibam_helper.SetBasePileup(context.getBasePileup().getBaseAndMappingFilteredPileup(MIN_BASE_QUALITY_SCORE, MIN_MAPPING_QUALITY_SCORE), ref);
    }


    @Output
    private PrintStream out;

    // control the various parameters to be used
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", required = false)
    public int MIN_BASE_QUALITY_SCORE = 10;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling", required = false)
    public int MIN_MAPPING_QUALITY_SCORE = 10;

    @Argument(fullName = "min_coverage", shortName = "cov", doc = "Minimum coverage", required = false)
    public int MIN_COVERAGE = 5;


    public List<PileupEvidence> ai_list;
    private MultiBAMHelper multibam_helper = null;
    private PileupEvidence DNA_Normal = new PileupEvidence();
    private PileupEvidence DNA_Tumor = new PileupEvidence();

    @Override
    public void initialize() {

        ai_list = new ArrayList<PileupEvidence>(Arrays.asList(DNA_Normal, DNA_Tumor));
        multibam_helper = new MultiBAMHelper(getToolkit(), false);

        multibam_helper.setMinPerSampleCoverage(MIN_COVERAGE);
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
        out.printf("%s:%d\n", context.getLocation().getContig(), context.getLocation().getStart());

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
