/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.commandline.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/27/11
 * Time: 11:38 AM
 * To change this template use File | Settings | File Templates.
 */
public class SeuratArgumentCollection {
    @Argument(fullName = "prior_alpha", shortName = "alpha", doc = "alpha parameter for the homozygosity beta distribution (default = 1)", required = false)
    public int beta_alpha = 1;

    @Argument(fullName = "prior_beta", shortName = "beta", doc = "beta parameter for the homozygosity beta distribution (default = 701)", required = false)
    public int beta_beta = 700;

    @Argument(fullName = "refnormal_only", shortName = "ref", required = false, doc = "Whether or not only reference-matching homozygous positions are allowed on the normal, for SNV discovery. Reduces false positives due to faulty alignments (default = true)")
    public boolean refnormal_only = true;

    @Argument(fullName = "both_strands", shortName = "both_strands", required = false, doc = "Whether or not variant evidence needs to appear on both strands on the tumor in order to be considered.")
    public boolean both_strands = false;

    @Argument(fullName = "coding_only", shortName = "coding", required = false, doc = "Process only coding regions.")
    public boolean coding_only = false;

    @Argument(fullName = "rna_snv", shortName = "rna_snv", required = false, doc = "Merge RNA evidence to DNA evidence when performing somatic mutation analysis")
    public boolean rna_snv = false;

    @Argument(fullName = "pileup_info", shortName = "pileup", required = false, doc = "Output pileup on each event record.")
    public boolean pileup_info = true;

    @Argument(fullName = "min_event_quality", shortName = "Q", doc = "Filter events with lower that this phred-scale quality.", required = false)
    public double quality = 10;

    @Argument(fullName = "expected_insert_size", shortName = "insert_size", required = false, doc = "Expected insert size for structural variation detection.")
    public int expected_insert_size = 1500;

    @Argument(fullName = "maximum_mismatches", shortName = "mm", required = false, doc = "Maximum number of mismatches against the reference that are allowed in a read. Set to a negative number to disable.")
    public int maximum_mismatches = 3;

    @Argument(fullName = "loh", required = false, doc = "Enable detection of loss of heterozygosity (LOH) loci in the tumor.")
    public boolean enable_loh = false;

    @Argument(fullName = "structvar", required = false, doc = "Enable detection of loss of somatic structural variation (large deletions/translocations).")
    public boolean enable_structvar = false;

    @Argument(fullName = "metrics", required = false, doc = "Provide additional metrics for calls (eg. median base quality, median mapping quality, median cycle).")
    public boolean enable_metrics = false;

    @Argument(fullName = "debug", required = false, doc = "Enable debugging mode/debugging information. Do not use if you're not familiar with the effects.")
    public boolean enable_debug = false;


    @Override
    public String toString() {
        return String.format("snv_alpha=%d;snv_beta=%d;expected_insert_size=%d;min_event_quality=%.1f;maximum_mismatches=%d", beta_alpha,
                beta_beta, expected_insert_size, quality, maximum_mismatches);

    }
}
