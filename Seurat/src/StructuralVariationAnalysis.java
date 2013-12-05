/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 1/11/12
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class StructuralVariationAnalysis extends RegionalAnalysisWalker {

    SeuratArgumentCollection arguments;
    double min_event_p;

    PileupEvidence normal_evidence;
    PileupEvidence tumor_evidence;

    //bunch of values here
    static final double prior_sv = 0.001;
    static final int cluster_size = 100000;
    static final int normal_window_size = 1000;
    static final int tumor_window_size = 1000;
    static final int sv_alpha = 1;
    static final int sv_beta = 1001;

    String current_contig = "";
    int commit_count = 0;


    HashMap<String, Event> uncommitted_events = new HashMap<String, Event>();
    List<Event> committed_events = new ArrayList<Event>();


    public void commitEvents(GenomeLoc pos) {
        List<String> just_committed = new ArrayList();

        if (!pos.getContig().equals(current_contig)) {
            //contig switch , committing all
            current_contig = pos.getContig();
            committed_events.addAll(uncommitted_events.values());
            uncommitted_events.clear();
            return;
        }

        for (Map.Entry<String, Event> entry : uncommitted_events.entrySet()) {
            Event e = entry.getValue();
            //go through all of the events, commit the ones past the window
            if (e.getLocation().getStart() < pos.getStart() - normal_window_size) {
                committed_events.add(e);
                just_committed.add(entry.getKey());
            }

        }

        for (String key : just_committed)
            uncommitted_events.remove(key);
    }


    class ClusterCount {
        public SortedMap<Integer, List<GATKSAMRecord>> normal_reads = new TreeMap<Integer, List<GATKSAMRecord>>();
        public SortedMap<Integer, List<GATKSAMRecord>> tumor_reads = new TreeMap<Integer, List<GATKSAMRecord>>();

        public int normal_count = 0;
        public int tumor_count = 0;


        public boolean progressTo(long position, int normalwindowsize, int tumorwindowsize) {
            List<Integer> to_remove = new ArrayList<Integer>();

            long normal_cut = position - normalwindowsize;
            long tumor_cut = position - tumorwindowsize;

            for (Map.Entry<Integer, List<GATKSAMRecord>> entry : normal_reads.entrySet()) {
                if (entry.getKey() < normal_cut) {
                    to_remove.add(entry.getKey());
                    normal_count = normal_count - entry.getValue().size();

                } else //since it's sorted, if we are in the window we can stop looking for things outside the window
                    break;
            }

            for (Integer i : to_remove) {
                normal_reads.remove(i);
            }

            to_remove.clear();

            for (Map.Entry<Integer, List<GATKSAMRecord>> entry : tumor_reads.entrySet()) {
                if (entry.getKey() < tumor_cut) {
                    to_remove.add(entry.getKey());
                    tumor_count = tumor_count - entry.getValue().size();
                } else //since it's sorted, if we are in the window we can stop looking for things outside the window
                    break;
            }

            for (Integer i : to_remove) {
                tumor_reads.remove(i);
            }

            if (tumor_reads.size() == 0) //no somatic abnormality anymore in this cluster, signal 'ready to disregard'
                return false;
            else
                return true;
        }
    }

    HashMap<String, ClusterCount> clusters = new HashMap<String, ClusterCount>();    //TODO: the map can be with key 'GenomeLoc' once the equals() and
    //hashCode() functions are verified to be okay
    ClusterCount normalreads_cluster = new ClusterCount();

    @Override
    public boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context) {

        arguments = argumentCollection;
        normal_evidence = Evidence.get("dna_normal");
        tumor_evidence = Evidence.get("dna_tumor");

        min_event_p = 1.0 - Math.pow(10.0, -argumentCollection.quality / 10.0);

        if (arguments.enable_structvar && normal_evidence != null && tumor_evidence != null)
            return true;
        else
            return false;
    }

    @Override
    public void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts) {
        if (commit_count == 0)
            commitEvents(context.getLocation());

        long pos = context.getPosition();
        //clear reads outside the window
        List<String> removed_cluster_keys = new ArrayList<String>();

        for (Map.Entry<String, ClusterCount> entry : clusters.entrySet()) {
            ClusterCount c = entry.getValue();
            if (!c.progressTo(pos, normal_window_size, tumor_window_size))
                removed_cluster_keys.add(entry.getKey());
        }

        for (String key : removed_cluster_keys)
            clusters.remove(key);

        normalreads_cluster.progressTo(pos, normal_window_size, tumor_window_size);

        int range = 3 * arguments.expected_insert_size;

        StringBuilder region_string_builder = null;

        for (PileupElement p : normal_evidence.Pileup) {
            if (p.getRead().getAlignmentStart() != pos || !p.getRead().getReadPairedFlag())
                continue;


            if (p.getRead().getMateUnmappedFlag())
                continue;

            int insert_size = Math.abs(p.getRead().getInferredInsertSize());

            if (insert_size < range) {
                List<GATKSAMRecord> normal_reads_loc = normalreads_cluster.normal_reads.get(p.getRead().getAlignmentStart());

                if (normal_reads_loc == null)
                    normal_reads_loc = new ArrayList<GATKSAMRecord>();

                normal_reads_loc.add(p.getRead());

                normalreads_cluster.normal_reads.put(p.getRead().getAlignmentStart(), normal_reads_loc);
                normalreads_cluster.normal_count++;
            } else {
                region_string_builder = new StringBuilder();
                region_string_builder.append(p.getRead().getMateReferenceName());
                region_string_builder.append('-');
                region_string_builder.append(p.getRead().getMateAlignmentStart() / cluster_size);

                ClusterCount cluster = clusters.get(region_string_builder.toString());

                if (cluster == null)
                    cluster = new ClusterCount();

                List<GATKSAMRecord> normal_reads_loc = cluster.normal_reads.get(p.getRead().getAlignmentStart());

                if (normal_reads_loc == null)
                    normal_reads_loc = new ArrayList<GATKSAMRecord>();

                normal_reads_loc.add(p.getRead());

                cluster.normal_reads.put(p.getRead().getAlignmentStart(), normal_reads_loc);
                cluster.normal_count++;

                clusters.put(region_string_builder.toString(), cluster);
            }

        }

        for (PileupElement p : tumor_evidence.Pileup) {
            if (p.getRead().getAlignmentStart() != pos)
                continue;

            if (p.getRead().getMateUnmappedFlag())
                continue;

            int insert_size = Math.abs(p.getRead().getInferredInsertSize());
            if (insert_size < range) {
                List<GATKSAMRecord> normal_reads_loc = normalreads_cluster.tumor_reads.get(p.getRead().getAlignmentStart());

                if (normal_reads_loc == null)
                    normal_reads_loc = new ArrayList<GATKSAMRecord>();

                normal_reads_loc.add(p.getRead());

                normalreads_cluster.tumor_reads.put(p.getRead().getAlignmentStart(), normal_reads_loc);
                normalreads_cluster.tumor_count++;
            } else {
                region_string_builder = new StringBuilder();
                region_string_builder.append(p.getRead().getMateReferenceName());
                region_string_builder.append('-');
                region_string_builder.append(p.getRead().getMateAlignmentStart() / cluster_size);

                ClusterCount cluster = clusters.get(region_string_builder.toString());

                if (cluster == null)
                    cluster = new ClusterCount();

                List<GATKSAMRecord> tumor_reads_loc = cluster.normal_reads.get(p.getRead().getAlignmentStart());

                if (tumor_reads_loc == null)
                    tumor_reads_loc = new ArrayList<GATKSAMRecord>();

                tumor_reads_loc.add(p.getRead());

                cluster.tumor_reads.put(p.getRead().getAlignmentStart(), tumor_reads_loc);
                cluster.tumor_count++;

                clusters.put(region_string_builder.toString(), cluster);
            }
        }

        //run test for each cluster
        for (Map.Entry<String, ClusterCount> entry : clusters.entrySet()) {
            ClusterCount cluster = entry.getValue();

            int K1 = cluster.normal_reads.size();
            int N1 = normalreads_cluster.normal_reads.size() + K1;
            int K2 = cluster.tumor_reads.size();
            int N2 = normalreads_cluster.tumor_reads.size() + K2;

            double p_sv = GetSomaticProbability(K1, N1, K2, N2, sv_alpha, sv_beta);

            if (p_sv > min_event_p) {
                String event_name = entry.getKey();

                if (uncommitted_events.get(event_name) == null || uncommitted_events.get(event_name).getProbability() < p_sv) {
                    //we got a better version of the event, or we haven't seen the event before
                    BreakpointInfo br = getBreakPointInfo(cluster);

                    //create the actual event!!
                    Event e = new Event(br.toString(), p_sv, context.getLocation());

                    for (GeneContext gc : ann_contexts) {
                        e.setAttribute(gc.name, null);
                    }

                    uncommitted_events.put(event_name, e);

                }
            }
        }

        if (++commit_count == tumor_window_size)
            commit_count = 0;
    }


    class BreakpointInfo {
        public int bp1_start;
        public int bp1_end;
        public int bp2_start;
        public int bp2_end;

        public String contig1;
        public String contig2;

        @Override
        public String toString() {
            StringBuilder b = new StringBuilder();
            b.append(contig1).append(':');
            b.append(bp1_start);
            b.append('-');
            b.append(bp1_end);
            b.append("<->");
            b.append(contig2).append(':');
            b.append(bp2_start);
            b.append('-');
            b.append(bp2_end);

            return b.toString();
        }

    }

    private BreakpointInfo getBreakPointInfo(ClusterCount cluster) {
        BreakpointInfo bp = new BreakpointInfo();

        //bp1 start is easy
        bp.bp1_start = cluster.tumor_reads.firstKey();

        bp.contig1 = cluster.tumor_reads.get(cluster.tumor_reads.firstKey()).get(0).getReferenceName();
        bp.contig2 = cluster.tumor_reads.get(cluster.tumor_reads.firstKey()).get(0).getMateReferenceName();

        for (List<GATKSAMRecord> readlist : cluster.tumor_reads.values()) {
            for (GATKSAMRecord r : readlist) {
                //bp1 end is the end of the longest-aligned read
                if (bp.bp1_end < r.getAlignmentEnd())
                    bp.bp1_end = r.getAlignmentEnd();

                if (bp.bp2_start == 0 || bp.bp2_start > r.getMateAlignmentStart())
                    bp.bp2_start = r.getMateAlignmentStart();

                //bp2 end is less accurate because we don't have the length of the reads.
                if (bp.bp2_end < r.getMateAlignmentStart())
                    bp.bp2_end = r.getMateAlignmentStart();

            }

        }

        return bp;
    }

    public static double GetSomaticProbability(int K1, int N1, int K2, int N2, double alpha, double beta) {
        double alpha_ref = alpha;
        double beta_ref = beta;

        double prior_nosv = 1 - prior_sv;
        double p_d_nosv = prior_nosv * SeuratMethods.BetaBinomialPdf(N2, K2, K1 + alpha_ref, N1 - K1 + beta_ref);
        double p_d_sv = prior_sv * SeuratMethods.BetaBinomialPdf(N2, K2, 1, 1);

        return (p_d_sv / (p_d_sv + p_d_nosv));

    }

    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return (normal_evidence.Pileup.size() > 5 && tumor_evidence.Pileup.size() > 5);
    }

    @Override
    public boolean isTranscriptionalAnalysis() {
        return false;
    }

    @Override
    public List<Event> reduce() {
        committed_events.addAll(uncommitted_events.values());

        for (Event e : committed_events) {
            //remove locations - they were used as part of the breakpoint evidence discovery process but isn't actually meaningful
            e.setLocation(null);
        }

        return committed_events;
    }
}
