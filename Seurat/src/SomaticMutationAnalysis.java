/*
 * Copyright (c) 2014 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/21/11
 * Time:12:13PM
 * To change this template use File|Settings|File Templates.
 */

public class SomaticMutationAnalysis extends RegionalAnalysisWalker {

    SeuratArgumentCollection arguments;
    double min_event_p;
    List<Event> temp_events = new ArrayList<Event>(50000);

    PileupEvidence normal_evidence;
    PileupEvidence tumor_evidence;

    PileupEvidence rna_normal_evidence;
    PileupEvidence rna_tumor_evidence;

    //priors
    double prior_hom = 0.9985;
    double prior_homvar = 0.0005;
    double prior_het = 0.001;

    double prior_snv = 0.01;
    double prior_sindel = 0.000001;

    public boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context) {
        arguments = argumentCollection;
        normal_evidence = Evidence.get("dna_normal");
        tumor_evidence = Evidence.get("dna_tumor");
        rna_normal_evidence = Evidence.get("rna_normal");
        rna_tumor_evidence = Evidence.get("rna_tumor");

        min_event_p = 1.0 - Math.pow(10.0, -argumentCollection.quality / 10.0);
        prior_snv = arguments.p_mutation;

        return (normal_evidence != null && tumor_evidence != null);
    }

    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return (normal_evidence.Pileup.size() > 5 && tumor_evidence.Pileup.size() > 5);
    }

    @Override
    public boolean isTranscriptionalAnalysis() {
        return false;
    }

    public double GetReferenceProbability(int K, int N, double alpha, double beta) {
        double alpha_ref = alpha;
        double beta_ref = beta;

        double alpha_nonref = beta;
        double beta_nonref = alpha;

        double p_d_hom = prior_hom * SeuratMethods.BetaBinomialPdf(N, K, alpha_ref, beta_ref);
        double p_d_het = prior_het * SeuratMethods.BetaBinomialPdf(N, K, 1, 1);
        double p_d_nonref = prior_homvar * SeuratMethods.BetaBinomialPdf(N, K, alpha_nonref, beta_nonref);

        double p = p_d_hom / (p_d_hom + p_d_het + p_d_nonref);
        return p;
    }

    public double GetNonReferenceProbability(int K, int N, double alpha, double beta) {
        double alpha_ref = beta;
        double beta_ref = alpha;

        double alpha_nonref = alpha;
        double beta_nonref = beta;

        double p_d_hom = prior_hom * SeuratMethods.BetaBinomialPdf(N, K, alpha_ref, beta_ref);
        double p_d_het = prior_het * SeuratMethods.BetaBinomialPdf(N, K, 1, 1);
        double p_d_nonref = prior_homvar * SeuratMethods.BetaBinomialPdf(N, K, alpha_nonref, beta_nonref);

        double p = p_d_nonref / (p_d_hom + p_d_het + p_d_nonref);
        return p;
    }


    public double GetSomaticProbability(int K1, int N1, int K2, int N2, double alpha, double beta, boolean refhom_only) {
        double alpha_ref = alpha;
        double beta_ref = beta;

        double alpha_nonref = beta;
        double beta_nonref = alpha;

        double prior_nosnv = 1 - prior_snv;
        double p_d_refnosnv = prior_nosnv * SeuratMethods.BetaBinomialPdf(N2, K2, K1 + alpha_ref, N1 - K1 + beta_ref);
        double p_d_snv = prior_snv * SeuratMethods.BetaBinomialPdf(N2, K2, 1, 1);

        double p_ref_snv = GetReferenceProbability(K1, N1, alpha_ref, beta_ref) * (p_d_snv / (p_d_snv + p_d_refnosnv));

        if (refhom_only)
            return p_ref_snv;
        else {
            double p_d_nonrefnosnv = prior_nosnv * SeuratMethods.BetaBinomialPdf(N2, K2, K1 + alpha_nonref, N1 - K1 + beta_nonref);
            double p_nonref_snv = GetNonReferenceProbability(K1, N1, alpha_nonref, beta_nonref) * (p_d_snv / (p_d_snv + p_d_nonrefnosnv));
            return p_ref_snv + p_nonref_snv;
        }
    }

    private class SomaticCall {
        public double p_snv = 0;
        public Byte alt_snv;
        public String ref_seq = "";
        public String alt_seq;
        public int K1 = 0, N1 = 0, K2 = 0, N2 = 0;
        public int rK1 = 0, rN1 = 0, rK2 = 0, rN2 = 0;

        public int dna_alt_allele_forward = 0, dna_alt_allele_reverse = 0, dna_alt_allele_total = 0;
        public double dna_alt_allele_forward_fraction = 0, dna_alt_allele_reverse_fraction = 0, dna_alt_allele_total_fraction = 0;
        public int dna_ref_allele_forward = 0, dna_ref_allele_reverse = 0, dna_ref_allele_total = 0;

        public SomaticCall(Byte alt_snv, String alt_seq) {
            this.alt_snv = alt_snv;
            this.alt_seq = alt_seq;
        }

        public SomaticCall(String alt_seq) {
            this.alt_snv = 'N';
            this.alt_seq = alt_seq;
        }
    }

    @Override
    public void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts) {
        GenomeLoc loc = context.getLocation();

        SomaticCall best_call = new SomaticCall("");
        best_call.ref_seq = String.format("%c", ref.getBase());


        for (Map.Entry<Byte, List<PileupElement>> tumor_var_pileup : tumor_evidence.NonReferencePileups.entrySet()) {
            SomaticCall current_call = new SomaticCall(tumor_var_pileup.getKey(), "");

            List<PileupElement> normal_var_pileup = normal_evidence.NonReferencePileups.get(current_call.alt_snv);

            current_call.K1 = normal_var_pileup == null ? 0 : normal_var_pileup.size();
            current_call.N1 = normal_evidence.RefCount + current_call.K1;

            current_call.K2 = tumor_var_pileup.getValue().size();
            current_call.N2 = tumor_evidence.RefCount + current_call.K2;

            if (rna_normal_evidence != null) {
                List<PileupElement> rna_normal_var_pileup = rna_normal_evidence.NonReferencePileups.get(current_call.alt_snv);
                current_call.rK1 = rna_normal_var_pileup == null ? 0 : rna_normal_var_pileup.size();
                current_call.rN1 = rna_normal_evidence.RefCount + current_call.rK1;

            }
            if (rna_tumor_evidence != null) {
                List<PileupElement> rna_tumor_var_pileup = rna_tumor_evidence.NonReferencePileups.get(current_call.alt_snv);
                current_call.rK2 = rna_tumor_var_pileup == null ? 0 : rna_tumor_var_pileup.size();
                current_call.rN2 = rna_tumor_evidence.RefCount + current_call.rK2;

            }


            if (current_call.alt_snv.equals(PileupEvidence.INSERTION_CHAR)) {
                HashMap<String, Integer> hs = new HashMap<String, Integer>();

                for (PileupElement p : tumor_var_pileup.getValue()) {
                    String seq = p.getEventBases();

                    if (hs.containsKey(seq)) {
                        hs.put(seq, hs.get(seq) + 1);
                    } else {
                        hs.put(seq, 1);
                    }

                }

                int cons_n = 0;
                for (Map.Entry<String, Integer> i : hs.entrySet()) {
                    if (i.getValue() > cons_n) {
                        cons_n = i.getValue();
                        current_call.alt_seq = i.getKey();
                    }
                }

            }
            if (current_call.alt_snv.equals(PileupEvidence.DELETION_CHAR)) {
                current_call.alt_seq = String.format("%c", ref.getBases()[99]);

            }

            if (arguments.merge_rna) {
                current_call.K1 += current_call.rK1;
                current_call.N1 += current_call.rN1;
                current_call.K2 += current_call.rK2;
                current_call.N2 += current_call.rN2;
            }


            /// call! ///
            current_call.p_snv = GetSomaticProbability(current_call.K1, current_call.N1, current_call.K2, current_call.N2, arguments.beta_alpha, arguments.beta_beta, arguments.refnormal_only);
            ///


            if (arguments.enable_debug) {
                System.err.printf("N=%d/%d -- T=%d/%d --- p_snv %f\n", current_call.K1, current_call.N1, current_call.K2, current_call.N2, current_call.p_snv);
            }

            boolean call_approved = true;

            do {
                if (current_call.p_snv < best_call.p_snv) {
                    call_approved = false;
                    break;
                }

                if (arguments.both_strands) {
                    boolean pos_strand = false;
                    boolean neg_strand = false;
                    for (PileupElement p : tumor_var_pileup.getValue()) {
                        if (p.getRead().getReadNegativeStrandFlag())
                            neg_strand = true;
                        else
                            pos_strand = true;
                        if (neg_strand && pos_strand) {
                            break;
                        }
                    }
                    if (!(neg_strand && pos_strand)) {
                        call_approved = false;
                        break;
                    }
                }

                if (call_approved) {
                    best_call = current_call;

                    //ar1 = ((float) K1) / N1;
                    //ar2 = ((float) K2) / N2;

                    if (arguments.enable_allele_metrics) {
                        for (PileupElement p : tumor_var_pileup.getValue()) {
                            if (p.getRead().getReadNegativeStrandFlag())
                                best_call.dna_alt_allele_forward++;
                            else
                                best_call.dna_alt_allele_reverse++;
                        }
                        best_call.dna_alt_allele_total = tumor_var_pileup.getValue().size();

                        best_call.dna_alt_allele_forward_fraction = ((double) best_call.dna_alt_allele_forward) / best_call.dna_alt_allele_total;
                        best_call.dna_alt_allele_reverse_fraction = ((double) best_call.dna_alt_allele_reverse) / best_call.dna_alt_allele_total;
                        best_call.dna_alt_allele_total_fraction = ((double) best_call.dna_alt_allele_total) / best_call.N2;

                        for (PileupElement p : tumor_evidence.ReferencePileup) {
                            if (p.getRead().getReadNegativeStrandFlag())
                                best_call.dna_ref_allele_forward++;
                            else
                                best_call.dna_ref_allele_reverse++;
                        }

                        best_call.dna_ref_allele_total = tumor_evidence.ReferencePileup.size();
                    }
                }
            } while (false);

        }

        if (best_call.p_snv > min_event_p) {   //event is over qual threshold

            String event_name;
            String alt;

            if (best_call.alt_snv.equals((PileupEvidence.DELETION_CHAR))) {
                event_name = "somatic_deletion";
                alt = best_call.alt_seq;
                loc = ref.getGenomeLocParser().createGenomeLoc(loc.getContig(), loc.getStart() - 1, loc.getStop() - 1);
                best_call.ref_seq = String.format("%c%c", ref.getBases()[99], ref.getBases()[100]);
            } else if (best_call.alt_snv.equals(PileupEvidence.INSERTION_CHAR)) {
                event_name = "somatic_insertion";
                alt = String.format("%c%s", normal_evidence.RefBase, best_call.alt_seq);
            } else {
                event_name = "somatic_SNV";
                alt = new String(new byte[]{best_call.alt_snv});
            }

            Event new_snv = new Event(event_name, best_call.p_snv, loc);

            new_snv.setAttribute("REF", best_call.ref_seq);
            new_snv.setAttribute("DP1", normal_evidence.Pileup.size());
            new_snv.setAttribute("DP2", tumor_evidence.Pileup.size());
            new_snv.setAttribute("AR1", ((float) best_call.K1) / best_call.N1);
            new_snv.setAttribute("AR2", ((float) best_call.K2) / best_call.N2);
            new_snv.setAttribute("LN", 1);

            new_snv.setAttribute("ALT", alt);

            //TODO: proper representation of insertions in ALT field
            if (best_call.alt_snv.equals(PileupEvidence.INSERTION_CHAR)) {
                new_snv.setAttribute("LN", best_call.alt_seq.length());
                new_snv.setAttribute("SEQ", best_call.alt_seq);
            }


            if (arguments.pileup_info) {
                new_snv.setAttribute("PILEUP1", normal_evidence.toString());
                new_snv.setAttribute("PILEUP2", tumor_evidence.toString());
            }

            if (arguments.enable_metrics) {

                EvidenceMetrics normal_metrics = new EvidenceMetrics(normal_evidence);
                EvidenceMetrics tumor_metrics = new EvidenceMetrics(tumor_evidence);

                normal_metrics.GenerateMetrics();
                tumor_metrics.GenerateMetrics();

                new_snv.setAttribute("MVC1", normal_metrics.VC);
                new_snv.setAttribute("MVBQ1", normal_metrics.VBQ);
                new_snv.setAttribute("MVMQ1", normal_metrics.VMQ);
                new_snv.setAttribute("MVC2", tumor_metrics.VC);
                new_snv.setAttribute("MVBQ2", tumor_metrics.VBQ);
                new_snv.setAttribute("MVMQ2", tumor_metrics.VMQ);
            }

            if (arguments.enable_allele_metrics) {
                new_snv.setAttribute("DNA_ALT_ALLELE_FORWARD", best_call.dna_alt_allele_forward);
                new_snv.setAttribute("DNA_ALT_ALLELE_FORWARD_FRACTION", best_call.dna_alt_allele_forward_fraction);
                new_snv.setAttribute("DNA_ALT_ALLELE_REVERSE", best_call.dna_alt_allele_reverse);
                new_snv.setAttribute("DNA_ALT_ALLELE_REVERSE_FRACTION", best_call.dna_alt_allele_reverse_fraction);
                new_snv.setAttribute("DNA_ALT_ALLELE_TOTAL", best_call.dna_alt_allele_total);
                new_snv.setAttribute("DNA_ALT_ALLELE_TOTAL_FRACTION", best_call.dna_alt_allele_total_fraction);
                new_snv.setAttribute("DNA_REF_ALLELE_FORWARD", best_call.dna_ref_allele_forward);
                new_snv.setAttribute("DNA_REF_ALLELE_REVERSE", best_call.dna_ref_allele_reverse);
                new_snv.setAttribute("DNA_REF_ALLELE_TOTAL", best_call.dna_ref_allele_total);
            }

            if (arguments.rna_call) {
                double rna_call;

                if (rna_normal_evidence != null) {
                    rna_call = GetSomaticProbability(best_call.rK1, best_call.rN1, best_call.rK2, best_call.rN2, arguments.beta_alpha, arguments.beta_beta, arguments.refnormal_only);
                } else {
                    rna_call = GetSomaticProbability(best_call.K1, best_call.N1, best_call.rK2, best_call.rN2, arguments.beta_alpha, arguments.beta_beta, arguments.refnormal_only);
                }

                new_snv.setAttribute("RNA_CALL", (float) (-10 * Math.log10(1 - rna_call)));
            }

            temp_events.add(new_snv);
        }

    }

    private static class EvidenceMetrics {
        PileupEvidence evidence;

        public double VBQ;
        public double VMQ;
        public double VC;

        public EvidenceMetrics(PileupEvidence evidence) {
            this.evidence = evidence;
        }


        public void GenerateMetrics() {
            //metrics helper objects
            DescriptiveStatistics VC_stats = new DescriptiveStatistics();
            DescriptiveStatistics VBQ_stats = new DescriptiveStatistics();
            DescriptiveStatistics VMQ_stats = new DescriptiveStatistics();

            ArrayList<Double> relative_cycles = new ArrayList<Double>();


            for (PileupElement element : evidence.Pileup) {

                VBQ_stats.addValue(element.getQual());
                VMQ_stats.addValue(element.getMappingQual());

                double average_cycle = (double) element.getRead().getReadLength() / 2.0;
                double relative_cycle = average_cycle - Math.abs(average_cycle - element.getOffset());

                VC_stats.addValue(relative_cycle);
            }

            VC = VC_stats.getPercentile(50);
            VBQ = VBQ_stats.getPercentile(50);
            VMQ = VMQ_stats.getPercentile(50);


        }


    }

    @Override
    public List<Event> reduce() {
        List<Event> final_list = new ArrayList<Event>();

        Event mergable_event = null;

        for (Event e : temp_events) {
            if (e.getName().equals("somatic_deletion")) {
                if (mergable_event == null)
                    mergable_event = e;
                else {
                    if (mergable_event.getName().equals(e.getName()) && mergable_event.getLocation().getStop() + 1 == e.getLocation().getStart()) {
                        if (mergable_event.getChildEvents().size() == 0) {
                            Event child = mergable_event;
                            mergable_event = new Event(child.getName(), child.getProbability(), child.getLocation());
                            mergable_event.addChildEvent(child);
                            mergable_event.setAttribute("ALT", child.getAttribute("ALT"));
                            mergable_event.setAttribute("REF", child.getAttribute("REF"));
                            mergable_event.setAttribute("LN", 1);
                        }

                        mergable_event.addChildEvent(e);
                        mergable_event.setProbability((mergable_event.getProbability() + e.getProbability()) / 2);
                        mergable_event.setLocation(mergable_event.getLocation().merge(e.getLocation()));

                        mergable_event.setAttribute("REF", mergable_event.getAttribute("REF").toString() + e.getAttribute("REF").toString().substring(1, 2));

                        mergable_event.setAttribute("LN", mergable_event.getLocation().size());
                    } else {
                        final_list.add(mergable_event);
                        mergable_event = null;
                    }
                }
            } else {
                if (mergable_event != null) {
                    final_list.add(mergable_event);
                    mergable_event = null;
                }
                final_list.add(e);
            }


        }
        if (mergable_event != null) {
            final_list.add(mergable_event);
        }

        return final_list;
    }

}
