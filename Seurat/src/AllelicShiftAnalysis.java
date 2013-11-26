/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/21/11
 * Time: 6:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class AllelicShiftAnalysis extends RegionalAnalysisWalker {

    SeuratArgumentCollection arguments;
    double min_event_p;
    List<Event> eventList = new ArrayList<Event>();

    PileupEvidence rna_normal_evidence;
    PileupEvidence rna_tumor_evidence;

    static final private double p_gene_ab_odds_prior = 10000;
    private double p_gene_ab_odds;

    GeneContext geneContext;

    @Override
    public boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context) {
        arguments = argumentCollection;

        if (Evidence.get("rna_normal") != null)
            rna_normal_evidence = Evidence.get("rna_normal");
        else
            rna_normal_evidence = Evidence.get("dna_tumor");

        rna_tumor_evidence = Evidence.get("rna_tumor");

        geneContext = context;

        min_event_p = 1.0 - Math.pow(10.0, -argumentCollection.quality / 10.0);
        p_gene_ab_odds = p_gene_ab_odds_prior;

        if (rna_normal_evidence != null && rna_tumor_evidence != null && context.context == GeneContext.GeneContextClass.Transcript)
            return true;
        else
            return false;
    }

    public double GetHetProbability(int K, int N, double alpha, double beta) {
        return 1 - GetReferenceProbability(K, N, alpha, beta) - GetNonReferenceProbability(K, N, beta, alpha);
    }

    public double GetReferenceProbability(int K, int N, double alpha, double beta) {
        double alpha_ref = alpha;
        double beta_ref = beta;

        double alpha_nonref = beta;
        double beta_nonref = alpha;

        double prior_hom = 0.985;
        double prior_homnonref = 0.005;
        double prior_het = 0.01;

        double p_d_hom = prior_hom * SeuratMethods.BetaBinomialPdf(N, K, alpha_ref, beta_ref);
        double p_d_het = prior_het * SeuratMethods.BetaBinomialPdf(N, K, 1, 1);
        double p_d_nonref = prior_homnonref * SeuratMethods.BetaBinomialPdf(N, K, alpha_nonref, beta_nonref);

        double p = p_d_hom / (p_d_hom + p_d_het + p_d_nonref);
        return p;
    }

    public double GetNonReferenceProbability(int K, int N, double alpha, double beta) {
        double alpha_ref = beta;
        double beta_ref = alpha;

        double alpha_nonref = alpha;
        double beta_nonref = beta;

        double prior_hom = 0.985;
        double prior_homnonref = 0.005;
        double prior_het = 0.01;

        double p_d_hom = prior_hom * SeuratMethods.BetaBinomialPdf(N, K, alpha_ref, beta_ref);
        double p_d_het = prior_het * SeuratMethods.BetaBinomialPdf(N, K, 1, 1);
        double p_d_nonref = prior_homnonref * SeuratMethods.BetaBinomialPdf(N, K, alpha_nonref, beta_nonref);

        double p = p_d_nonref / (p_d_hom + p_d_het + p_d_nonref);
        return p;
    }

    public double getAlellicBalanceOdds(int N1, int K1, int N2, int K2) {
        double alpha_m1 = K1 + 1;
        double beta_m1 = N1 - K1 + 1;

        double alpha_m2 = 1;
        double beta_m2 = 1;

        double log_likelihood_m1 = Math.log(SeuratMethods.BetaBinomialPdf(N2, K2, alpha_m1, beta_m1));
        //double het_p = GetHetProbability(K1, N1, alpha, beta);
        double het_p = 0.01;
        double log_likelihood_m2 = Math.log(het_p * SeuratMethods.BetaBinomialPdf(N2, K2, alpha_m2, beta_m2)
                + (1 - het_p) * SeuratMethods.BetaBinomialPdf(N2, K2, alpha_m1, beta_m1));

        return Math.exp(log_likelihood_m1 - log_likelihood_m2);
    }


    @Override
    public void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts) {

        int refcount_one = rna_normal_evidence.ReferencePileup.size();
        int refcount_two = rna_tumor_evidence.ReferencePileup.size();

        for (Map.Entry<Byte, List<PileupElement>> nonref_pileup : rna_normal_evidence.NonReferencePileups.entrySet()) {
            Byte c_alt = nonref_pileup.getKey();

            if (c_alt == 'I' || c_alt == 'D')
                continue;

            List<PileupElement> two_nonref_pileup = rna_tumor_evidence.NonReferencePileups.get(c_alt);

            int K1 = nonref_pileup.getValue().size();
            int N1 = refcount_one + K1;

            int K2 = two_nonref_pileup == null ? 0 : two_nonref_pileup.size();
            int N2 = refcount_two + K2;

            double p1 = ((double) K1) / N1;
            double p2 = ((double) K2) / N2;

            if (N1 >= 5 && N2 >= 5 && p2 > p1) {
                double p_ab_odds = getAlellicBalanceOdds(N1, K1, N2, K2);
                p_gene_ab_odds = p_gene_ab_odds * p_ab_odds;

                double prob_ab = 1 - (p_ab_odds / (1 + p_ab_odds));

                //if (prob_ab > min_event_p) //strong Evidence for allelic imbalance
                //{
                Event new_ail = new Event("allelic_imbalance_locus", prob_ab, context.getLocation());
                new_ail.setAttribute("REF", String.format("%c", ref.getBase()));
                new_ail.setAttribute("DP1", rna_normal_evidence.Pileup.size());
                new_ail.setAttribute("DP2", rna_tumor_evidence.Pileup.size());

                new_ail.setAttribute("ALT", new String(new byte[]{c_alt}));

                if (arguments.pileup_info) {
                    new_ail.setAttribute("PILEUP", rna_normal_evidence.toString() + '/' + rna_tumor_evidence.toString());
                }

                eventList.add(new_ail);
                //}
            }
        }


    }

    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return (rna_normal_evidence.Pileup.size() > 5 && rna_tumor_evidence.Pileup.size() > 5);
    }

    @Override
    public boolean isTranscriptionalAnalysis() {
        return true;
    }

    @Override
    public List<Event> reduce() {
        //emit event if detected
        Event allelic_loss_event;
        double p_noal = (p_gene_ab_odds / (1 + p_gene_ab_odds));
        allelic_loss_event = new Event("allelic_imbalance", 1 - p_noal, geneContext);
        eventList.add(allelic_loss_event);
        return eventList;
    }
}
