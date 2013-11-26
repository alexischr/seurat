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
public class LOHAnalysis extends RegionalAnalysisWalker {


    SeuratArgumentCollection arguments;
    double min_event_p;
    List<Event> temp_events = new ArrayList<Event>();

    PileupEvidence normal_evidence;
    PileupEvidence tumor_evidence;

    //priors
    double prior_hom = 0.985;
    double prior_homvar = 0.005;
    double prior_het = 0.01;

    double prior_loh = 0.01;

    @Override
    public boolean initialize(Map<String, PileupEvidence> Evidence, SeuratArgumentCollection argumentCollection, GeneContext context) {
        arguments = argumentCollection;
        normal_evidence = Evidence.get("dna_normal");
        tumor_evidence = Evidence.get("dna_tumor");

        min_event_p = 1.0 - Math.pow(10.0, -argumentCollection.quality / 10.0);

        if (arguments.enable_loh && normal_evidence != null && tumor_evidence != null)
            return true;
        else
            return false;
    }

    public double GetHetProbability(int K, int N, double alpha, double beta) {
        double alpha_ref = alpha;
        double beta_ref = beta;

        double alpha_nonref = beta;
        double beta_nonref = alpha;

        double p_d_hom = prior_hom * SeuratMethods.BetaBinomialPdf(N, K, alpha_ref, beta_ref);
        double p_d_het = prior_het * SeuratMethods.BetaBinomialPdf(N, K, 1, 1);
        double p_d_nonref = prior_homvar * SeuratMethods.BetaBinomialPdf(N, K, alpha_nonref, beta_nonref);

        double p = p_d_het / (p_d_hom + p_d_het + p_d_nonref);
        return p;

    }

    public double GetLOHProbability(int K1, int N1, int K2, int N2, double alpha, double beta, boolean refhom_only) {
        double alpha_ref = alpha;
        double beta_ref = beta;

        double alpha_nonref = beta;
        double beta_nonref = alpha;

        double prior_noloh = 1 - (2 * prior_loh);

        double p_d_loh = prior_loh * SeuratMethods.BetaBinomialPdf(N2, K2, alpha_ref, beta_ref);
        double p_d_nonref_loh = prior_loh * SeuratMethods.BetaBinomialPdf(N2, K2, alpha_nonref, beta_nonref);
        double p_d_noloh = prior_noloh * SeuratMethods.BetaBinomialPdf(N2, K2, K1 + 1, N1 - K1 + 1);


        double p_nonref_loh = (p_d_loh + p_d_nonref_loh) / (p_d_loh + p_d_nonref_loh + p_d_noloh);

        return GetHetProbability(K1, N1, alpha, beta) * p_nonref_loh;
    }


    @Override
    public void map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, List<GeneContext> ann_contexts) {
        Byte alt_loh = 'N';
        double max_p_loh = 0;

        int refcount_one = normal_evidence.ReferencePileup.size();
        int refcount_two = tumor_evidence.ReferencePileup.size();

        if (refcount_one < 5 && refcount_two < 5)
            return;

        for (Map.Entry<Byte, List<PileupElement>> nonref_pileup : normal_evidence.NonReferencePileups.entrySet()) {
            Byte c_alt = nonref_pileup.getKey();

            if (c_alt == 'I' || c_alt == 'D')
                continue;

            List<PileupElement> two_nonref_pileup = tumor_evidence.NonReferencePileups.get(c_alt);

            int K1 = nonref_pileup.getValue().size();
            int N1 = normal_evidence.Pileup.size();

            int K2 = two_nonref_pileup == null ? 0 : two_nonref_pileup.size();
            int N2 = tumor_evidence.Pileup.size();

            if (N1 >= 5 && N2 >= 5) {
                //double p_loh = (1 - GetReferenceProbability(K1, N1, 1, N1 + 1)) * GetReferenceProbability(K2, N2, 1, N2 + 1)
                //      + (1 - GetNonReferenceProbability(K1, N1, N1 + 1, 1)) * GetNonReferenceProbability(K2, N2, N2 + 1, 1);

                double p_loh = GetLOHProbability(K1, N1, K2, N2, arguments.beta_alpha, arguments.beta_beta, false);
                //double p_loh = (1 - GetReferenceProbability(K1, N1, arguments.beta_alpha, arguments.beta_beta)) * GetReferenceProbability(K2, N2, arguments.beta_alpha, arguments.beta_beta);
                /*(1 - GetNonReferenceProbability(K1, N1, arguments.beta_beta, arguments.beta_alpha)) * GetNonReferenceProbability(K2, N2, arguments.beta_beta, arguments.beta_alpha); */

                if (p_loh > max_p_loh) {
                    max_p_loh = p_loh;
                    alt_loh = c_alt;
                }
            }
        }

        if (max_p_loh > min_event_p) {
            Event new_loh = new Event("somatic_LOH", max_p_loh, context.getLocation());
            new_loh.setAttribute("REF", String.format("%c", ref.getBase()));
            new_loh.setAttribute("DP1", normal_evidence.Pileup.size());
            new_loh.setAttribute("DP2", tumor_evidence.Pileup.size());
            new_loh.setAttribute("ALT", new String(new byte[]{alt_loh}));

            if (arguments.pileup_info) {
                new_loh.setAttribute("PILEUP", normal_evidence.toString() + '/' + tumor_evidence.toString());
            }

            temp_events.add(new_loh);

        }
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
        return temp_events;
    }
}
