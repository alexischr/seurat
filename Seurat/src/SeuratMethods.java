/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.apache.commons.math.special.Beta;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/28/11
 * Time: 5:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class SeuratMethods {

    public static double BetaBinomialPdf(int N, int K, double alpha, double beta) {

        //return Math.exp(Beta.logBeta(K + alpha, N - K + beta)) / Math.exp(Beta.logBeta(alpha, beta));    this is needlessly high accuracy, and slower
        return Math.exp(Beta.logBeta(K + alpha, N - K + beta) - Beta.logBeta(alpha, beta));
    }


}
