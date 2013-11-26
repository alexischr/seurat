/*
 * Copyright (c) 2012 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/30/11
 * Time: 12:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneContext {
    enum GeneContextClass {
        Nongenic,
        Exon,
        Intron,
        UTR,
        Transcript,
        Unknown
    }

    public GeneContextClass context;
    public String name;
    public String gene_name = "";
    public boolean coding = false;

    public GeneContext(String name, GeneContextClass context) {
        this.name = name;
        this.context = context;
    }

    public GeneContext(String gene_name, String name, GeneContextClass context) {
        this.gene_name = gene_name;
        this.name = name;
        this.context = context;
    }

}
