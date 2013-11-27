/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Alexis
 * Date: 12/17/10
 * Time: 4:10 AM
 * To change this template use File | Settings | File Templates.
 */
public class PileupEvidence {

    public List<PileupElement> Pileup, ReferencePileup;
    public Map<Byte, List<PileupElement>> NonReferencePileups, AllelePileups;


    public static final Byte DELETION_CHAR = 'D';
    public static final Byte INSERTION_CHAR = 'I';

    final public Byte[] classes = {'A', 'C', 'G', 'T', DELETION_CHAR, INSERTION_CHAR};


    public int VariantCount;
    public int RefCount;

    public byte RefBase;

    public double GetAllelicRatio() {
        if (Pileup.size() == 0)
            return 0;
        else return (double) VariantCount / Pileup.size();
    }


    public PileupEvidence() {
        Pileup = new ArrayList<PileupElement>();
        NonReferencePileups = new HashMap<Byte, List<PileupElement>>();
        ReferencePileup = new ArrayList<PileupElement>();
        AllelePileups = new HashMap<Byte, List<PileupElement>>();

        //add one allele pileup per class

        for (Byte class_char : classes) {
            AllelePileups.put(class_char, new ArrayList<PileupElement>());
        }

        VariantCount = 0;
        RefCount = 0;
    }

    public void Reset(byte ref_base) {
        Pileup.clear();
        ReferencePileup.clear();

        for (List<PileupElement> l : NonReferencePileups.values()) {
            l.clear();
        }

        NonReferencePileups.clear();
        VariantCount = 0;
        RefCount = 0;
        RefBase = ref_base;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (PileupElement p : Pileup) {
            byte base = p.getBase();

            if (p.isDeletion())
                base = DELETION_CHAR;

            if (p instanceof ExtendedEventPileupElement) {
                if (((ExtendedEventPileupElement) p).isInsertion())
                    base = INSERTION_CHAR;
            }

            if (p.getRead().getReadNegativeStrandFlag())
                s.append(Character.toLowerCase((char) base));
            else
                s.append((char) base);
        }
        return s.toString();
    }

    public void Add(PileupElement p) {

        Byte base = p.getBase();

        if (p.isDeletion())
            base = DELETION_CHAR;

        if (p.isDeletion() && p instanceof ExtendedEventPileupElement)
            return;

        if (p instanceof ExtendedEventPileupElement) {
            if (((ExtendedEventPileupElement) p).isInsertion()) {
                base = INSERTION_CHAR;
            } else { //we don't know how to handle any other events at the moment!!

            }
        }

        Pileup.add(p);

        if (base == RefBase) {
            RefCount++;
            ReferencePileup.add(p);
        } else {
            VariantCount++;

            List<PileupElement> nrl = NonReferencePileups.get(base);
            if (nrl == null)
                nrl = AllelePileups.get(base);

            if (nrl == null)
                return;          //not one of the 'real' alleles (e.g a code , N or M), disregard

            nrl.add(p);
            NonReferencePileups.put(base, nrl);
        }
    }

}
