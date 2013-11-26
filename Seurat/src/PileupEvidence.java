/*
 * Copyright (c) 2012 by The Translational Genomics Research Institute.
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

    public List<PileupElement> Pileup;
    public List<PileupElement> ReferencePileup;
    public Map<Byte, List<PileupElement>> NonReferencePileups;

    public Map<Byte, List<PileupElement>> AllelePileups;

    //public Set<Byte> NonReferenceBases;
    public int NonrefCount;
    public byte RefBase;

    public double GetAllelicRatio() {
        if (Pileup.size() == 0)
            return 0;
        else return (double) NonrefCount / Pileup.size();
    }


    public PileupEvidence() {
        Pileup = new ArrayList<PileupElement>();
        NonReferencePileups = new HashMap<Byte, List<PileupElement>>();
        ReferencePileup = new ArrayList<PileupElement>();
        AllelePileups = new HashMap<Byte, List<PileupElement>>();

        //add one allele pileup per nucleotide
        AllelePileups.put((byte) 'A', new ArrayList<PileupElement>());
        AllelePileups.put((byte) 'C', new ArrayList<PileupElement>());
        AllelePileups.put((byte) 'G', new ArrayList<PileupElement>());
        AllelePileups.put((byte) 'T', new ArrayList<PileupElement>());

        //extended pileup 'nucleotides'
        AllelePileups.put((byte) 'D', new ArrayList<PileupElement>());
        AllelePileups.put((byte) 'I', new ArrayList<PileupElement>());

        NonrefCount = 0;
    }

    public void Reset(byte ref_base) {
        Pileup.clear();
        ReferencePileup.clear();

        for (List<PileupElement> l : NonReferencePileups.values()) {
            l.clear();
        }

        NonReferencePileups.clear();
        NonrefCount = 0;
        RefBase = ref_base;
    }

    @Override
    public String toString() {
        StringBuffer s = new StringBuffer();
        for (PileupElement p : Pileup) {
            byte base = p.getBase();

            if (p.isDeletion())
                base = (byte) 'D';

            if (p instanceof ExtendedEventPileupElement) {
                if (((ExtendedEventPileupElement) p).isInsertion())
                    base = (byte) 'I';
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
            base = (byte) 'D';

        if (p.isDeletion() && p instanceof ExtendedEventPileupElement)
            return;

        if (p instanceof ExtendedEventPileupElement) {
            if (((ExtendedEventPileupElement) p).isInsertion()) {
                base = (byte) 'I';
            } else { //we don't know how to handle any other events at the moment!!

            }
        }

        Pileup.add(p);


        if (base == RefBase)
            ReferencePileup.add(p);
        else {
            List<PileupElement> nrl = NonReferencePileups.get(base);
            if (nrl == null)
                nrl = AllelePileups.get(base);

            if (nrl == null)
                return;          //not one of the 'real' alleles (e.g a code , N or M), disregard

            nrl.add(p);

            NonReferencePileups.put(base, nrl);

            NonrefCount++;
        }
    }

}
