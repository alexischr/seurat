/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;


public class Event {

    public Event(String name, double probability, GeneContext context) {
        this.name = name;
        this.probability = probability;
        this.context = context;
    }

    public Event(String name, double probability, GenomeLoc location) {
        this.name = name;
        this.probability = probability;
        this.location = location;
        this.context = new GeneContext("unknown", GeneContext.GeneContextClass.Nongenic);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }


    public GeneContext getGeneContext() {
        return context;
    }


    private String name;

    private GenomeLoc location = null;
    private GeneContext context = null;

    private List<Event> childEvents = new ArrayList<Event>();

    private Map<String, Object> attributes = new LinkedHashMap<String, Object>();

    private boolean transcriptional = false;

    public void setAttribute(String key, Object value) {
        attributes.put(key, value);
    }

    public Object getAttribute(String key) {
        return attributes.get(key);
    }

    public String attributeString() {
        StringBuilder s = new StringBuilder();

        for (Map.Entry<String, Object> entry : attributes.entrySet()) {
            if (entry.getValue() != null) {
                s.append(entry.getKey()).append('=');
                s.append(entry.getValue());
                s.append(';');
            } else {
                s.append(entry.getKey());
                s.append(';');
            }
        }

        s.deleteCharAt(s.length() - 1);

        return s.toString();
    }


    public String attributeStringVCF() {
        StringBuilder s = new StringBuilder();

        for (Map.Entry<String, Object> entry : attributes.entrySet()) {
            if (entry.getKey().equals("REF") || entry.getKey().equals("ALT"))
                continue;
            else {
                if (entry.getValue() != null) {
                    s.append(entry.getKey()).append('=');

                    String formatted_value;

                    if (entry.getValue().getClass() == Float.class)
                        formatted_value = String.format("%.3f", entry.getValue());
                    else
                        formatted_value = entry.getValue().toString();


                    s.append(formatted_value);
                    s.append(';');
                } else {
                    s.append(entry.getKey());
                    s.append(';');
                }
            }
        }

        s.deleteCharAt(s.length() - 1);

        return s.toString();
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();

        double phredq = -10 * Math.log10(1 - probability);

        if (phredq > 255)
            phredq = 255; //TODO remove hardcoded values
        s.append(name);
        s.append('\t');
        s.append(String.format("%.1f", phredq));
        s.append('\t');
        if (location != null) {
            s.append(location.toString());
            s.append('\t');
        }

        s.append(attributeString());

        /*for (Event e : childEvents) {
            s.append(e.toString());
            s.append('\t');
        } */

        return s.toString();
    }

    public void addChildEvent(Event childEvent) {
        childEvents.add(childEvent);
    }

    public List<Event> getChildEvents() {
        return childEvents;
    }

    public GenomeLoc getLocation() {
        return location;
    }

    public void setLocation(GenomeLoc location) {
        this.location = location;
    }

    public double getProbability() {
        return probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

    public double probability;

}
