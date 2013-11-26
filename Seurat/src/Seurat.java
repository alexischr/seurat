/*
 * Copyright (c) 2013 by The Translational Genomics Research Institute.
 */

package org.broadinstitute.sting.gatk.walkers.tgen;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.SeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.codecs.refseq.RefSeqCodec;
import org.broadinstitute.sting.utils.codecs.refseq.Transcript;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * A sequence analysis program for somatic mutation discovery
 * and allelic imbalance in paired tumor and normal genome and transcriptome
 * data.
 */
@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@By(DataSource.REFERENCE)
@Reference(window = @Window(start = -100, stop = 100))
public class Seurat extends LocusWalker<Integer, Long> {

    @Output
    private PrintStream out;

    @Argument(fullName = "go", shortName = "go", doc = "Large event output file")
    private String gene_out;

    @ArgumentCollection
    private SeuratArgumentCollection arguments = new SeuratArgumentCollection();

    // control the various parameters to be used
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling.", required = false)
    private int MIN_BASE_QUALITY_SCORE = 10;

    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for calling.", required = false)
    private int MIN_MAPPING_QUALITY_SCORE = 10;

    @Argument(fullName = "min_coverage", shortName = "mcv", doc = "Minimum total coverage for a locus to be used.", required = false)
    private int MIN_COVERAGE = 6;

    @Argument(fullName = "output_all", shortName = "all", doc = "Output all loci with data.", required = false)
    public Boolean OUTPUT_ALL = false;

    @Argument(fullName = "refseq", shortName = "refseq",
            doc = " Name of RefSeq transcript annotation file. If specified, gene-wide and exon-wide events can be detected, and SNVs/LOH events will be annotated with the gene name.", required = false)
    String RefseqFileName = null;

    @Argument(fullName = "indels", required = false, doc = "Indel support")
    boolean DETECT_INDELS = false;

    private MultiBAMHelper multibam_helper = null;

    public static String SEURAT_VERSION = "2.6";
    public static final int QUALITY_CAP = 255;

    private SeekableRODIterator refseqIterator;

    PrintStream gene_out_file;
    Map<String, List<RegionalAnalysisWalker>> analyzers = new HashMap<String, List<RegionalAnalysisWalker>>();
    Class[] analyses = {SomaticMutationAnalysis.class, ExonDeletionAnalysis.class, LOHAnalysis.class/*, StructuralVariationAnalysis.class */, AllelicShiftAnalysis.class};

    List<RegionalAnalysisWalker> nt_walkers = new ArrayList<RegionalAnalysisWalker>();

    @Override
    public void initialize() {
        //prepare multi-BAM support
        multibam_helper = new MultiBAMHelper(getToolkit());
        multibam_helper.setMinPerSampleCoverage(MIN_COVERAGE);
        multibam_helper.setMaxMismatches(arguments.maximum_mismatches);

        try {
            gene_out_file = new PrintStream(gene_out);
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        //annotations
        if (RefseqFileName != null) {
            logger.info("Using RefSeq annotations from " + RefseqFileName);

            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                    getToolkit().getGenomeLocParser(),
                    getToolkit().getArguments().unsafe);
            RMDTrack refseq = builder.createInstanceOfTrack(RefSeqCodec.class, new File(RefseqFileName));

            refseqIterator = new SeekableRODIterator(refseq.getHeader(),
                    refseq.getSequenceDictionary(), getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                    getToolkit().getGenomeLocParser(),
                    refseq.getIterator());

            if (refseqIterator == null) {
                logger.info("No gene annotations available");
                if (arguments.coding_only)
                    throw new UserException.BadInput("A gene annotation file must be provided for the \"coding_only\" option to be allowed.");

            }
        }

        //non-transcriptional (global) walkers
        nt_walkers = new ArrayList<RegionalAnalysisWalker>();

        //initialize walkers for the global context
        for (Class analysis_class : analyses) {

            RegionalAnalysisWalker r = null;
            try {
                r = (RegionalAnalysisWalker) analysis_class.newInstance();
            } catch (InstantiationException e) {
                e.printStackTrace();
            } catch (IllegalAccessException e) {
                e.printStackTrace();
            }
            if (!r.isTranscriptionalAnalysis() && r.initialize(multibam_helper.Evidence, arguments, null)) {
                nt_walkers.add(r);
                if (arguments.enable_debug)
                    System.err.printf("Enabled walker: %s\n", analysis_class.getSimpleName());
            }
        }

        writeVCFHeaders();


        gene_out_file.println("#context\tgene\tevent\tquality\tinfo");
    }

    //output file headers
    private void writeVCFHeaders() {

        out.println("##fileformat=VCFv4.1");
        out.println("##source=Seurat-" + SEURAT_VERSION);
        out.println("##seuratarguments=" + arguments.toString());

        //INFO headers
        out.println(VCFInfo("TYPE", "1", "String", "The type of somatic change detected"));
        out.println(VCFInfo("PILEUP1", "1", "String", "The pileup for the normal"));
        out.println(VCFInfo("PILEUP2", "1", "String", "The pileup for the tumor"));

        out.println(VCFInfo("AR1", "1", "Float", "Allele frequency of ALT allele in normal"));
        out.println(VCFInfo("AR2", "1", "Float", "Allele frequency of ALT allele in tumor"));

        out.println("##INFO=<ID=DP1,Number=1,Type=Integer,Description=\"The depth of coverage in normal\">");
        out.println("##INFO=<ID=DP2,Number=1,Type=Integer,Description=\"The depth of coverage in tumor\">");
        out.println("##INFO=<ID=SEQ,Number=1,Type=String,Description=\"The bases inserted\">");
        out.println("##INFO=<ID=LN,Number=1,Type=Integer,Description=\"The length of a change\">");
        out.println(VCFInfo("MVC1", "1", "String", "The median for the variant evidence distance from the end of the read (normal)"));
        out.println(VCFInfo("MVBQ1", "1", "String", "The median for the base quality of variant evidence in the normal sample"));
        out.println(VCFInfo("MVMQ1", "1", "String", "The median for the mapping quality of variant evidence in the normal sample"));
        out.println(VCFInfo("MVC2", "1", "String", "The median for the variant evidence distance from the end of the read (tumor)"));
        out.println(VCFInfo("MVBQ2", "1", "String", "The median for the base quality of variant evidence in the tumor sample"));
        out.println(VCFInfo("MVMQ2", "1", "String", "The median for the mapping quality of variant evidence in the tumor sample"));

        out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    }

    private static String VCFInfo(String id, String number, String type, String description) {
        return String.format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">", id, number, type, description);
    }

    @Override
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //skip slow BAM separation if minimum coverage for each sample is impossible to be covered

        if (arguments.enable_debug) {
            System.err.printf("%s - %d\n", ref.getLocus().toString(), context.hasExtendedEventPileup() ? context.getExtendedEventPileup().depthOfCoverage() : context.getBasePileup().depthOfCoverage());
        }

        if (context.hasExtendedEventPileup()) {
            return multibam_helper.SetBasePileup(context.getExtendedEventPileup().getBaseAndMappingFilteredPileup(MIN_BASE_QUALITY_SCORE, MIN_MAPPING_QUALITY_SCORE), ref);
        } else {
            return multibam_helper.SetBasePileup(context.getBasePileup().getBaseAndMappingFilteredPileup(MIN_BASE_QUALITY_SCORE, MIN_MAPPING_QUALITY_SCORE), ref);
        }
    }

    /**
     * The map function runs once per single-base locus, and accepts a 'context', a
     * data structure consisting of the reads which overlap the locus, the sites over
     * which they fall, and the base from the reference that overlaps.
     *
     * @param tracker The accessor for reference metadata.
     * @param ref     The reference base that lines up with this locus.
     * @param context Information about reads aligning to this locus.
     * @return In this case, returns a count of how many loci were seen at this site (1).
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        //get refseq annotation
        RODRecordList annotationList = (refseqIterator == null ? null : refseqIterator.seekForward(context.getLocation()));

        List<GeneContext> geneContexts = getAnnotationContext(annotationList);
        List<RegionalAnalysisWalker> active_walkers = new ArrayList<RegionalAnalysisWalker>();

        boolean any_coding = false;

        //regional walkers
        for (GeneContext geneContext : geneContexts) {
            if (geneContext.context == GeneContext.GeneContextClass.Transcript && geneContext.coding == false && arguments.coding_only == true)
                continue;

            if (geneContext.coding)
                any_coding = true;

            List<RegionalAnalysisWalker> walkers = analyzers.get(geneContext.name);
            if (walkers == null) {
                walkers = new ArrayList<RegionalAnalysisWalker>();
                analyzers.put(geneContext.name, walkers);
                //initialize walkers for these contexts
                for (Class analysis_class : analyses) {
                    RegionalAnalysisWalker r = null;
                    try {
                        r = (RegionalAnalysisWalker) analysis_class.newInstance();
                    } catch (InstantiationException e) {
                        e.printStackTrace();
                    } catch (IllegalAccessException e) {
                        e.printStackTrace();
                    }
                    if (r.isTranscriptionalAnalysis() && r.initialize(multibam_helper.Evidence, arguments, geneContext))
                        walkers.add(r);
                }
            }
            active_walkers.addAll(walkers);
        }

        if (!arguments.coding_only || any_coding)
            active_walkers.addAll(nt_walkers);

        for (RegionalAnalysisWalker walker : active_walkers) {
            if (walker.filter(tracker, ref, context))
                walker.map(tracker, ref, context, geneContexts);
        }


        return 0;

    }


    private void WriteEvents(List<Event> events, String context_name) {
        if (events == null)
            return;

        for (Event e : events) {

            double phred = -10 * Math.log10(1 - e.getProbability());

            if (phred > QUALITY_CAP)
                phred = QUALITY_CAP;

            if (phred < arguments.quality && !(OUTPUT_ALL))
                continue;

            if (e.getLocation() == null) {
                String gene_name = (e.getGeneContext().gene_name == null ? "" : e.getGeneContext().gene_name);
                gene_out_file.println(context_name + '\t' + gene_name + '\t' + e.toString());
            } else {
                String str_gene = "";
                if (!context_name.equals("unknown"))
                    str_gene = "CX=" + context_name + ';';

                out.printf("%s\t%d\t.\t%s\t%s\t%.1f\tPASS\tTYPE=%s;%s%s\n",
                        e.getLocation().getContig(), e.getLocation().getStart(),
                        e.getAttribute("REF"), e.getAttribute("ALT"), phred,
                        e.getName(), str_gene, e.attributeStringVCF());

            }

        }
    }


    /**
     * Provides an initial value for the reduce function.  Hello walker counts loci,
     * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
     *
     * @return 0.
     */
    @Override
    public Long reduceInit() {
        return 0L;
    }

    @Override
    public boolean includeReadsWithDeletionAtLoci() {
        return DETECT_INDELS;
    }

    @Override
    public boolean generateExtendedEvents() {
        return DETECT_INDELS;
    }

    /**
     * Combines the result of the latest map with the accumulator.  In inductive terms,
     * this represents the step loci[x + 1] = loci[x] + 1
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return The total count of loci processed so far.
     */
    @Override
    public Long reduce(Integer value, Long sum) {
        return sum + value;
    }

    /**
     * Retrieves the final result of the traversal.
     *
     * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
     *               by the reduce function.
     */
    @Override
    public void onTraversalDone(Long result) {
        //global walkers output
        for (RegionalAnalysisWalker walker : nt_walkers)
            WriteEvents(walker.reduce(), "unknown");


        //context walkers output                                                                                                          \
        for (Map.Entry<String, List<RegionalAnalysisWalker>> walkerList : analyzers.entrySet()) {
            for (RegionalAnalysisWalker walker : walkerList.getValue())
                WriteEvents(walker.reduce(), walkerList.getKey());
        }


        //hack to prevent wrong exit codes later on, even though we're done
        System.exit(0);
    }

    public static List<GeneContext> getAnnotationContext(RODRecordList annotationList) {
        List<GeneContext> contexts = new ArrayList<GeneContext>();

        if (annotationList != null && annotationList.size() != 0) {
            for (GATKFeature feature : annotationList) {
                Transcript transcript = (Transcript) feature.getUnderlyingObject();

                GeneContext gene_cx = new GeneContext(transcript.getGeneName(), transcript.getGeneName() + '_' + transcript.getTranscriptId(), GeneContext.GeneContextClass.Transcript);

                List<GenomeLoc> exons = transcript.getExons();

                boolean in_any_exon = false;
                boolean in_coding = transcript.overlapsCodingP((annotationList.getLocation()));

                if (in_coding) {
                    for (int x = 1; x <= exons.size(); x++) {
                        if (exons.get(x - 1).overlapsP(annotationList.getLocation())) {
                            contexts.add(new GeneContext(transcript.getGeneName(), (transcript.getGeneName() + '_' + transcript.getTranscriptId() + "_exon_" + x), GeneContext.GeneContextClass.Exon));
                            in_any_exon = true;
                            gene_cx.coding = true;
                            break;
                        }

                    }
                    if (!in_any_exon) //intron, but which intron?
                    {
                        for (int x = 1; x <= exons.size(); x++) {
                            if (annotationList.getLocation().isBefore(exons.get(x))) {

                                contexts.add(new GeneContext(transcript.getGeneName(), (transcript.getGeneName() + '_' + transcript.getTranscriptId() + "_intron_" + x), GeneContext.GeneContextClass.Intron));
                                break;
                            }
                        }

                    }
                } else {
                    contexts.add(new GeneContext(transcript.getGeneName(), (transcript.getGeneName() + transcript.getTranscriptId() + "_UTR"), GeneContext.GeneContextClass.Intron));
                }

                contexts.add(gene_cx);
            }
        } else
            contexts.add(new GeneContext("unknown", GeneContext.GeneContextClass.Nongenic));

        return contexts;
    }

}
