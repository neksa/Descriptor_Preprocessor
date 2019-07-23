//--------------------------------------------------------------------------------
// <copyright file="ChIPSeqPeakBias.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using Genomics;
    using Shared;
    using Tools;

    /// <summary>
    /// ChIP seq peak bias.
    /// </summary>
    public class ChIPSeqPeakBias 
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.ChIPSeqPeakBias"/> class.
        /// </summary>
        public ChIPSeqPeakBias()
        {
        }

        /// <summary>
        /// Gets or sets the ChIP seq file.
        /// </summary>
        /// <value>The ChIP seq file.</value>
        public string ChIPSeqFile { get; set; }

        /// <summary>
        /// Gets or sets the motif match file.
        /// </summary>
        /// <value>The motif match file.</value>
        public string MotifMatchFile { get; set; }

        /// <summary>
        /// Gets or sets the segmentation file.
        /// </summary>
        /// <value>The segmentation file.</value>
        public string SegmentationFile { get; set; }

        /// <summary>
        /// Execute this instance.
        /// </summary>
        public void Execute()
        {
            var peaks = new BedFile(this.ChIPSeqFile, BedFile.Bed6Plus4Layout);
            var motifMatches = new BedFile(this.MotifMatchFile, BedFile.Bed6Layout);
            var segmentation = new BedFile(this.SegmentationFile, BedFile.Bed6Layout);

            HashSet<string> peaksOverlappingMotifs   = this.GetOverlappingRegions(motifMatches, peaks);
            HashSet<string> peaksOverlappingSegment  = this.GetOverlappingRegions(segmentation, peaks);
            HashSet<string> motifsOverlappingSegment = this.GetOverlappingRegions(segmentation, motifMatches);

            int peaksAtMotifs   = peaksOverlappingMotifs.Count;
            int peaksInSegment  = peaksOverlappingSegment.Count(x => peaksOverlappingMotifs.Contains(x));
            int motifsInSegment = motifsOverlappingSegment.Count;

            int motifs = motifMatches.Locations.Count;
            Console.WriteLine("First peak overlapping segment " + peaksOverlappingSegment.First());
            Console.WriteLine("First peak overlapping motif " + peaksOverlappingMotifs.First());

            Console.WriteLine("Peaks\tPeaksOverlappingMotif\tMotifs\tPeaksInSegment\tPeaksInSegmentOverlappingMotif\tMotifsInSegment");
            Console.WriteLine(string.Join("\t", new List<int> { peaks.Locations.Count, peaksAtMotifs, motifs, peaksOverlappingSegment.Count, peaksInSegment, motifsInSegment }));
        }

        /// <summary>
        /// Gets the overlapping regions.
        /// </summary>
        /// <returns>The overlapping regions.</returns>
        /// <param name="regions">Regions to be overlapped</param>
        /// <param name="overlay">Regions to overlay</param>
        private HashSet<string> GetOverlappingRegions(BedFile regions, BedFile overlay)
        {
            return new HashSet<string>(regions.Locations
                .Where(x => overlay.MaxLocationSize.ContainsKey(x.Value.Chromosome))
                .Select(x => TRFScorer.GetOverlaps(
                    x.Value, 
                    overlay.ChromosomeIndexedLocations, 
                    BedFile.IndexSize, 
                    overlay.MaxLocationSize[x.Value.Chromosome]))
                .SelectMany(x => x)
                .Select(x => x.Name));
        }

        /// <summary>
        /// Execute the ChIP seq peak bias analysis.
        /// </summary>
        public class Executor : IAnalysisExecutor<ChIPSeqPeakBias, ChIPSeqPeakBias.Executor.Arguments>
        {
            /// <summary>
            /// The options.
            /// </summary>
            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { Arguments.ChIPSeqFile, "ChIP seq peaks to classify" },
                        { Arguments.MotifMatchFile, "MEME motif matches in BedLocus format" },
                        { Arguments.SegmentationFile, "Genome segmentation from Segway. Must match cell line of ChIP seq." },
                    };
                }
            }

            /// <summary>
            /// Commany line arguments to the analysis
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The TF ChIP seq file measured from a given cell line.
                /// </summary>
                ChIPSeqFile,

                /// <summary>
                /// The FIMO motif match file Locus bed file for the TF's motif.
                /// </summary>
                MotifMatchFile,

                /// <summary>
                /// The segmentation file from Segway for the cell line analyzed.
                /// </summary>
                SegmentationFile,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Calculates a score for the enrichment of TF ChIP-seq peaks that overlap the\n" +
                        "binding motif for that TF in promoter and enhancer genome segments\n" +
                        "given the expected distribution of peaks in these genome segments\n" +
                        "based on the distribution of MEME motif matches across the segments";
                }
            }

            /// <summary>
            /// Execute the analysis
            /// </summary>
            /// <param name = "commandArgs">Command line arguments</param>
            public override void Execute(Args commandArgs)
            {
                var analysis = new ChIPSeqPeakBias();

                analysis.ChIPSeqFile      = commandArgs.StringEnumArgs[Arguments.ChIPSeqFile];
                analysis.MotifMatchFile   = commandArgs.StringEnumArgs[Arguments.MotifMatchFile];
                analysis.SegmentationFile = commandArgs.StringEnumArgs[Arguments.SegmentationFile];

                analysis.Execute();
            }
        }
    }
}