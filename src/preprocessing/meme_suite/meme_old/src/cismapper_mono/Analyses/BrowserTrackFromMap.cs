//--------------------------------------------------------------------------------
// <copyright file="BrowserTrackFromMap.cs" 
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
    using Data;

    /// <summary>
    /// Creates a pretty genome browser track from a map
    /// </summary>
    public class BrowserTrackFromMap
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.BrowserTrackFromMap"/> class.
        /// </summary>
        /// <param name="xmlFile">Xml file.</param>
        public BrowserTrackFromMap()
        {
        }

        /// <summary>
        /// Gets or sets the map file.
        /// </summary>
        /// <value>The map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the output file.
        /// </summary>
        /// <value>The name of the output file.</value>
        public string OutputFile { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Analyses.BrowserTrackFromMap"/> use genes.
        /// </summary>
        /// <value><c>true</c> if use genets; otherwise, <c>false</c>.</value>
        public bool UseGenes { get; set; }
       
        /// <summary>
        /// Gets or sets the threshold.
        /// </summary>
        /// <value>The threshold.</value>
        public double Threshold { get; set; }

        /// <summary>
        /// Execute this instance.
        /// </summary>
        public void Execute()
        {
            var filter = new MapLinkFilter { };
            if (this.Threshold >= 0)
            {
                filter.ConfidenceThreshold = this.Threshold;
            }

            var corrMap = TssRegulatoryMap.LoadMap(this.MapFileName, filter );
            if (corrMap.Count == 0) return;
            if (this.UseGenes)
            {
                corrMap = corrMap.ConvertToGenes();
            }

            var scores = corrMap.Links.OrderBy(x => x.ConfidenceScore).ToArray();

            const int thresholdCount = 10; 
            double binSize = (double)scores.Length / thresholdCount;
            var binEdges = new double[] { 0.0 }
                .Concat(Enumerable.Range(1, 9).Select(x => scores[(int)Math.Floor(x * binSize)].ConfidenceScore))
                .Concat(new double[] { 1 })
                .Reverse()
                .ToArray();

            var binColors = new string[] { "192,192,192" }
                .Concat(binEdges
                    .Select((x, i) => string.Format(
                                "{0},0,{1}", 
                                (int)((double)i / (thresholdCount - 1) * 255),
                                255 - (int)((double)i / (thresholdCount - 1) * 255))))
                .ToList();

            var links = corrMap.Links.Select(x => 
            {
                var upstreamLocation   = x.LocusLocation.Start < x.TssLocation.DirectionalStart ? x.LocusLocation : x.TssLocation;
                var downstreamLocation = x.LocusLocation.End > x.TssLocation.DirectionalStart ? x.LocusLocation : x.TssLocation;

                int upstreamLength   = x.LocusLocation.Start < x.TssLocation.DirectionalStart ? x.LocusSize : 1;
                int downstreamLength = x.LocusLocation.End > x.TssLocation.DirectionalStart ? x.LocusSize : 1;

                int regionStart = upstreamLocation.Start - (upstreamLength == 1 ? 1 : 0);
                int regionEnd = downstreamLocation.End + (downstreamLength == 1 ? 1 : 0);

                var blockCount = "2";
                var blockSizes = upstreamLength + "," + downstreamLength;
                var blockStarts = x.LinkLength > 0 ? 
                    "0," + (regionEnd - regionStart - 1) :
                    "0," + (regionEnd - regionStart - x.LocusSize);

                if (x.LocusLocation.OverlapsDirectionalStart(x.TssLocation))
                {
                    blockCount = "1";
                    blockStarts = "0";
                    blockSizes = x.LocusSize.ToString();
                }

                return new string[] 
                {
                    x.LocusLocation.Chromosome,
                    regionStart.ToString(),
                    regionEnd.ToString(),
                    x.LocusLocation.Name + "_" + x.TssLocation.Name,
                    x.ConfidenceScore.ToString(),
                    upstreamLocation == x.LocusLocation ? "+" : "-",
                    regionStart.ToString(),
                    regionEnd.ToString(),
                    binColors[Bin(x.ConfidenceScore, binEdges)],
                    blockCount,
                    blockSizes,
                    blockStarts
                };
            })
                .ToList();


            Tables.ToNamedTsvFile(
                this.OutputFile,
                links,
                new string[] {"track name=regulatory_map description=\"Regulatory Map\" itemRgb=\"On\"" });
        }

        /// <summary>
        /// Bin the specified pvalue.
        /// </summary>
        /// <param name="pvalue">Pvalue to bin.</param>
        private int Bin(double pvalue, double[] binEdges)
        {
            if (pvalue == 1) return 0;
            for (int i = 1; i < binEdges.Length; i++)
            {
                if (pvalue < binEdges[i - 1] && pvalue >= binEdges[i])
                {
                    return i - 1;
                }
            }

            throw new Exception("Unexpected pvalue for binning: " + pvalue);
        }

        /// <summary>
        /// Executes the browser track maker.
        /// </summary>
        public class Executor : IAnalysisExecutor<BrowserTrackFromMap, BrowserTrackFromMap.Executor.Arguments>
        {
            /// <summary>
            /// Command line arguments for the analysis
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The map file.
                /// </summary>
                MapFileName,

                /// <summary>
                /// The name of the output file.
                /// </summary>
                OutputFile,

                /// <summary>
                /// The use genes.
                /// </summary>
                UseGenes,

                /// <summary>
                /// The threshold.
                /// </summary>
                Threshold,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Produces pretty tracks of maps for the genome browser";
                }
            }

            /// <summary>
            /// The options.
            /// </summary>
            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { Arguments.MapFileName, "The map file to convert" },
                        { Arguments.OutputFile, "The browser track output" },
                        { Arguments.UseGenes, "Optional flag to force track to use gene targets" },
                        { Arguments.Threshold, "Minimum link confidence to track" },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis
            /// </summary>
            /// <param name = "commandArgs">Command line arguments</param>
            public override void Execute(Args commandArgs)
            {
                var analysis = new BrowserTrackFromMap();

                this.ReflectStringArgs(analysis, new Arguments[] { Arguments.MapFileName, Arguments.OutputFile });

                analysis.UseGenes = commandArgs.Flags.Contains(Arguments.UseGenes.ToString());
                analysis.Threshold = commandArgs.StringEnumArgs.ContainsKey(Arguments.Threshold) ?
                    double.Parse(commandArgs.StringEnumArgs[Arguments.Threshold]) :
                    -1;
                
                analysis.Execute();
            }
        }
    }
}
