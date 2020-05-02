//--------------------------------------------------------------------------------
// <copyright file="ConvertMapToGenes.cs" 
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
    using System.Text.RegularExpressions;
    using Genomics;
    using Shared;
    using Tools;

    /// <summary>
    /// Converts a map from transcript targets to gene targets.
    /// </summary>
    public class ConvertMapToGenes
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.ConvertMapToGenes"/> class.
        /// </summary>
        public ConvertMapToGenes()
        {
        }

        /// <summary>
        /// Gets or sets the map file name.
        /// </summary>
        /// <value>The map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the annotation file name.
        /// </summary>
        /// <value>The annotation file.</value>
        public string AnnotationFileName { get; set; }

        /// <summary>
        /// Gets or sets the rna source used to build the map.
        /// </summary>
        /// <value>The rna source.</value>
        public string RnaSource { get; set; }

        /// <summary>
        /// Gets or sets the name of the histone used to build the map.
        /// </summary>
        /// <value>The name of the histone.</value>
        public string HistoneName { get; set; }

        /// <summary>
        /// Gets or sets the max range of links to include.
        /// </summary>
        /// <value>The max range.</value>
        public int MaxRange { get; set; }

        /// <summary>
        /// Gets or sets the pvalue maximum threshold to include.
        /// </summary>
        /// <value>The pvalue threshold.</value>
        public double PvalueThreshold { get; set; }

        /// <summary>
        /// Convert the map to genes.
        /// </summary>
        public void Convert()
        {
            MapLinkFilter filter = new MapLinkFilter
            {
                MaximumLinkLength = this.MaxRange,
                ConfidenceThreshold = this.PvalueThreshold,
                LinkTypeFilter = MapLinkFilter.LinkType.Any,
            };

            TssRegulatoryMap map = TssRegulatoryMap.LoadMap(this.MapFileName, filter);

            var expression = new GtfExpressionFile(
                GtfExpressionFile.ExpressionTypeFromString(this.RnaSource), 
                this.AnnotationFileName);

            var geneMap = map.ConvertToGenes(IUnknown.QueryInterface<IExpressionData>(expression));

            foreach (var link in geneMap.Links)
            {
                var geneLocation = expression.Transcripts[link.TranscriptName];

                string[] lineData = new string[]
                {
                    geneLocation.Chromosome,
                    geneLocation.Start.ToString(), 
                    geneLocation.End.ToString(),
                    link.GeneName,
                    "NA",
                    geneLocation.Strand,
                    link.LocusName,
                    link.Correlation.ToString(),
                    link.ConfidenceScore.ToString(),
                    link.LinkLength.ToString(),
                    this.HistoneName,
                    link.TranscriptName,
                };

                Console.WriteLine(string.Join("\t", lineData));
            }
        }

        /// <summary>
        /// Factory class to create, configure, and execute the experiment
        /// </summary>
        public class Executor : IAnalysisExecutor<ConvertMapToGenes, ConvertMapToGenes.Executor.Arguments>
        {
            /// <summary>
            /// Command line arguments used
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The map file name.
                /// </summary>
                MapFileName,

                /// <summary>
                /// The gff annotation file name.
                /// </summary>
                AnnotationFileName,

                /// <summary>
                /// The rna source used to build the map.
                /// </summary>
                RnaSource,

                /// <summary>
                /// The name of the histone used to build the map.
                /// </summary>
                HistoneName,

                /// <summary>
                /// The max range.
                /// </summary>
                MaxRange,

                /// <summary>
                /// The pvalue threshold.
                /// </summary>
                PvalueThreshold,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Converts a map from transcript to gene targets taking the link with the minimum p-value";
                }
            }

            /// <summary>
            /// Gets the command line option descriptions.
            /// </summary>
            /// <value>The options data.</value>
            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { Arguments.MapFileName,        "Map to be converted" },
                        { Arguments.AnnotationFileName, "Gene annotation file used to generate the map" },
                        { Arguments.RnaSource,          "The RNA source used to build the map." },
                        { Arguments.HistoneName,        "The name of the histone used to build the map." },
                        { Arguments.MaxRange,           "Maximum range of links to include" },
                        { Arguments.PvalueThreshold,    "Maximum p-value threshold to include" },
                    };
                }
            }

            /// <summary>
            /// Execute the experiment
            /// </summary>
            /// <param name = "commandArgs">Command line arguments</param>
            public override void Execute(Args commandArgs)
            {
                commandArgs.Validate(typeof(Arguments));

                var converter = new ConvertMapToGenes();

                converter.MapFileName        = commandArgs.StringEnumArgs[Arguments.MapFileName];
                converter.AnnotationFileName = commandArgs.StringEnumArgs[Arguments.AnnotationFileName];
                converter.RnaSource          = commandArgs.StringEnumArgs[Arguments.RnaSource];
                converter.HistoneName        = commandArgs.StringEnumArgs[Arguments.HistoneName];
                converter.MaxRange           = commandArgs.IntArgs[Arguments.MaxRange.ToString()];
                converter.PvalueThreshold    = double.Parse(commandArgs.StringEnumArgs[Arguments.PvalueThreshold]);

                converter.Convert();
            }
        }
    }
}