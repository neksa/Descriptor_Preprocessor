//--------------------------------------------------------------------------------
// <copyright file="ConvertMapTOScoredRegionPairs.cs" 
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

    public class ConvertMapToScoredRegionPairs
    {
        /// <summary>
        /// The map.
        /// </summary>
        private TssRegulatoryMap map;

        /// <summary>
        /// The expression.
        /// </summary>
        private IExpressionData expression;

        /// <summary>
        /// The Locus file.
        /// </summary>
        private BedFile locusFile;

        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.ConvertMapToScoredRegionPairs"/> class.
        /// </summary>
        public ConvertMapToScoredRegionPairs()
        {
        }

        /// <summary>
        /// Gets or sets the name of the annotation file.
        /// </summary>
        /// <value>The name of the annotation file.</value>
        public string AnnotationFileName { get; set; }

        /// <summary>
        /// Gets or sets the rna source.
        /// </summary>
        /// <value>The rna source.</value>
        public string RnaSource { get; set; }

        /// <summary>
        /// Gets or sets the name of the map file.
        /// </summary>
        /// <value>The name of the map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the Locus file.
        /// </summary>
        /// <value>The name of the Locus file.</value>
        public string LocusFileName { get; set; }

        /// <summary>
        /// Gets the expression.
        /// </summary>
        /// <value>The expression.</value>
        public IExpressionData Expression
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.expression,
                    () => IUnknown.QueryInterface<IExpressionData>(new GtfExpressionFile(
                        GtfExpressionFile.ExpressionTypeFromString(this.RnaSource), 
                        this.AnnotationFileName)));
            }
        }

        /// <summary>
        /// Gets the map.
        /// </summary>
        /// <value>The map.</value>
        public TssRegulatoryMap Map
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.map,
                    () => TssRegulatoryMap.LoadMap(
                        this.MapFileName, 
                        new MapLinkFilter
                        {
                            LinkTypeFilter = MapLinkFilter.LinkType.Any,
                        }));
            }
        }

        /// <summary>
        /// Gets the Locus file.
        /// </summary>
        /// <value>The Locus file.</value>
        public BedFile LocusFile
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.locusFile, 
                    () => new BedFile(this.LocusFileName, BedFile.Bed3Layout));
            }   
        }

        /// <summary>
        /// Convert this instance.
        /// </summary>
        public void Convert()
        {
            Console.WriteLine(string.Join("\n", this.Map.Links.Select(x => new 
                {
                    TssLocation = this.Expression.Transcripts[x.TranscriptName],
                    LocusLocation = this.LocusFile.Locations[x.LocusName],
                    Score = x.ConfidenceScore,
                })
                .Select(x => string.Join("\t", new string[]
                {
                    x.TssLocation.Chromosome,
                    x.TssLocation.DirectionalStart.ToString(),
                    x.TssLocation.DirectionalStart.ToString(),
                    x.LocusLocation.Chromosome.ToString(),
                    x.LocusLocation.Start.ToString(),
                    x.LocusLocation.End.ToString(),
                    x.Score.ToString()
                }))));
        }

        /// <summary>
        /// Factory class to create, configure, and execute the experiment
        /// </summary>
        public class Executor : IAnalysisExecutor<ConvertMapToScoredRegionPairs, ConvertMapToScoredRegionPairs.Executor.Arguments>
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
                /// The name of the Locus file.
                /// </summary>
                LocusFileName,

                /// <summary>
                /// The rna source.
                /// </summary>
                RnaSource,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Converts a map to a region pair and score tsv file";
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
                        { Arguments.LocusFileName,          "Locus regions mapped." },
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

                var converter = new ConvertMapToScoredRegionPairs();

                converter.MapFileName        = commandArgs.StringEnumArgs[Arguments.MapFileName];
                converter.AnnotationFileName = commandArgs.StringEnumArgs[Arguments.AnnotationFileName];
                converter.RnaSource          = commandArgs.StringEnumArgs[Arguments.RnaSource];
                converter.LocusFileName         = commandArgs.StringEnumArgs[Arguments.LocusFileName];

                converter.Convert();
            }
        }
    }
}

