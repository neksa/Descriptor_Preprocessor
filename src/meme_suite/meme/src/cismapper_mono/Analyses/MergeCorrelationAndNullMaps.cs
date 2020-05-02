//--------------------------------------------------------------------------------
// <copyright file="MergeCorrelationAndNullMaps.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2015. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using Analyses;
    using Data;
    using Genomics;
    using Shared;
    using Tools;

    public class MergeCorrelationAndNullMaps : BaseAnalysis
    {
        /// <summary>
        /// The map.
        /// </summary>
        private TssRegulatoryMap map;

        /// <summary>
        /// The null map.
        /// </summary>
        private TssRegulatoryMap nullMap;

        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.MergeCorrelationAndNullMaps"/> class.
        /// </summary>
        public MergeCorrelationAndNullMaps()
            : base(".", ".", ".")
        {
        }

        /// <summary>
        /// Gets or sets the name of the map file.
        /// </summary>
        /// <value>The name of the map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the null map file.
        /// </summary>
        /// <value>The name of the null map file.</value>
        public string NullMapFileName { get; set; }

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
        /// Gets the null map.
        /// </summary>
        /// <value>The null map.</value>
        public TssRegulatoryMap NullMap
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.nullMap,
                    () => TssRegulatoryMap.LoadMap(
                        this.NullMapFileName, 
                        new MapLinkFilter
                    {
                        LinkTypeFilter = MapLinkFilter.LinkType.Any,
                    }));
            }
        }

        public void Execute()
        {
            foreach (var link in this.NullMap.Select(tss => tss.Value.Values.Select(x => new MapLink
            {
                TranscriptName = x.TranscriptName,
                LocusName = x.LocusName,
                GeneName = x.GeneName,
                LinkLength = x.LinkLength,
                TssName = x.TssName,
                Strand = x.Strand,
                TssPosition = x.TssPosition,
                Chromosome = x.Chromosome,
                HistoneName = this.Map.Links.Contains(x) ? 
                    this.Map[x.TranscriptName][x.LocusName].HistoneName : 
                    "None",
                Correlation = this.Map.Links.Contains(x) ? 
                    this.Map[x.TranscriptName][x.LocusName].Correlation : 
                    double.NaN,
                ConfidenceScore = this.Map.Links.Contains(x) ? 
                    this.Map[x.TranscriptName][x.LocusName].ConfidenceScore : 
                    2,
            })).SelectMany(x => x))
            {
                Console.WriteLine(
                    "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", 
                    link.Chromosome, 
                    link.TssPosition,
                    link.TssPosition,
                    link.TranscriptName, 
                    "NA", ////scoreData.CorrectedPvalue, 
                    link.Strand, 
                    link.LocusName, 
                    link.Correlation, 
                    link.ConfidenceScore, 
                    link.LinkLength, 
                    link.HistoneName, 
                    link.GeneName);
            }
        }

        /// <summary>
        /// Executor.
        /// </summary>
        public class Executor : IAnalysisExecutor<MergeCorrelationAndNullMaps, MergeCorrelationAndNullMaps.Executor.Arguments>
        {
            /// <summary>
            /// Command line arguments for the analysis
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The map files to use
                /// </summary>
                MapFileName,

                /// <summary>
                /// The name of the null map file.
                /// </summary>
                NullMapFileName,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Merges the links from a null and a correlation map to produce a single map.\n" +
                        "Items unscored in the correlation map are given a score of 2";
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
                        { Arguments.MapFileName, "Correlation map file" },
                        { Arguments.NullMapFileName, "Null map file" },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis
            /// </summary>
            /// <param name = "commandArgs">Command line arguments</param>
            public override void Execute(Args commandArgs)
            {
                var analysis = new MergeCorrelationAndNullMaps();

                this.CommandArgs = commandArgs;

                this.ReflectStringArgs(analysis, new Arguments[]
                {
                    Arguments.MapFileName,
                    Arguments.NullMapFileName,
                });

                analysis.Execute();
            }
        }
    }
}

