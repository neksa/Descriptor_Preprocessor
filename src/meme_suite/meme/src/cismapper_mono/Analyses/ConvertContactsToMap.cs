//--------------------------------------------------------------------------------
// <copyright file="ConvertContactsToMap.cs" 
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

    public class ConvertContactsToMap
    {
        public ConvertContactsToMap()
        {
        }

        /// <summary>
        /// Gets or sets the name of the contact file.
        /// </summary>
        /// <value>The name of the contact file.</value>
        public string ContactFileName { get; set; }

        /// <summary>
        /// Gets or sets the annotation file name.
        /// </summary>
        /// <value>The annotation file.</value>
        public string AnnotationFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the Locus file.
        /// </summary>
        /// <value>The name of the Locus file.</value>
        public string LocusFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the map file.
        /// </summary>
        /// <value>The name of the map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the max range of links to include.
        /// </summary>
        /// <value>The max range.</value>
        public int MaxRange { get; set; }


        public void Convert()
        {
            var contacts = Helpers.GetFileDataLines(this.ContactFileName, true).Select(line =>
            {
                var fields = line.Split('\t');
                var pair = fields[6];
                var ids = pair.Split('_');

                return new
                {
                    Pair = pair,
                    Id1 = ids[0],
                    Chr1 = fields[0],
                    Start1 = int.Parse(fields[1]),
                    End1 = int.Parse(fields[2]),
                    Id2 = ids[1],
                    Chr2 = fields[3],
                    Start2 = int.Parse(fields[4]),
                    End2 = int.Parse(fields[5]),
                    Score = double.Parse(fields[7]),
                };
            }).ToList();

            var locations = contacts.Select(x => new Location[]
            {
                new Location
                {
                    Name = x.Id1, Chromosome = x.Chr1, Start = x.Start1, End = x.End1, AlternateName = x.Id2, Score = x.Score 
                },
                new Location
                {
                    Name = x.Id2, Chromosome = x.Chr2, Start = x.Start2, End = x.End2, AlternateName = x.Id1, Score = x.Score 
                },
            }).SelectMany(x => x)
                .ToList();

            var uniqueLocations = locations
                .ToLookup(x => x.Name)
                .ToDictionary(x => x.Key, x => x.First());

            var contactMap = locations
                .ToLookup(x => x.Name, x => x)
                .ToDictionary(x => x.Key, x => x
                    .ToDictionary(y => y.AlternateName, y => uniqueLocations[y.AlternateName]));

            var indexedLocations = new BedFile(uniqueLocations);

            var LocusFile = new BedFile(this.LocusFileName, BedFile.Bed3Layout);
            var annotationSet = new GtfExpressionFile(GtfExpressionFile.ExpressionType.Cage, this.AnnotationFileName);

            var LocusOverlaps = indexedLocations.Locations.ToDictionary(
                x => x.Key, 
                x => TRFScorer.GetOverlaps(
                    x.Value, 
                    LocusFile.ChromosomeIndexedLocations, 
                    BedFile.IndexSize, 
                    LocusFile.MaxLocationSize[x.Value.Chromosome]));

            var tssOverlaps = indexedLocations.Locations.ToDictionary(
                x => x.Key, 
                x => TRFScorer.GetOverlaps(
                    x.Value, 
                    annotationSet.ChromosomeIndexedLocations, 
                    BedFile.IndexSize, 
                    annotationSet.MaxLocationSize[x.Value.Chromosome],
                    Location.OverlapsDirectionalStart));

            var links = contactMap
                .Where(x => tssOverlaps[x.Key] != null)
                .ToDictionary(x => x.Key, x => x.Value.Where(y => LocusOverlaps[y.Key] != null))
                .Where(x => x.Value.Any())
                .Select(x => tssOverlaps[x.Key].Select(y => new {
                    Id1 = x.Key,
                    Tss = y.Name, 
                    LocusList = x.Value
                        .Select(z => LocusOverlaps[z.Key]
                            .Select(a => new { Name = a.Name, Id2 = z.Key }))
                        .SelectMany(z => z).ToList()
                })
                .Select(y => y.LocusList.Select(z =>
                {
                    var LocusLocation = LocusFile.Locations[z.Name];
                    var tssStart = annotationSet.Transcripts[y.Tss].DirectionalStart;
                    var linkLength = LocusLocation.End < tssStart ? LocusLocation.End - tssStart : LocusLocation.Start - tssStart;
                    return new MapLink 
                    { 
                        TranscriptName = y.Tss,
                        TssName = y.Tss, 
                        LocusName = z.Name,
                        ConfidenceScore = contactMap[y.Id1][z.Id2].Score,
                        LinkLength = linkLength,
                    };
                }))
                .SelectMany(y => y))
                .SelectMany(x => x)
                .ToList();

            var convertedMap = new TssRegulatoryMap(links);

            NullMapBuilder.WriteMap(convertedMap, annotationSet.Locations, "Contact", this.MapFileName);
        }

        /// <summary>
        /// Factory class to create, configure, and execute the experiment
        /// </summary>
        public class Executor : IAnalysisExecutor<ConvertContactsToMap, ConvertContactsToMap.Executor.Arguments>
        {
            /// <summary>
            /// Command line arguments used
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The name of the contact file.
                /// </summary>
                ContactFileName,

                /// <summary>
                /// The Locus file name.
                /// </summary>
                LocusFileName,

                /// <summary>
                /// Output map file name
                /// </summary>
                MapFileName,

                /// <summary>
                /// The gff annotation file name.
                /// </summary>
                AnnotationFileName,

                /// <summary>
                /// The max range.
                /// </summary>
                MaxRange,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Converts a contact map to a map between gene TSSs and Loci";
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
                        { Arguments.ContactFileName,    "Contact map to be converted" },
                        { Arguments.AnnotationFileName, "Gene annotation file used to generate the map" },
                        { Arguments.LocusFileName,        "Loci to be used in the map." },
                        { Arguments.MaxRange,           "Maximum range of links to include" },
                        { Arguments.MapFileName,        "Name of the map file output" },
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

                var converter = new ConvertContactsToMap();

                converter.ContactFileName    = commandArgs.StringEnumArgs[Arguments.ContactFileName];
                converter.AnnotationFileName = commandArgs.StringEnumArgs[Arguments.AnnotationFileName];
                converter.LocusFileName        = commandArgs.StringEnumArgs[Arguments.LocusFileName];
                converter.MaxRange           = commandArgs.IntArgs[Arguments.MaxRange.ToString()];
                converter.MapFileName        = commandArgs.StringEnumArgs[Arguments.MapFileName];

                converter.Convert();
            }
        }
    }
}

