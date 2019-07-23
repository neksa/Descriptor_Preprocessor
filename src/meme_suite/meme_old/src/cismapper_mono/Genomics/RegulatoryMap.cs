//--------------------------------------------------------------------------------
// <copyright file="RegulatoryMap.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using Shared;
    using Tools;

    /// <summary>
    /// Map encapsulation
    /// </summary>
    public abstract class RegulatoryMap<TOrigin, TDestination, TImpl, TReverseImpl> : Dictionary<TOrigin, Dictionary<TDestination, MapLink>> 
        where TOrigin : IEqualityComparer<TOrigin>
        where TDestination : IEqualityComparer<TDestination>
    {
        /// <summary>
        /// The transcripts in the map
        /// </summary>
        protected HashSet<MapLink.Tss> transcripts;

        /// <summary>
        /// The loci in the map
        /// </summary>
        protected HashSet<MapLink.Locus> loci;

        /// <summary>
        /// The tss link degree.
        /// </summary>
        protected Dictionary<MapLink.Tss, int> tssDegree;

        private Dictionary<TOrigin, int> originDegree;

        /// <summary>
        /// The links.
        /// </summary>
        private HashSet<MapLink> links;

        /// <summary>
        /// The directional link lengths.
        /// </summary>
        private List<int> locusPositions;

        /// <summary>
        /// The reversed map.
        /// </summary>
        //private TReverseImpl reversedMap;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.RegulatoryMap"/> class.
        /// </summary>
        public RegulatoryMap()
            : base()
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="RegulatoryMap"/> class from a dictionary
        /// </summary>
        /// <param name="source">Source links.</param>
        public RegulatoryMap(Dictionary<TOrigin, Dictionary<TDestination, MapLink>> source)
            : base(source)
        {
        }

        /// <summary>
        /// Gets the link count.
        /// </summary>
        /// <value>The link count.</value>
        public int LinkCount
        {
            get
            {
                return this.Sum(x => x.Value.Count);
            }
        }

        public abstract HashSet<MapLink.Tss> Transcripts { get; }

        public abstract HashSet<MapLink.Locus> Loci { get; }

        /// <summary>
        /// Gets the origin degree.
        /// </summary>
        /// <value>The origin degree.</value>
        public Dictionary<TOrigin, int> OriginDegree
        {
            get 
            {
                return Helpers.CheckInit(
                    ref this.originDegree,
                    () => this.ToDictionary(x => x.Key, x => x.Value.Count));
            }
        }

        /// <summary>
        /// Gets distinct transcript-Locus links
        /// </summary>
        /// <value>The links.</value>
        public HashSet<MapLink> Links
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.links,
                    () => new HashSet<MapLink>(this
                            .Select(x => x.Value.Values)
                            .SelectMany(x => x), 
                        new MapLink.SameTranscriptAndLocus()));
            }
        }

        /// <summary>
        /// Gets the Locus positions.
        /// </summary>
        /// <value>The Locus positions.</value>
        public List<int> LocusPositions
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.locusPositions,
                    () => this.Links.Select(x => x.LocusPosition).ToList());
            }
        }

        /// <summary>
        /// Loads the map.
        /// </summary>
        /// <returns>The map.</returns>
        /// <param name="mapFileName">Map file name.</param>
        /// <param name="filter">Filter of map links.</param>
        public static TssRegulatoryMap LoadMap(
            string mapFileName,
            MapLinkFilter filter)
        {
            using (TextReader tr = new StreamReader(mapFileName))
            {
                string line = null;

                TssRegulatoryMap tssMap = new TssRegulatoryMap();
                LocusRegulatoryMap LocusMap = new LocusRegulatoryMap();

                while ((line = tr.ReadLine()) != null)
                {
                    var fields = line.Split('\t');


                    var transcriptName = fields[3];
                    var confidence     = double.Parse(fields[8]);
                    var correlation    = double.Parse(fields[7]);
                    var distance       = int.Parse(fields[9]);
                    var LocusName        = fields[6];
                    var strand         = fields[5];

                    var link = new MapLink
                    {
                        ConfidenceScore = confidence,
                        Correlation = correlation,
                        LinkLength = distance,
                        TranscriptName = transcriptName,
                        TssName = transcriptName,
                        LocusName = LocusName,
                        Strand = strand,
                        Chromosome = fields[0],
                        TssPosition = int.Parse(fields[1]),
                        HistoneName = fields[10],
                        GeneName = fields[11],
                    };

                    if ((filter.TranscriptSet == null || (filter.TranscriptSet != null && filter.TranscriptSet.Contains(transcriptName))) &&
                        (filter.LocusSet == null || (filter.LocusSet != null && filter.LocusSet.Contains(LocusName))))
                    {
                        AddLink(ref tssMap, ref LocusMap, filter, link);
                    }
                }

                return filter.PostProcessLinkFilterType(tssMap, LocusMap);
            }
        }

        /// <summary>
        /// Converts Loci of each tss into a flat list
        /// </summary>
        /// <returns>The dictionary list.</returns>
        public Dictionary<TOrigin, List<TDestination>> ToDictionaryList()
        {
            return this.ToDictionary(x => x.Key, x => x.Value.Keys.ToList());
        }


        /// <summary>
        /// Creates a new map with only the specified transcripts
        /// </summary>
        /// <returns>The transcripts.</returns>
        /// <param name="transcripts">Transcript set to be in map.</param>
        public TImpl RestrictTranscripts(HashSet<string> transcripts)
        {
            return (TImpl)Activator.CreateInstance(typeof(TImpl), new object[] 
            {
                this.Links.Where(x => transcripts.Contains(x.TranscriptName))
            });
        }

        /// <summary>
        /// Creates a new map with only the specified transcripts
        /// </summary>
        /// <returns>The transcripts.</returns>
        /// <param name="transcripts">Transcript set to be in map.</param>
        public TImpl RestrictTranscripts(HashSet<MapLink.Tss> transcripts)
        {
            return (TImpl)Activator.CreateInstance(typeof(TImpl), new object[] 
            {
                this.Links.Where(x => transcripts.Contains(x.TranscriptName))
            });
        }


        /// <summary>
        /// Gets the best neighbor map for the best n confidence scores
        /// </summary>
        /// <returns>The best neighbor map.</returns>
        /// <param name="n">Number of best neighbors.</param>
        public TImpl GetBestNeighborMap(int n)
        {
            return (TImpl)Activator.CreateInstance(typeof(TImpl), new object[]
            {
                this.ToDictionary(
                    x => x.Key,
                    x => x.Value.Values
                    .OrderBy(y => y.ConfidenceScore).Take(n).ToDictionary(y => y.LocusName, y => y))
            });
        }

        public TImpl GetWorstNeighborMap(int n)
        {
            return (TImpl)Activator.CreateInstance(typeof(TImpl), new object[]
            {
                this.ToDictionary(
                    x => x.Key,
                    x => x.Value.Values
                    .OrderBy(y => -y.ConfidenceScore).Take(n).ToDictionary(y => y.LocusName, y => y))
            });
        }

        public TImpl GetMedianNeighborMap()
        {
            return (TImpl)Activator.CreateInstance(typeof(TImpl), new object[]
            {
                this.ToDictionary(
                    x => x.Key,
                    x => 
                    {
                        var sortedLoci = x.Value.Values
                            .OrderBy(y => y.ConfidenceScore)
                            .ToList();

                        return sortedLoci.Take(sortedLoci.Count / 2 + 1).Reverse().Take(1).ToDictionary(y => y.LocusName, y => y);
                    })
            });
        }

        public TImpl GetBestNthNeighborMap(int n)
        {
            return (TImpl)Activator.CreateInstance(typeof(TImpl), new object[]
            {
                this.ToDictionary(
                    x => x.Key,
                    x => x.Value.Values
                    .OrderBy(y => y.ConfidenceScore).Take(n).Reverse().Take(1).ToDictionary(y => y.LocusName, y => y))
            });
        }
       
        /// <summary>
        /// Gets the link lengths orientated by TSS strand
        /// </summary>
        /// <returns>The link lengths.</returns>
        /// <param name="expressionData">Expression data.</param>
        public List<int> DirectionalLinkLengths(IExpressionData expressionData)
        {
            return this.Links.Select(x => expressionData.Transcripts.ContainsKey(x.TranscriptName) ? 
                (expressionData.Transcripts[x.TranscriptName].Strand == "+" ? 
                    x.LinkLength : -x.LinkLength) : x.LinkLength).ToList();
        }

        /// <summary>
        /// Adds a link to the map if it passes the link filter criteria.
        /// </summary>
        /// <param name="tssMap">Tss map.</param>
        /// <param name="LocusMap">Locus map.</param>
        /// <param name="filter">Filter.</param>
        /// <param name="link">Link.</param>
        private static void AddLink(
            ref TssRegulatoryMap tssMap,
            ref LocusRegulatoryMap LocusMap,
            MapLinkFilter filter,
            MapLink link)
        {
            if (filter.IsValidLink(link.Strand, link.LinkLength, link.ConfidenceScore, link.Correlation))
            {
                filter.ApplyLinkFilterType(link, tssMap, LocusMap);
            }
        }

       
    }
}
