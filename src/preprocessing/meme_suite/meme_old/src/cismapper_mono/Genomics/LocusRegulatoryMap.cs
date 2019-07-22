
namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Shared;

    public class LocusRegulatoryMap : RegulatoryMap<MapLink.Locus, MapLink.Tss, LocusRegulatoryMap, TssRegulatoryMap>
    {
        public LocusRegulatoryMap()
        {
        }

        public LocusRegulatoryMap(Dictionary<MapLink.Locus, Dictionary<MapLink.Tss, MapLink>> source)
            : base (source)
        {
        }


        /// <summary>
        /// Initializes a new instance of the <see cref="RegulatoryMap"/> class.
        /// </summary>
        /// <param name="source">Source links.</param>
        public LocusRegulatoryMap(IEnumerable<MapLink> source)
            : base(source.ToLookup(x => x.LocusName, x => x)
                .ToDictionary(x => x.Key, x => x.ToLookup(y => y.TranscriptName, y => y).ToDictionary(y => y.Key, y => y.First())))
        {
        }


        /// <summary>
        /// Gets the transcripts in the map
        /// </summary>
        /// <value>The transcripts.</value>
        public override HashSet<MapLink.Tss> Transcripts
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcripts, 
                    () => new HashSet<MapLink.Tss>(this.Values.Select(x => x.Keys).SelectMany(x => x)));
            }
        }


        /// <summary>
        /// Gets the Loci in the map
        /// </summary>
        /// <value>The Loci.</value>
        public override HashSet<MapLink.Locus> Loci
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.loci,
                    () => new HashSet<MapLink.Locus>(this.Keys));
            }
        }

        /// <summary>
        /// Gets the nearest neighbor map for the nearest N neighbors
        /// </summary>
        /// <returns>The nearest neighbor map.</returns>
        /// <param name="n">Number of nearest neighbors.</param>
        public LocusRegulatoryMap GetNearestNeighborMap(int n)
        {
            return new LocusRegulatoryMap(this.Select(
                x => x.Value.Values.OrderBy(y => y.AbsLinkLength).Take(n))
                .SelectMany(x => x));
        }

        /// <summary>
        /// Reverse this instance.
        /// </summary>
        /// <returns>A regulatory map with Loci as keys</returns>
        public TssRegulatoryMap Invert()
        {
            return new TssRegulatoryMap(this.Links
                .ToLookup(x => x.TranscriptName, x => x)
                .ToDictionary(x => x.Key, x => x.ToDictionary(y => y.LocusName, y => y)));
        }
    }
}

