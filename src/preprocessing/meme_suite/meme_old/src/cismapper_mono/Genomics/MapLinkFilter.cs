//--------------------------------------------------------------------------------
// <copyright file="MapLinkFilter.cs" 
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
    /// Map link filter.
    /// </summary>
    public class MapLinkFilter
    {
        /// <summary>
        /// The use promoter range.
        /// </summary>
        private bool usePromoterRange;

        /// <summary>
        /// The promoter upstream range.
        /// </summary>
        private int promoterUpstreamRange;

        /// <summary>
        /// The promoter downstream range.
        /// </summary>
        private int promoterDownstreamRange;

        /// <summary>
        /// The confidence threshold.
        /// </summary>
        private double confidenceThreshold;

        /// <summary>
        /// The correlation threshold.
        /// </summary>
        private double correlationThreshold;

        /// <summary>
        /// The absolute correlation threshold.
        /// </summary>
        private double absoluteCorrelationThreshold;

        /// <summary>
        /// The maximum length of the link.
        /// </summary>
        private int maximumLinkLength;

        private TssRegulatoryMap auxiliaryMap;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.MapLinkFilter"/> class.
        /// </summary>
        public MapLinkFilter()
        {
        }

        /// <summary>
        /// Link type.
        /// </summary>
        public enum LinkType
        {
            /// <summary>
            /// Any.
            /// </summary>
            Any,

            /// <summary>
            /// The link must be to the nearest Locus
            /// </summary>
            NearestLocus,

            /// <summary>
            /// The link for each Locus must be to the nearest gene
            /// </summary>
            NearestGene,

            /// <summary>
            /// The link for each Locus must be to the nearest gene 
            /// AND that link must be the highest confidence Locus linked
            /// </summary>
            NearestBestGene,

            /// <summary>
            /// The link for each gene must be to its closest Locus or the link for each Locus must
            /// be to the closest gene
            /// </summary>
            NearestLocusOrGene,

            /// <summary>
            /// The nearest Locus and gene.
            /// </summary>
            NearestLocusAndGene,

            /// <summary>
            /// The Locus that is the closest of the nearest gene links
            /// </summary>
            NearestLocusOfNearestGene,

            /// <summary>
            /// The link must be the best Locus
            /// </summary>
            BestLocusLink,

            /// <summary>
            /// The worst locus link.
            /// </summary>
            WorstLocusLink,

            /// <summary>
            /// The Locus is the best for that gene and that gene is the 
            /// nearest to that Locus
            /// </summary>
            BestLocusNearestGene,

            /// <summary>
            /// The Locus is the best for that gene and that gene is the 
            /// best for that Locus
            /// </summary>
            BestLocusBestGene,

            /// <summary>
            /// The link for each Locus must be to the best gene
            /// </summary>
            BestGeneLink,

            /// <summary>
            /// Excludes the closest Locus to the gene
            /// </summary>
            NotNearestLocus,

            /// <summary>
            /// Excludes the gene closest to each Locus
            /// </summary>
            NotNearestGene,

            /// <summary>
            /// Excludes both the closest Locus to each gene and closest gene to each Locus
            /// </summary>
            NotNearestLocusOrGene,

            /// <summary>
            /// Links between Loci and genes that are not nearest neighbors by 
            /// either heuristic, BUT it is the best Locus linked to that gene
            /// </summary>
            BestLocusNotNearestLocusOrGene,

            /// <summary>
            /// The best distal Locus.
            /// </summary>
            BestDistalLocus,
        }

        /// <summary>
        /// Gets or sets the promoter upstream range.
        /// </summary>
        /// <value>The promoter upstream range.</value>
        public int PromoterUpstreamRange 
        {
            get
            {
                return this.promoterUpstreamRange;
            }

            set
            {
                this.usePromoterRange = true;
                this.promoterUpstreamRange = value;
            }
        }

        /// <summary>
        /// Gets or sets the promoter downstream range.
        /// </summary>
        /// <value>The promoter downstream range.</value>
        public int PromoterDownstreamRange
        {
            get
            {
                return this.promoterDownstreamRange;
            }

            set
            {
                this.usePromoterRange = true;
                this.promoterDownstreamRange = value;
            }
        }

        /// <summary>
        /// Gets or sets the maximum length of the link.
        /// </summary>
        /// <value>The maximum length of the link.</value>
        public int MaximumLinkLength 
        {
            get
            {
                return this.maximumLinkLength;
            }

            set
            {
                this.UseMaximumLinkLength = true;
                this.maximumLinkLength = value;
            }
        }

        /// <summary>
        /// Gets or sets the confidence threshold.
        /// </summary>
        /// <value>The confidence threshold.</value>
        public double ConfidenceThreshold 
        {
            get
            {
                return this.confidenceThreshold;
            }

            set
            {
                this.UseConfidenceThreshold = true;
                this.confidenceThreshold = value;
            }
        }

        /// <summary>
        /// Gets or sets the correlation threshold.
        /// </summary>
        /// <value>The correlation threshold.</value>
        public double CorrelationThreshold 
        {
            get
            {
                return this.correlationThreshold;
            }

            set
            {
                if (this.UseAbsoluteCorrelationThreshold)
                {
                    throw new Exception("Error: cannot use both a correlation and an absolute correlation threshold");
                }

                this.UseCorrelationThreshold = true;
                this.correlationThreshold = value;
            }
        }

        /// <summary>
        /// Gets or sets the absolute correlation threshold.
        /// </summary>
        /// <value>The absolute correlation threshold.</value>
        public double AbsoluteCorrelationThreshold
        {
            get
            {
                return this.absoluteCorrelationThreshold;
            }

            set
            {
                if (this.UseCorrelationThreshold)
                {
                    throw new Exception("Error: cannot use both a correlation and an absolute correlation threshold");
                }

                this.UseAbsoluteCorrelationThreshold = true;
                this.absoluteCorrelationThreshold = value;
            }
        }

        /// <summary>
        /// Gets a value indicating whether this <see cref="Genomics.MapLinkFilter"/> use promoter range.
        /// </summary>
        /// <value><c>true</c> if use promoter range; otherwise, <c>false</c>.</value>
        public bool UsePromoterRange 
        {
            get
            {
                return this.usePromoterRange && this.PromoterUpstreamRange != 0 && this.PromoterDownstreamRange != 0;
            }
        }

        /// <summary>
        /// Gets a value indicating whether this <see cref="Genomics.MapLinkFilter"/> use maximum link length.
        /// </summary>
        /// <value><c>true</c> if use maximum link length; otherwise, <c>false</c>.</value>
        public bool UseMaximumLinkLength { get; private set; }

        /// <summary>
        /// Gets a value indicating whether this <see cref="Genomics.MapLinkFilter"/> use confidence threshold.
        /// </summary>
        /// <value><c>true</c> if use confidence threshold; otherwise, <c>false</c>.</value>
        public bool UseConfidenceThreshold { get; private set; }

        /// <summary>
        /// Gets a value indicating whether this <see cref="Genomics.MapLinkFilter"/> use correlation threshold.
        /// </summary>
        /// <value><c>true</c> if use correlation threshold; otherwise, <c>false</c>.</value>
        public bool UseCorrelationThreshold { get; private set; }

        /// <summary>
        /// Gets a value indicating whether this <see cref="Genomics.MapLinkFilter"/> use absolute correlation threshold.
        /// </summary>
        /// <value><c>true</c> if use absolute correlation threshold; otherwise, <c>false</c>.</value>
        public bool UseAbsoluteCorrelationThreshold { get; private set; }

        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Genomics.MapLinkFilter"/> keep best link only.
        /// </summary>
        /// <value><c>true</c> if keep best link only; otherwise, <c>false</c>.</value>
        public LinkType LinkTypeFilter { get; set; }

        /// <summary>
        /// Gets or sets the locus set.
        /// </summary>
        /// <value>The locus set.</value>
        public HashSet<string> LocusSet { get; set; }

        /// <summary>
        /// Optional expression data set for TSS location inforamation
        /// </summary>
        /// <value>The expression data.</value>
        public IExpressionData ExpressionData { get; set; }

        /// <summary>
        /// Gets or sets the transcript set to restrict the map to.
        /// </summary>
        /// <value>The transcript set.</value>
        public HashSet<string> TranscriptSet { get; set; }

        /// <summary>
        /// Determines whether this instance is valid link the specified strand position confidence correlation.
        /// </summary>
        /// <returns><c>true</c> if this instance is valid link the specified strand position confidence correlation; otherwise, <c>false</c>.</returns>
        /// <param name="strand">Transcript/gene strand.</param>
        /// <param name="position">Position of the Locus relative to the transcript/gene.</param>
        /// <param name="confidence">Confidence of the link.</param>
        /// <param name="correlation">Correlation of the link.</param>
        public bool IsValidLink(string strand, int position, double confidence, double correlation)
        {
            if (strand == "-")
            {
                position = -position;
            }

            if (this.UseMaximumLinkLength && this.MaximumLinkLength >= 0 && Math.Abs(position) > this.MaximumLinkLength)
            {
                return false;
            }

            bool inPromoter = position > 0 ? 
                position <= this.PromoterDownstreamRange :
                -position <= this.PromoterUpstreamRange;

            if (this.UsePromoterRange && inPromoter)
            {
                return false;
            }

            bool isConfident = confidence <= this.ConfidenceThreshold;

            if (this.UseConfidenceThreshold && !isConfident)
            {
                return false;
            }

            bool isCorrelated = correlation >= this.CorrelationThreshold;

            if (this.UseCorrelationThreshold && !isCorrelated)
            {
                return false;
            }

            bool isAbsoluteCorrelated = Math.Abs(correlation) >= this.CorrelationThreshold;

            if (this.UseAbsoluteCorrelationThreshold && !isAbsoluteCorrelated)
            {
                return false;
            }

            return true;
        }

        public void EnsureKey(MapLink link, TssRegulatoryMap tssMap)
        {
            if (!tssMap.ContainsKey(link.TranscriptName))
            {
                tssMap.Add(link.TranscriptName, new Dictionary<MapLink.Locus, MapLink>());
            }

            if (tssMap[link.TranscriptName].ContainsKey(link.LocusName))
            {
                throw new Exception(string.Format("Duplicate Locus-TSS link {0}-{1} in map", link.LocusName, link.TranscriptName));
            }
        }

        public void EnsureKey(MapLink link, LocusRegulatoryMap LocusMap)
        {
            if (!LocusMap.ContainsKey(link.LocusName))
            {
                LocusMap.Add(link.LocusName, new Dictionary<MapLink.Tss, MapLink>());
            }

            if (LocusMap[link.LocusName].ContainsKey(link.TranscriptName))
            {
                throw new Exception(string.Format("Duplicate Locus-TSS link {0}-{1} in map", link.LocusName, link.TranscriptName));
            }
        }

        /// <summary>
        /// Determines whether this instance is valid link type the specified link map.
        /// </summary>
        /// <returns><c>true</c> if this instance is valid link type the specified link map; otherwise, <c>false</c>.</returns>
        /// <param name="link">Link.</param>
        /// <param name="map">Map.</param>
        public void ApplyLinkFilterType(MapLink link, TssRegulatoryMap tssMap, LocusRegulatoryMap LocusMap)
        {
            this.EnsureKey(link, tssMap);

            if ((this.LinkTypeFilter == LinkType.BestLocusNotNearestLocusOrGene || 
                this.LinkTypeFilter == LinkType.BestDistalLocus)
                && this.auxiliaryMap == null)
            {
                this.auxiliaryMap = new TssRegulatoryMap();
            }

            switch (this.LinkTypeFilter)
            {
                case MapLinkFilter.LinkType.Any:
                    tssMap[link.TranscriptName].Add(link.LocusName, link);
                    break;

                case MapLinkFilter.LinkType.NearestLocus:
                    this.AddReplaceLink(link, tssMap, this.IsClosestLocus);
                    break;

                case MapLinkFilter.LinkType.NearestGene:
                case MapLinkFilter.LinkType.NearestLocusOfNearestGene:
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    break;

                case MapLinkFilter.LinkType.NearestBestGene:
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    break;

                case MapLinkFilter.LinkType.BestLocusNearestGene:
                    this.AddReplaceLink(link, tssMap, this.IsBestLocus);
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    break;

                case MapLinkFilter.LinkType.BestLocusBestGene:
                    this.AddReplaceLink(link, tssMap, this.IsBestLocus);
                    this.AddReplaceLink(link, LocusMap, this.IsBestGene);
                    break;

                case MapLinkFilter.LinkType.NearestLocusOrGene:
                    this.AddReplaceLink(link, tssMap, this.IsClosestLocus);
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    break;

                case MapLinkFilter.LinkType.NearestLocusAndGene:
                    this.AddReplaceLink(link, tssMap, this.IsClosestLocus);
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    break;

                case MapLinkFilter.LinkType.BestLocusLink:
                    this.AddReplaceLink(link, tssMap, this.IsBestLocus);
                    break;

                case MapLinkFilter.LinkType.WorstLocusLink:
                    this.AddReplaceLink(link, tssMap, this.IsWorstLocus);
                    break;

                case MapLinkFilter.LinkType.BestGeneLink:
                    this.AddReplaceLink(link, LocusMap, this.IsBestGene);
                    break;

                case MapLinkFilter.LinkType.NotNearestLocus:
                    break;

                case MapLinkFilter.LinkType.NotNearestGene:
                    tssMap[link.TranscriptName].Add(link.LocusName, link);
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    break;

                case MapLinkFilter.LinkType.NotNearestLocusOrGene:
                    tssMap[link.TranscriptName].Add(link.LocusName, link);
                    this.EnsureKey(link, LocusMap);
                    LocusMap[link.LocusName].Add(link.TranscriptName, link);
                    break;

                case MapLinkFilter.LinkType.BestLocusNotNearestLocusOrGene:
                    tssMap[link.TranscriptName].Add(link.LocusName, link);
                    this.EnsureKey(link, LocusMap);
                    LocusMap[link.LocusName].Add(link.TranscriptName, link);
                    this.EnsureKey(link, this.auxiliaryMap);
                    this.AddReplaceLink(link, this.auxiliaryMap, this.IsBestLocus);
                    break;

                case MapLinkFilter.LinkType.BestDistalLocus:
                    tssMap[link.TranscriptName].Add(link.LocusName, link);
                    this.EnsureKey(link, LocusMap);
                    this.AddReplaceLink(link, LocusMap, this.IsClosestGene);
                    this.EnsureKey(link, this.auxiliaryMap);
                    this.AddReplaceLink(link, this.auxiliaryMap, this.IsClosestLocus);
                    break;
            }
        }

        public TssRegulatoryMap PostProcessLinkFilterType(TssRegulatoryMap tssMap, LocusRegulatoryMap LocusMap)
        {
            switch (this.LinkTypeFilter)
            {
                case MapLinkFilter.LinkType.NearestGene:
                    return LocusMap.Invert();

                case MapLinkFilter.LinkType.NearestBestGene:
                    return new TssRegulatoryMap(LocusMap.Links
                        .ToLookup(x => x.TranscriptName, x => x)
                        .Select(x => x.OrderBy(y => y.ConfidenceScore).First()));

                case MapLinkFilter.LinkType.NearestLocusOfNearestGene:
                    return new TssRegulatoryMap(LocusMap.Links
                        .ToLookup(x => x.TranscriptName, x => x)
                        .Select(x => x.OrderBy(y => y.AbsLinkLength).First()));

                case MapLinkFilter.LinkType.BestGeneLink:
                    return LocusMap.Invert();

                case MapLinkFilter.LinkType.NearestLocusOrGene:
                    return new TssRegulatoryMap(tssMap.Links.Concat(LocusMap.Links));

                case MapLinkFilter.LinkType.NearestLocusAndGene:
                    return new TssRegulatoryMap(tssMap.Links.Where(LocusMap.Links.Contains));

                case MapLinkFilter.LinkType.BestLocusNearestGene:
                    return new TssRegulatoryMap(tssMap.Links.Where(LocusMap.Links.Contains));

                case MapLinkFilter.LinkType.BestLocusBestGene:
                    return new TssRegulatoryMap(tssMap.Links.Where(LocusMap.Links.Contains));

                case MapLinkFilter.LinkType.NotNearestGene:
                    return new TssRegulatoryMap(tssMap.Links.Where(link => LocusMap.Links.Contains(link)));

                case MapLinkFilter.LinkType.NotNearestLocusOrGene:
                    var Locus2nnmap = LocusMap.GetNearestNeighborMap(5);
                    return new TssRegulatoryMap(Locus2nnmap.Links.Where(link => !this.IsClosestLocus(link, tssMap) && !this.IsClosestGene(link, LocusMap)));

                case MapLinkFilter.LinkType.BestLocusNotNearestLocusOrGene:
                    return new TssRegulatoryMap(LocusMap.Links
                        .Where(link => !this.IsClosestLocus(link, tssMap) && 
                            !this.IsClosestGene(link, LocusMap) &&
                            this.auxiliaryMap.Links.Contains(link)));

                case MapLinkFilter.LinkType.BestDistalLocus:
                    {
                        var distalLinks = tssMap.Links
                            .Where(link => !this.IsClosestLocus(link, this.auxiliaryMap) &&
                                          !this.IsClosestGene(link, LocusMap));

                        TssRegulatoryMap bestDistalMap = new TssRegulatoryMap();
                        foreach (var link in distalLinks)
                        {
                            this.EnsureKey(link, bestDistalMap);
                            this.AddReplaceLink(link, bestDistalMap, IsBestLocus);
                        }

                        return bestDistalMap;
                    }

                default:
                    return tssMap;
            }
        }

        public void AddReplaceLink(
            MapLink link, 
            TssRegulatoryMap map,
            Func<MapLink, TssRegulatoryMap, bool> test)
        {
            if (test(link, map))
            {
                map[link.TranscriptName].Clear();
                map[link.TranscriptName].Add(link.LocusName, link);
            }
        }

        public void AddReplaceLink(
            MapLink link, 
            LocusRegulatoryMap map,
            Func<MapLink, LocusRegulatoryMap, bool> test)
        {
            this.EnsureKey(link, map);
            if (test(link, map))
            {
                map[link.LocusName].Clear();
                map[link.LocusName].Add(link.TranscriptName, link);
            }
        }

        /// <summary>
        /// Determines whether this given link is closest link the specified link map.
        /// </summary>
        /// <returns><c>true</c> if this instance is closest link the specified link map; otherwise, <c>false</c>.</returns>
        /// <param name="link">Link.</param>
        /// <param name="map">Map.</param>
        public bool IsClosestLocus(MapLink link, TssRegulatoryMap tssMap)
        {
            return !tssMap[link.TranscriptName].Any(x => x.Value.AbsLinkLength < link.AbsLinkLength);
        }

        /// <summary>
        /// Determines whether this given link is most confident link the specified link map.
        /// </summary>
        /// <returns><c>true</c> if this instance is closest link the specified link map; otherwise, <c>false</c>.</returns>
        /// <param name="link">Link.</param>
        /// <param name="map">Map.</param>
        private bool IsBestLocus(MapLink link, TssRegulatoryMap tssMap)
        {
            return !tssMap[link.TranscriptName].Any(x => x.Value.ConfidenceScore < link.ConfidenceScore);
        }

        /// <summary>
        /// Determines whether this instance is worst Locus the specified link tssMap.
        /// </summary>
        /// <returns><c>true</c> if this instance is worst Locus the specified link tssMap; otherwise, <c>false</c>.</returns>
        /// <param name="link">Link.</param>
        /// <param name="tssMap">Tss map.</param>
        private bool IsWorstLocus(MapLink link, TssRegulatoryMap tssMap)
        {
            return !tssMap[link.TranscriptName].Any(x => x.Value.ConfidenceScore > link.ConfidenceScore);
        }

        /// <summary>
        /// Determines whether this given link is closest gene the specified link map.
        /// </summary>
        /// <returns><c>true</c> if this instance is closest gene the specified link map; otherwise, <c>false</c>.</returns>
        /// <param name="link">Link.</param>
        /// <param name="map">Map.</param>
        public bool IsClosestGene(MapLink link, LocusRegulatoryMap LocusMap)
        {
            return !LocusMap[link.LocusName].Any(x => x.Value.AbsLinkLength < link.AbsLinkLength);
        }

        /// <summary>
        /// Determines whether this given link is closest gene the specified link map.
        /// </summary>
        /// <returns><c>true</c> if this instance is closest gene the specified link map; otherwise, <c>false</c>.</returns>
        /// <param name="link">Link.</param>
        /// <param name="map">Map.</param>
        private bool IsBestGene(MapLink link, LocusRegulatoryMap LocusMap)
        {
            return !LocusMap[link.LocusName].Any(x => x.Value.ConfidenceScore < link.ConfidenceScore);
        }
    }
}