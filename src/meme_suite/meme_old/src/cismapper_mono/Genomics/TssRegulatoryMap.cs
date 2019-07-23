
namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Shared;
    using Tools;

    public class TssRegulatoryMap : RegulatoryMap<MapLink.Tss, MapLink.Locus, TssRegulatoryMap, LocusRegulatoryMap>
    {

        public TssRegulatoryMap()
        {
        }

        public TssRegulatoryMap(Dictionary<MapLink.Tss, Dictionary<MapLink.Locus, MapLink>> source)
            : base (source)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="RegulatoryMap"/> class.
        /// </summary>
        /// <param name="source">Source links.</param>
        public TssRegulatoryMap(IEnumerable<MapLink> source)
            : base(source.ToLookup(x => x.TranscriptName, x => x)
                .ToDictionary(x => x.Key, x => x
                    .ToLookup(y => y.LocusName, y => y)
                    .ToDictionary(y => y.Key, y => y.First())))
        {
        }


        /// <summary>
        /// Convert the specified source dictionary to a map
        /// </summary>
        /// <returns>A regulatory map</returns>
        /// <param name="source">Source links.</param>
        public static TssRegulatoryMap Convert(Dictionary<MapLink.Tss, Dictionary<MapLink.Locus, MapLink>> source)
        {
            return new TssRegulatoryMap(source);
        }


        /// <summary>
        /// Convert the specified source key-value pairs to a map
        /// </summary>
        /// <returns>A regulatory map</returns>
        /// <param name="source">Source links.</param>
        public static TssRegulatoryMap Convert(IEnumerable<KeyValuePair<MapLink.Tss, Dictionary<MapLink.Locus, MapLink>>> source)
        {
            return new TssRegulatoryMap(source.ToDictionary(x => x.Key, x => x.Value));
        }

        /// <summary>
        /// Converts a flat list of links into a regulatory map.
        /// </summary>
        /// <param name="source">Source links.</param>
        public static implicit operator TssRegulatoryMap(List<MapLink> source)
        {
            return new TssRegulatoryMap(source);
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
                    () => new HashSet<MapLink.Tss>(this.Keys));
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
                    () => new HashSet<MapLink.Locus>(this.Values.Select(x => x.Keys).SelectMany(x => x)));
            }
        }

        /// <summary>
        /// Gets the tss degree.
        /// </summary>
        /// <value>The tss degree.</value>
        public Dictionary<MapLink.Tss, int> TssDegree
        {
            get 
            {
                return Helpers.CheckInit(
                    ref this.tssDegree,
                    () => this.ToDictionary(x => x.Key, x => x.Value.Count));
            }
        }

        /// <summary>
        /// Converts to genes.
        /// </summary>
        /// <returns>The to genes.</returns>
        /// <param name="expression">Expression data.</param>
        public TssRegulatoryMap ConvertToGenes(IExpressionData expression)
        {
            return Convert(expression.Transcripts
                .Where(x => this.ContainsKey((MapLink.Tss)x.Key))
                .ToLookup(
                    x => new MapLink.Tss(x.Value.AlternateName),
                    x => this[x.Value.Name].Values)
                .ToDictionary(
                    x => x.Key, 
                    x => new Dictionary<MapLink.Locus, MapLink>(
                        x.SelectMany(y => y)    //// x.Key is gene name
                        .ToLookup(y => y.LocusName, y => y)
                        .ToDictionary(
                            y => y.Key, 
                            y => 
                        {
                            var transcriptLink = y.OrderBy(z => z.ConfidenceScore).First();
                            return new MapLink
                            {
                                ConfidenceScore = transcriptLink.ConfidenceScore,
                                Correlation = transcriptLink.Correlation,
                                TranscriptName = transcriptLink.TranscriptName,
                                TssName = x.Key,
                                GeneName = x.Key,
                                LinkLength = transcriptLink.LinkLength,
                                LocusName = transcriptLink.LocusName,
                                Strand = expression.Genes[x.Key].Strand,
                            };
                        }))));
        }

        public TssRegulatoryMap ConvertToGenes()
        {
            return Convert(this.Links.ToLookup(x => new MapLink.Tss(x.GeneName), x => x)
                .ToDictionary(
                    x => x.Key, 
                    x => new Dictionary<MapLink.Locus, MapLink>(x.ToLookup(y => y.LocusName, y => y).ToDictionary(
                            y => y.Key, 
                            y => 
                            {
                                var transcriptLink = x.OrderBy(z => z.ConfidenceScore).First();
                                return new MapLink
                                {
                                    ConfidenceScore = transcriptLink.ConfidenceScore,
                                    Correlation = transcriptLink.Correlation,
                                    TranscriptName = transcriptLink.TranscriptName,
                                    TssName = x.Key,
                                    GeneName = x.Key,
                                    LinkLength = transcriptLink.LinkLength,
                                    LocusName = transcriptLink.LocusName,
                                    Strand = transcriptLink.Strand,
                                };
                    }))));
        }

        /// <summary>
        /// Gets the TSS link degree distribution
        /// </summary>
        /// <returns>The degree distribution.</returns>
        public List<int> TssDegreeDistribution()
        {
            ////var distribution = this.Links
            ////    .GroupBy(x => x.TranscriptName) // Get links by transcript
            ////    .Select(x => x.Count())         // Count the transcript links
            ////    .GroupBy(x => x, x => x)        // Group the link counts
            ////    .ToDictionary(x => x.Key, x => x.Count());  // Count the frequency of link counts

            var distribution = this.TssDegree.ToLookup(x => x.Value, x => x).ToDictionary(x => x.Key, x => x.Count());
            var maxDegree = this.TssDegree.Values.Max();

            ////int maxDegree = distribution.Keys.Max();

            var degreeDistribution = new List<int>();
            for (int i = 0; i <= maxDegree; i++)
            {
                degreeDistribution.Add(distribution.ContainsKey(i) ? distribution[i] : 0);
            }

            return degreeDistribution;
        }

        /// <summary>
        /// Removes the highly connected tsses.
        /// </summary>
        /// <returns>The highly connected tsses.</returns>
        /// <param name="minCount">Minimum count.</param>
        public TssRegulatoryMap RemoveHighlyConnectedTsses(int minCount)
        {
            var tsses = new HashSet<MapLink.Tss>(this.TssDegree.Where(x => x.Value < minCount).Select(x => x.Key));

            return Convert(this.Where(x => tsses.Contains(x.Key)));
        }

        /// <summary>
        /// Gets the map for a given maximum range
        /// </summary>
        /// <returns>The range.</returns>
        /// <param name="maxRange">Max range.</param>
        public TssRegulatoryMap ApplyRange(int maxRange)
        {
            return Convert(this.Select(x => new 
            {
                Key = x.Key,
                Value = x.Value.Where(y => y.Value.AbsLinkLength <= maxRange).ToDictionary(y => y.Key, y => y.Value)
            })
                .Where(x => x.Value.Count > 0)
                .ToDictionary(x => x.Key, x => x.Value));
        }

        /// <summary>
        /// Creates a map with links shorter than the range removed
        /// </summary>
        /// <returns>The range.</returns>
        /// <param name="minRange">Minimum range.</param>
        public TssRegulatoryMap RemoveRange(int minRange)
        {
            return Convert(this.Select(x => new 
            {
                Key = x.Key,
                Value = x.Value.Where(y => y.Value.AbsLinkLength >= minRange).ToDictionary(y => y.Key, y => y.Value)
            })
                .Where(x => x.Value.Count > 0)
                .ToDictionary(x => x.Key, x => x.Value));
        }

        /// <summary>
        /// Removes Loci fromthe range relative to the TSS. May be a fixed distance or a relative upstream/downstream distance
        /// </summary>
        /// <returns>The range.</returns>
        /// <param name="minRange">Minimum range.</param>
        public TssRegulatoryMap RemoveRange(List<int> minRange)
        {
            if (minRange.Count == 1)
            {
                return this.RemoveRange(minRange[0]);
            }
            else if (minRange.Count == 2)
            {
                return Convert(this.Select(x => new 
                {
                    Key = x.Key,
                    Value = x.Value.Where(y => 
                        (y.Value.LocusPosition < 0 && y.Value.LocusPosition < minRange[0]) ||
                        (y.Value.LocusPosition >= 0 && y.Value.LocusPosition > minRange[1]))
                        .ToDictionary(y => y.Key, y => y.Value)
                    })
                    .Where(x => x.Value.Count > 0)
                    .ToDictionary(x => x.Key, x => x.Value));
            }
            else
            {
                throw new Exception("Invalid number genomic range: " + string.Join(", ", minRange));
            }
        }

        /// <summary>
        /// Gets the nearest neighbor map for the nearest N neighbors
        /// </summary>
        /// <returns>The nearest neighbor map.</returns>
        /// <param name="n">Number of nearest neighbors.</param>
        public TssRegulatoryMap GetNearestNeighborMap(int n)
        {
            return new TssRegulatoryMap(this.Select(
                x => x.Value.Values.OrderBy(y => y.AbsLinkLength).Take(n))
                .SelectMany(x => x));
        }


        /// <summary>
        /// Gets the map for given threshold.
        /// </summary>
        /// <returns>The map for threshold.</returns>
        /// <param name="threshold">Confidence threshold.</param>
        public TssRegulatoryMap ApplyThreshold(double threshold)
        {
            return Convert(this.Select(x => new
            {
                Key = x.Key,
                Value = x.Value.Where(y => y.Value.ConfidenceScore <= threshold).ToDictionary(y => y.Key, y => y.Value)
            })
                .Where(x => x.Value.Count > 0)
                .ToDictionary(x => x.Key, x => x.Value));
        }

        /// <summary>
        /// Gets the map for given threshold.
        /// </summary>
        /// <returns>The map for threshold.</returns>
        /// <param name="threshold">Confidence threshold.</param>
        public TssRegulatoryMap ApplyThresholdExclusive(double threshold)
        {
            return Convert(this.Select(x => new
            {
                Key = x.Key,
                Value = x.Value.Where(y => y.Value.ConfidenceScore < threshold).ToDictionary(y => y.Key, y => y.Value)
            })
                .Where(x => x.Value.Count > 0)
                .ToDictionary(x => x.Key, x => x.Value));
        }

        /// <summary>
        /// Applies the correlation.
        /// </summary>
        /// <returns>A map with a minimum correlation.</returns>
        /// <param name="correlation">Correlation threshols.</param>
        public TssRegulatoryMap ApplyCorrelation(double correlation)
        {
            return Convert(this.Select(x => new
            {
                Key = x.Key,
                Value = x.Value.Where(y => y.Value.Correlation > correlation).ToDictionary(y => y.Key, y => y.Value)
            }) 
                .Where(x => x.Value.Count > 0)
                .ToDictionary(x => x.Key, x => x.Value));
        }

        /// <summary>
        /// Gets the map for given threshold.
        /// </summary>
        /// <returns>The map for threshold.</returns>
        /// <param name="threshold">Confidence threshold.</param>
        public TssRegulatoryMap RemoveThreshold(double threshold)
        {
            return Convert(this.Select(x => new
            {
                Key = x.Key,
                Value = x.Value.Where(y => y.Value.ConfidenceScore >= threshold).ToDictionary(y => y.Key, y => y.Value)
            })
                .Where(x => x.Value.Count > 0)
                .ToDictionary(x => x.Key, x => x.Value));
        }


        /// <summary>
        /// Creates a Locus-centric NN map included only TSSes which are targeted by a Locus
        /// </summary>
        /// <returns>The centric NN map.</returns>
        public LocusRegulatoryMap LocusCentricNNMap()
        {
            return new LocusRegulatoryMap(this.Reverse().GetNearestNeighborMap(1).Links);
        }


        /// <summary>
        /// Reverse this instance.
        /// </summary>
        /// <returns>A regulatory map with Loci as keys</returns>
        public LocusRegulatoryMap Reverse()
        {
            return new LocusRegulatoryMap(this.Links
                    .ToLookup(x => x.LocusName, x => x)
                    .ToDictionary(x => x.Key, x => x.ToDictionary(y => y.TranscriptName, y => y)));
        }

        /// <summary>
        /// Samples the map. If there are 5000 sample failures, then the number of length bins is reduced by 2 and more 
        /// sampling is redone to add more maps.
        /// </summary>
        /// <returns>The map.</returns>
        /// <param name="basis">Map basis to mimic in the sampled map.</param>
        /// <param name="lengthBins">Length bins.</param>
        /// <param name="sampleCount">Number of sampled maps to generate</param>
        public TssRegulatoryMap[] SampleMap(TssRegulatoryMap basis, int lengthBins, int sampleCount)
        {
            List<TssRegulatoryMap> sampleMaps = new List<TssRegulatoryMap>();

            while (lengthBins > 0 && sampleMaps.Count < sampleCount)
            {
                var newRound = this.SampleMapInternal(basis, lengthBins, sampleCount - sampleMaps.Count);
                sampleMaps.AddRange(newRound);
                lengthBins -= 2;
            }

            if (sampleMaps.Count != sampleCount)
            {
                throw new Exception(string.Format("Tried to sample {0} maps but only obtained {1}", sampleCount, sampleMaps.Count));
            }

            return sampleMaps.ToArray();
        }


        /// <summary>
        /// Samples the map.
        /// </summary>
        /// <returns>The map.</returns>
        /// <param name="basis">Basis map to emulate in sampling.</param>
        /// <param name="lengthBins">Length bins.</param>
        /// <param name="sampleCount">Number of samples to create.</param>
        private List<TssRegulatoryMap> SampleMapInternal(TssRegulatoryMap basis, int lengthBins, int sampleCount)
        {
            var tssBasisDegree = basis.ToDictionary(x => x.Key, x => x.Value.Count);

            var basisLinksByLength = basis.Links
                .OrderBy(x => x.LocusPosition).ToList();

            ////var basisUpstreamLinks = basisLinksByLength.Where(x => x.LocusPosition < 0).OrderBy(x => -x.LocusPosition).ToList();
            ////var basisDownstreamLinks = basisLinksByLength.Where(x => x.LocusPosition >= 0).ToList();

            double basisBinSize = (double)basisLinksByLength.Count / lengthBins;
            var boundaryIndices = Stats.Sequence(basisBinSize, basisLinksByLength.Count - basisBinSize, basisBinSize);

            var boundaryLengths = boundaryIndices
                .Select(i => (basisLinksByLength[i - 1].LocusPosition + basisLinksByLength[i].LocusPosition) / 2)
                .ToList();

            ////Console.WriteLine("Boundary count = {0};\nBin count = {1}\nBoundaries = ({2})", boundaryLengths.Count, lengthBins, string.Join(", ", boundaryLengths));

            Func<int, int> getBin = (length) =>
            {
                if (length < boundaryLengths[0])
                {
                    return 0;
                }

                if (length >= boundaryLengths.Last())
                {
                    return lengthBins - 1;
                }

                for (int i = 1; i < boundaryLengths.Count; i++)
                {
                    if (length >= boundaryLengths[i - 1] && length < boundaryLengths[i])
                    {
                        return i;
                    }
                }

                throw new Exception("Could not find a bin for length " + length);
            };

            // Bin all links to a tss in the basis map
            var binnedLinks = this.Links
                .Where(x => tssBasisDegree.ContainsKey(x.TranscriptName))
                .Select(x => new 
                {
                    Bin = getBin(x.LocusPosition),
                    Link = x
                })
                .ToLookup(x => x.Bin, x => x.Link)
                .ToDictionary(x => x.Key, x => x.ToList());

            // Construct bins of all links for each tss separately
            var perTssBinnedBasisLinks = binnedLinks.Select(x => x.Value
                .Select(y => new
                {
                    Bin = x.Key,
                    TranscriptName = y.TranscriptName,
                    Link = y
                }))
                .SelectMany(x => x)
                .ToLookup(x => x.TranscriptName, x => x)
                .ToDictionary(
                    x => x.Key, 
                    x => x.ToLookup(y => y.Bin, y => y.Link).ToDictionary(y => y.Key, y => y.ToList()));

            // Bin the basis map
            var binnedBasisLinks = basis.Links
                .Where(x => tssBasisDegree.ContainsKey(x.TranscriptName))
                .Select(x => new 
                {
                    Bin = getBin(x.LocusPosition),
                    Link = x
                })
                .ToLookup(x => x.Bin, x => x.Link)
                .ToDictionary(x => x.Key, x => x.ToList());

            ////Console.WriteLine("Sampleable bin counts\n{0}", string.Join("\n", binnedLinks.OrderBy(x => x.Key).Select(x => x.Key + "\t" + x.Value.Count)));
            ////Console.WriteLine("Basis bin counts\n{0}", string.Join("\n", binnedBasisLinks.OrderBy(x => x.Key).Select(x => x.Key + "\t" + x.Value.Count)));

            Func<TssRegulatoryMap> makeMap = () =>
            {
                // Ensure at least one link per tss sampled evenly from the equal-occupancy bins
                var tssSampledLinkData = perTssBinnedBasisLinks.Select(x => 
                {
                    var bins = x.Value.Keys.ToList();
                    var bin = bins[Stats.RNG.Next(bins.Count)];
                    return new 
                    {
                        Bin = bin,
                        Link = x.Value[bin][Stats.RNG.Next(x.Value[bin].Count)]
                    };
                });

                var tssSampledLinks = new HashSet<MapLink>(tssSampledLinkData.Select(x => x.Link));
                var tssSampledBinCounts = tssSampledLinkData.ToLookup(x => x.Bin, x => x).ToDictionary(x => x.Key, x => x.Count());

                // Set aside tsses with only one link
                var tssSingleLinks = tssSampledLinks.Take(basis.OriginDegree.Count(x => x.Value == 1));
                var singleLinkTsses = new HashSet<MapLink.Tss>(tssSingleLinks.Select(x => x.TranscriptName));

                // Construct bins of all links for each tss separately
                /*var filteredPerTssBinnedBasisLinks = binnedLinks.Select(x => x.Value
                    .Where(y => !singleLinkTsses.Contains(y.TranscriptName))
                    .Select(y => new
                    {
                        Bin = x.Key,
                        TranscriptName = y.TranscriptName,
                        Link = y
                    }))
                    .SelectMany(x => x)
                    .ToLookup(x => x.TranscriptName, x => x)
                    .ToDictionary(
                        x => x.Key, 
                        x => x.ToLookup(y => y.Bin, y => y.Link).ToDictionary(y => y.Key, y => y.ToList()));

                var tssDoubleSampledLinkData = filteredPerTssBinnedBasisLinks.Select(x => 
                {
                    var bins = x.Value.Keys.ToList();
                    var bin = bins[Stats.RNG.Next(bins.Count)];
                    return new 
                    {
                        Bin = bin,
                        Link = x.Value[bin][Stats.RNG.Next(x.Value[bin].Count)]
                    };
                });

                var tssDoubleSampledLinks = new HashSet<LinkData>(tssDoubleSampledLinkData.Select(x => x.Link));
                var tssDoubleSampledBinCounts = tssDoubleSampledLinkData.ToLookup(x => x.Bin, x => x).ToDictionary(x => x.Key, x => x.Count());

                // Set aside tsses with only one link
                var tssDoubleLinks = tssDoubleSampledLinks.Take(basis.TssDegree.Count(x => x.Value == 2));
                var doubleLinkTsses = new HashSet<string>(tssDoubleLinks.Select(x => x.TranscriptName));*/

                // Remove pre-sampled links single links from bins as well as all links to TSSes that will have only one link
                var filteredBinnedLinks = binnedLinks
                    .ToDictionary(
                        x => x.Key, 
                        x => x.Value.Where(y => 
                            !tssSampledLinks.Contains(y) && !singleLinkTsses.Contains(y.TranscriptName)) 
                        ////&& !tssDoubleSampledLinks.Contains(y) && !doubleLinkTsses.Contains(y.TranscriptName)
                        .ToList());

                // Filter from bins
                var sampledLinks = binnedBasisLinks
                    .Select(x => Stats.SampleWithoutReplacement(
                        x.Value.Count  - 
                        (tssSampledBinCounts.ContainsKey(x.Key) ? tssSampledBinCounts[x.Key] : 0),
                        ////+ (tssDoubleSampledBinCounts.ContainsKey(x.Key) ? tssDoubleSampledBinCounts[x.Key] : 0)
                        filteredBinnedLinks[x.Key].Count)
                        .Select(i => filteredBinnedLinks[x.Key][i]))
                    .SelectMany(x => x)
                    .ToList();

                //// Remove from bin samples enough to account for pre-sampled tss links
                ////var filteredSampledLinks = Stats.SampleWithoutReplacement(basis.LinkCount - basis.Count, sampledLinks.Count).Select(i => sampledLinks[i]);

                var sampledMap = new TssRegulatoryMap(tssSampledLinks
                    ////.Concat(tssDoubleSampledLinks)
                    .Concat(sampledLinks));

                ////var sampleTssDegree = sampledMap.ToDictionary(x => x.Key, x => x.Value.Count);

                /*Console.WriteLine("Map tss and link count: {0}, {1} ", basis.Count, basis.LinkCount);
                Console.WriteLine("Tss degree = {{\n{0}\n}}", string.Join("\n", tssBasisDegree
                    .ToLookup(x => x.Value, x => x)
                    .ToDictionary(x => x.Key, x => x.Count())
                    .OrderBy(x => x.Key)
                    .Select(x => "\t" + x.Key + "\t" + x.Value)));
                Console.WriteLine("Sampledtss and link count: {0}, {1}", sampledMap.Count, sampledMap.LinkCount);
                Console.WriteLine("Sampled Tss degree = {{\n{0}\n}}", string.Join("\n", sampleTssDegree
                    .ToLookup(x => x.Value, x => x)
                    .ToDictionary(x => x.Key, x => x.Count())
                    .OrderBy(x => x.Key)
                    .Select(x => "\t" + x.Key + "\t" + x.Value)));*/

                return sampledMap;
            };

            var sampledMaps = new List<TssRegulatoryMap>();
            int count = 0;
            for (int i = 0; i < sampleCount; i++)
            {
                count++;
                TssRegulatoryMap map = null;
                bool success = true;
                try
                {
                    map = makeMap();
                }
                catch
                {
                    i--;
                    success = false;
                }

                if (count % 100 == 0)
                {
                    Console.Write(".");
                }

                if (count > 1000)
                {
                    Console.WriteLine("Could not create sample maps in 1,000 tries");
                    return sampledMaps;
                }

                if (success && map != null)
                {
                    sampledMaps.Add(map);
                }
            }

            return sampledMaps;
        }
    }
}

