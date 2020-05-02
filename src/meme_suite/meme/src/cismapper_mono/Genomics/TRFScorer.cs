//--------------------------------------------------------------------------------
// File: TRFScorer.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.IO;
    using Shared;
    using Tools;
    using Data;

    /// <summary>
    /// Utility methods for scoring or transforming ChIP-seq peak or other data at collections of genomic locations.
    /// </summary>
    public class TRFScorer
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.TRFScorer"/> class.
        /// </summary>
        public TRFScorer()
        {
        }

        /// <summary>
        /// Wrapper of indexed scoring data
        /// </summary>
        public class IndexedScores
        {
            /// <summary>
            /// Gets or sets scoring index.
            /// </summary>
            /// <value>The index.</value>
            public Dictionary<string, int> Index { get; set; }

            /// <summary>
            /// Gets or sets the scores.
            /// </summary>
            /// <value>The scores.</value>
            public List<double[]> Scores { get; set; }
        }

        public class ScoreEntry
        {
            public string Tss { get; set; }
            public List<double> PromoterScores { get; set; }
            public List<double> Locicores { get; set; }
            public int LocusCount { get; set; }
            public double Exp { get; set; }
        }

        /// <summary>
        /// Score type.
        /// </summary>
        public enum ScoreType
        {
            /// <summary>
            /// Scores at Locus regions
            /// </summary>
            Locus,

            /// <summary>
            /// Scores in the promoter region near the TSS
            /// </summary>
            Tss
        }


        /// <summary>
        /// Missing score action.
        /// </summary>
        public enum MissingScoreAction
        {
            /// <summary>
            /// Use raw scores.
            /// </summary>
            None,

            /// <summary>
            /// Fill missing scores with expected value
            /// </summary>
            Expected,
        }

        /// <summary>
        /// Scoring function.
        /// </summary>
        public enum ScoringFunction
        {
            /// <summary>
            /// Raw summation of all associated peak heights.
            /// </summary>
            Sum,

            /// <summary>
            /// The shuffled sum.
            /// </summary>
            ShuffledSum,

            /// <summary>
            /// Weighted summation of all associated peak heights.
            /// </summary>
            CorrelationWeightedSum,

            /// <summary>
            /// The shuffled correlation weighted sum.
            /// </summary>
            ShuffledCorrelationWeightedSum,

            /// <summary>
            /// The confidence weighted sum.
            /// </summary>
            ConfidenceWeightedSum,

            /// <summary>
            /// The shuffled confidence weighted sum.
            /// </summary>
            ShuffledConfidenceWeightedSum,

            /// <summary>
            /// The distance weighted function.
            /// </summary>
            DistanceWeightedSum,
        }
      
        /// <summary>
        /// Evaluates a genomic range.
        /// </summary>
        /// <returns>The range.</returns>
        /// <param name="items">Items.</param>
        /// <param name="start">Start.</param>
        /// <param name="end">End.</param>
        /// <param name="initialValue">Initial value.</param>
        /// <param name="PreStartTester">Pre start tester.</param>
        /// <param name="PreEndTester">Pre end tester.</param>
        /// <param name="Evaluator">Evaluator.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public static T EvaluateRange<T>(
            List<Genomics.Location> items, 
            int start, 
            int end, 
            T initialValue, 
            Func<Genomics.Location, int, bool> PreStartTester,
            Func<Genomics.Location, int, bool> PreEndTester,
            Func<T, int, Genomics.Location, T> Evaluator)
        {
            int i = FindRangeStart(items, start, PreStartTester);

            T val = initialValue;

            if (i > 0 && !PreStartTester(items[i-1], start))
            {
                Console.WriteLine("Pos {0}; start {1}; end {2}", items[i-1], start, end);
            }

            int length = end - start + 1;

            while (i < items.Count && PreEndTester(items[i], end))
            {
                //Console.WriteLine("Pos {0}; start {1}; end {2}", items[i].Mid, start, end);
                val = Evaluator(val, length, items[i]);
                i++;
            }
            
            return val;
        }

        /// <summary>
        /// Gets the closest location to a given position outside the given range
        /// </summary>
        /// <returns>The closes exclude range.</returns>
        /// <param name="items">Items.</param>
        /// <param name="position">Position.</param>
        /// <param name="range">Range.</param>
        public static Genomics.Location GetClosestExcludeRange(List<Genomics.Location> items, int position, int range)
        {
            List<Genomics.Location> filteredItems = items
                .Where(x => Math.Abs(x.Start - position) > range && Math.Abs(x.End - position) > range)
                .ToList();

            return GetClosestMid(filteredItems, position);
        }

        /// <summary>
        /// Intersects the two location sets returing the left items that overlap anything from the right
        /// </summary>
        /// <returns>The left.</returns>
        /// <param name="f1">F1.</param>
        /// <param name="f2">F2.</param>
        /// <param name="chr">Chr.</param>
		public static List<Genomics.Location> IntersectLeft(BedFile f1, IAnnotation f2, string chr, int buffer)
        {
            int maxFeatureSize = Math.Max(f1.MaxLocationSize.Max(x => x.Value), f2.MaxFeatureSize);

            int numFlankingBins = maxFeatureSize / (2 * BedFile.IndexSize) + 1;

            if (!f1.ChromosomeIndexedLocations.ContainsKey(chr) || 
                !f2.IndexedTranscriptLocations.ContainsKey(chr))
            {
                //Console.WriteLine("Invalid chromosome: " + chr);
                return new List<Location>();
            }

            return f1.ChromosomeIndexedLocations [chr].Select(x =>
                {
                    return x.Value
                        .Where(loc => Stats.Sequence(x.Key - numFlankingBins, x.Key + numFlankingBins, 1)
                        .Any(bin => f2.IndexedTranscriptLocations[chr].ContainsKey(bin) &&
                            f2.IndexedTranscriptLocations[chr][bin]
                                .Any(l2 => loc.Overlaps(buffer == 0 ? l2 : new Genomics.Location{ Start = l2.Start - buffer, End = l2.End + buffer}))))
                        .ToList();
                })
                .SelectMany(x => x)
                .ToList();
        }



        /// <summary>
        /// Intersects the two location sets returing the left items that overlap anything from the right
        /// </summary>
        /// <returns>The left.</returns>
        /// <param name="f1">F1.</param>
        /// <param name="f2">F2.</param>
        /// <param name="chr">Chr.</param>
        public static List<Genomics.Location> IntersectLeft(IAnnotation f1, BedFile f2, string chr, int buffer)
        {
            int maxFeatureSize = Math.Max(f2.MaxLocationSize.Max(x => x.Value), f1.MaxFeatureSize);

            int numFlankingBins = maxFeatureSize / (2 * BedFile.IndexSize) + 1;

            if (!f1.IndexedTranscriptLocations.ContainsKey(chr) || 
                !f2.ChromosomeIndexedLocations.ContainsKey(chr))
            {
                //Console.WriteLine("Invalid chromosome: " + chr);
                return new List<Location>();
            }

            return f1.IndexedTranscriptLocations [chr].Select(x =>
            {
                return x.Value
                    .Where(loc => Stats.Sequence(x.Key - numFlankingBins, x.Key + numFlankingBins, 1)
                        .Any(bin => f2.ChromosomeIndexedLocations[chr].ContainsKey(bin) &&
                            f2.ChromosomeIndexedLocations[chr][bin]
                            .Any(l2 => loc.Overlaps(buffer == 0 ? l2 : new Genomics.Location{ Start = l2.Start - buffer, End = l2.End + buffer}))))
                    .ToList();
            })
                .SelectMany(x => x)
                .ToList();
        }

        /// <summary>
        /// Intersects the two location sets returing the left items that overlap anything from the right
        /// </summary>
        /// <returns>The left.</returns>
        /// <param name="f1">F1.</param>
        /// <param name="f2">F2.</param>
        /// <param name="chr">Chr.</param>
        public static List<Genomics.Location> IntersectLeft(BedFile f1, BedFile f2, string chr, int buffer)
        {
            int maxFeatureSize = Math.Max(f1.MaxLocationSize.Max(x => x.Value), f2.MaxLocationSize.Max(x => x.Value));

            int numFlankingBins = maxFeatureSize / (2 * BedFile.IndexSize) + 1;

            if (!f1.ChromosomeIndexedLocations.ContainsKey(chr) || 
                !f2.ChromosomeIndexedLocations.ContainsKey(chr))
            {
                //Console.WriteLine("Invalid chromosome: " + chr);
                return new List<Location>();
            }

            return f1.ChromosomeIndexedLocations [chr].Select(x =>
            {
                return x.Value
                    .Where(loc => Stats.Sequence(x.Key - numFlankingBins, x.Key + numFlankingBins, 1)
                        .Any(bin => f2.ChromosomeIndexedLocations[chr].ContainsKey(bin) &&
                                    f2.ChromosomeIndexedLocations[chr][bin]
                            .Any(l2 => loc.Overlaps(buffer == 0 ? l2 : new Genomics.Location{ Start = l2.Start - buffer, End = l2.End + buffer}))))
                    .ToList();
            })
                .SelectMany(x => x)
                .ToList();
        }

		/// <summary>
		/// Gets locations that overlap a given item
		/// </summary>
		/// <returns>The overlaps.</returns>
		/// <param name="item">Item.</param>
		/// <param name="indexedLocations">Indexed locations.</param>
		/// <param name="binSize">Bin size.</param>
		/// <param name="maxFeatureSize">Max feature size.</param>
		public static List<Genomics.Location> GetOverlaps(
            Location item, 
            Dictionary<string, Dictionary<int, List<Genomics.Location>>> indexedLocations, 
            int binSize, 
            int maxFeatureSize,
            Func<Location, Location, bool> overlap = null)
		{
			maxFeatureSize = Math.Max(maxFeatureSize, item.End - item.Start);
			int numFlankingBins = maxFeatureSize / (2 * BedFile.IndexSize) + 1;
            int itemStartBin = item.Start / binSize;
            int itemEndBin = item.End / binSize;

			if (!indexedLocations.ContainsKey(item.Chromosome))
			{
				return null;
			}

            var overlappingItems = Stats.Sequence(itemStartBin - numFlankingBins, itemEndBin + numFlankingBins, 1)
				.Where(bin => indexedLocations[item.Chromosome].ContainsKey(bin))
				.Select(bin => indexedLocations[item.Chromosome][bin]
                    .Where(x => overlap == null ? item.Overlaps(x) : overlap(item, x)))
				.Where(x => x.Any())
				.SelectMany(x => x)
				.ToList();

            return overlappingItems;
		}

        /// <summary>
        /// Gets the closest.
        /// </summary>
        /// <returns>The closest.</returns>
        /// <param name="items">Items.</param>
        /// <param name="selector">Selector.</param>
        /// <param name="position">Position.</param>
        public static Genomics.Location GetClosest(List<Genomics.Location> items, Func<Genomics.Location, int> selector, int position)
        {
            int index = FindRangeStart(items, position, (item, start) => selector(item) < start);

            if (index == 0)
            {
                return items [0];
            }

            if (index >= items.Count - 1)
            {
                return items [items.Count - 1];
            }

            int leftDist = Math.Abs(selector(items [index - 1]) - position);
            int centerDist = Math.Abs(selector(items [index]) - position);
            int rightDist = Math.Abs(selector(items [index + 1]) - position);

            if (leftDist < centerDist && leftDist < rightDist)
            {
                return items [index - 1];
            }

            if (centerDist < rightDist)
            {
                return items[index];
            }

            return items[index + 1];
        }

        /// <summary>
        /// Gets the closest upstream or downstream location.
        /// Items must be sorted by start or end if searching downstream or upstream, respectively.
        /// </summary>
        /// <returns>The closest upstream and downstream.</returns>
        /// <param name="position">Position.</param>
        /// <param name="items">Items.</param>
        /// <param name="strand">Strand.</param>
        public static Location GetClosestDirectional(
            int position, 
            List<Location> items, 
            bool upstream)
        {
            int index = FindRangeStart(items, position, (item, start) => item.Start < start);

            if (index == 0)
            {
                return upstream ? null : items[0];
            }

            if (index >= items.Count - 1)
            {
                return upstream ? items [items.Count - 1] : null;
            }

            return upstream ? items[index - 1] : items[index];
        }

        /// <summary>
        /// Get closes location on a chromosome to a given position. 
        /// The items must be sorted by position.
        /// </summary>
        /// <returns>The closest.</returns>
        /// <param name="items">Items.</param>
        /// <param name="position">Position.</param>
        public static Genomics.Location GetClosestMid(List<Genomics.Location> items, int position)
        {
            return GetClosest(items, x => x.Mid, position);
        }

        /// <summary>
        /// Get closes location to a given position
        /// </summary>
        /// <returns>The closest.</returns>
        /// <param name="items">Items.</param>
        /// <param name="position">Position.</param>
        public static Genomics.Location GetClosestMinDistance(List<Genomics.Location> items, int position, int minDistance)
        {
            int index = FindRangeStart(items, position, (item, start) => item.Mid < start);
            
            if (index == 0)
            {
                return items [0];
            }
            
            if (index >= items.Count - 1)
            {
                return items [items.Count - 1];
            }
            
            int leftDist = Math.Abs(items [index - 1].Mid - position);
            int centerDist = Math.Abs(items [index].Mid - position);
            int rightDist = Math.Abs(items [index + 1].Mid - position);
            
            if (leftDist < centerDist && leftDist < rightDist)
            {
                return items [index - 1];
            }
            
            if (centerDist < rightDist)
            {
                return items[index];
            }
            
            return items[index + 1];
        }

		/// <summary>
		/// Creates the tss scores.
		/// </summary>
		/// <returns>The tss scores.</returns>
		/// <param name="tssData">Tss data.</param>
		/// <param name="peakFile">Peak file.</param>
		/// <param name="D0">D0.</param>
		/// <param name="maxRange">Max range.</param>
        public static Dictionary<string, double> CreateTssScores(
            Dictionary<string, Genomics.Location> tssData,
            BedFile peakFile,
            double D0,
            List<int> maxRange)
        {
            if (maxRange.Count > 1 && D0 != 0)
            {
                throw new Exception("Scoring with exponential decay only allows single promoter range value");
            }

            int upstreamDistance = maxRange[0];
            int dowstreamDistance = maxRange[0];
            if (maxRange.Count == 2)
            {
                dowstreamDistance = maxRange[1];
            }

            return CreateScores(
                tssData,
                peakFile,
                (items, location) => TRFScorer.EvaluateRange<double>(
                    items,
                    location.Start - (location.Strand == "+" ? upstreamDistance : dowstreamDistance),
                    location.End   + (location.Strand == "+" ? dowstreamDistance : upstreamDistance),
                    0,
                    (item, start) => item.Mid < start,
                    (item, end) => item.Mid <= end,
                    (aggregateScore, length, peak) => aggregateScore + (D0 == 0 ? 
                        peak.Score : 
                        Math.Exp(-(double)Math.Abs(peak.Mid - location.Start) / D0) * peak.Score))
				);
        }

		/// <summary>
		/// Creates the tss scores.
		/// </summary>
		/// <returns>The tss scores.</returns>
		/// <param name="tssData">Tss data.</param>
		/// <param name="peakFile">Peak file.</param>
		/// <param name="D0">D0.</param>
		/// <param name="maxRange">Max range.</param>
		public static Dictionary<string, double> CreateTssChangeScores(
			Dictionary<string, Genomics.Location> tssData,
			BedFile peakFile1,
			BedFile peakFile2,
			double D0,
            List<int> maxRange)
		{
            if (maxRange.Count > 1 && D0 != 0)
            {
                throw new Exception("Scoring with exponential decay only allows single promoter range value");
            }

			var scoreset1 = CreateTssScores(tssData, peakFile1, D0, maxRange);
			var scoreset2 = CreateTssScores(tssData, peakFile2, D0, maxRange);

			return scoreset1.Keys.Select(x => new
				{
					Key = x,
					Value = scoreset1[x] / scoreset2[x]
				})
				.ToDictionary(x => x.Key, x => x.Value);
		}
        
		/// <summary>
		/// Creates the Locus scores for all Loci
		/// </summary>
		/// <returns>The Locus scores.</returns>
		/// <param name="LocusFile">Locus file.</param>
		/// <param name="peakFile">Peak file.</param>
		/// <param name="normalize">If set to <c>true</c> normalize.</param>
        public static Dictionary<string, double> CreateLocicores(
			BedFile LocusFile,
            BedFile peakFile,
            bool normalize)
        {
            var temp = CreateScores(
                LocusFile.Locations,
                peakFile,
                (items, location) => TRFScorer.EvaluateRange<double>(
                    items, 
                    location.Start, 
                    location.End - 1,
                    0,
                    (item, start) => item.End < start,
                    (item, end) => item.Start <= end,
                    (aggregateScore, length, peak) => aggregateScore + peak.Score / (normalize ? length : 1)));
            return temp;
           
        }

        /*public static Dictionary<string, double> CreateLocicores(
            BedFile LocusFile,
            BedFile peakFile)
        {
            CreateScores(
                LocusFile.Locations,
                peakFile,
                (items, location) => GetOverlaps(location, peakFile.ChromosomeIndexedLocations, BedFile.IndexSize, peakFile.MaxLocationSize).Sum(x => x.Score));
        }*/

		/// <summary>
		/// Creates the Locus scores for all Loci
		/// </summary>
		/// <returns>The Locus scores.</returns>
		/// <param name="LocusFile">Locus file.</param>
        /// <param name = "peakFile1"></param>
        /// <param name = "peakFile2"></param>
		/// <param name="normalize">If set to <c>true</c> normalize.</param>
		public static Dictionary<string, double> CreateLocusChangeScores(
			BedFile LocusFile,
			BedFile peakFile1,
			BedFile peakFile2,
			bool normalize)
		{
			var scoreset1 = CreateLocicores(LocusFile, peakFile1, normalize);
			var scoreset2 = CreateLocicores(LocusFile, peakFile2, normalize);

			return scoreset1.Keys.Select(x => new
				{
					Key = x,
					Value = scoreset1[x] / scoreset2[x]
				})
				.ToDictionary(x => x.Key, x => x.Value);
		}

		/// <summary>
		/// Creates a single Locus score.
		/// </summary>
		/// <returns>The Locus score.</returns>
		/// <param name="location">Location.</param>
		/// <param name="peakFile">Peak file.</param>
		public static double CreateLocicore(
			Genomics.Location location,
			BedFile peakFile)
		{
			if (peakFile.ChromosomeOrderedLocations.ContainsKey(location.Chromosome))
			{
                int n = 0;
				return TRFScorer.EvaluateRange<double>(
					peakFile.ChromosomeOrderedLocations[location.Chromosome], 
					location.Start, 
					location.End - 1,
					0,
					(item, start) => item.End < start,
					(item, end) => item.Start <= end,
					(aggregateScore, length, peak) =>
					{
                        return (aggregateScore * (n-1) + peak.Score) / n;
					});
			}

			return 0.0;
		}
        
		/// <summary>
		/// Creates the scores.
		/// </summary>
		/// <returns>The scores.</returns>
		/// <param name="locations">Locations.</param>
		/// <param name="peakFile">Peak file.</param>
		/// <param name="Scorer">Scorer.</param>
        public static Dictionary<string, double> CreateScores(
            Dictionary<string, Genomics.Location> locations,
            BedFile peakFile,
            Func<List<Genomics.Location>, Genomics.Location, double> Scorer)
        {
            var peaks = peakFile.GetOrderedLocations(x => x.Mid);
            
            return
                (from location in locations
                 where peaks.ContainsKey(location.Value.Chromosome)
                 select new 
                 {
                    Name = location.Key,
                    Score = Scorer(peaks [location.Value.Chromosome], location.Value)
                }).ToDictionary(x => x.Name, x => x.Score);
        }

        /// <summary>
        /// Creates index of TRF score files for efficient loading
        /// </summary>
        /// <returns>The index.</returns>
        /// <param name="locations">Locations.</param>
        /// <param name="type">Type.</param>
        /// <param name="tissue">Cell line.</param>
        public static Dictionary<string, int> CreateIndex(Dictionary<string, Genomics.Location> locations, ScoreType scoreType, string dataSource)
        {
            var data = locations.Keys.Select((x, i) => new { Name = x, Index = i });

            using (TextWriter tw = Helpers.CreateStreamWriter(string.Format("../temp/GenomicFeatures/{0}.{1}.index", scoreType, dataSource)))
            {
                tw.WriteLine(string.Join("\n", data.Select(x => string.Format("{0}\t{1}\n", x.Name, x.Index))));
            }

            return data.ToDictionary(x => x.Name, x => x.Index);
        }

        /// <summary>
        /// Retrieves a score index for efficient loading of DRM or TSS scores
        /// </summary>
        /// <returns>The index.</returns>
        /// <param name="type">Type.</param>
        /// <param name="tissue">Cell line.</param>
        public static Dictionary<string, int> ReadIndex(ScoreType scoreType, string tissue)
        {
            using (TextReader tr = new StreamReader(string.Format("../temp/GenomicFeatures/{0}.{1}.index", scoreType, tissue)))
            {
                return tr.ReadToEnd().Split('\n')
                    .Where(line => !string.IsNullOrWhiteSpace(line))
                    .Select(line => 
                    {
                        var fields = line.Split('\t').Select(y => y.Trim()).ToArray();
                        return new
                        {
                            Name = fields[0],
                            Index = int.Parse(fields[1])
                        };
                    })
                    .ToDictionary(x => x.Name, x => x.Index);
            }
        }

		/// <summary>
		/// Writes the score name value pairs.
		/// </summary>
		/// <param name="fileName">File name.</param>
		/// <param name="scores">Scores.</param>
		/// <param name="usedItems">Used items.</param>
		/// <typeparam name="ValueType">The 1st type parameter.</typeparam>
        public static void WriteScoreNameValuePairs<TValueType>(string fileName, Dictionary<string, TValueType> scores, HashSet<string> usedItems)
        {
            using (TextWriter tw = Helpers.CreateStreamWriter(fileName))
            {
                if (usedItems == null)
                {
                    tw.WriteLine(scores.Count);
                    tw.WriteLine(string.Join("\n", scores.Select(x => string.Format("{0}\t{1}", x.Key, x.Value))));
                }
                else
                {
                    var data = scores.Where(x => usedItems.Contains(x.Key));
                    tw.WriteLine(data.Count());
                    tw.WriteLine(string.Join("\n", data.Select(x => string.Format("{0}\t{1}", x.Key, x.Value))));
                }
            }
        }

        /// <summary>
        /// Reads the score name value pairs.
        /// </summary>
        /// <returns>The score name value pairs.</returns>
        /// <param name="index">Index.</param>
        /// <param name="fileName">File name.</param>
        /// <param name="ValueParser">Value parser.</param>
        /// <typeparam name="TValueType">The 1st type parameter.</typeparam>
        public static TValueType[] ReadScoreNameValuePairs<TValueType>(Dictionary<string, int> index, string fileName, Func<string, TValueType> ValueParser)
        {
            using (TextReader tr = new StreamReader(fileName))
            {
                tr.ReadLine();

                int itemCount = index.Count();
                GC.Collect();

                TValueType[] array = new TValueType[itemCount];

                int tabIndex = 0;
                int newLineIndex = 0;
                int lineStartIndex = 0;
                string lines = tr.ReadToEnd();
                int i;

                do
                {
                    tabIndex = lines.IndexOf('\t', tabIndex + 1);

                    if (tabIndex > 0)
                    {
                        newLineIndex = lines.IndexOf('\n', tabIndex);

                        // string t = lines.Substring(lineStartIndex, tabIndex - lineStartIndex );
                        i = index[lines.Substring(lineStartIndex, tabIndex - lineStartIndex )];
                        array[i] = ValueParser(lines.Substring(tabIndex + 1, newLineIndex - tabIndex - 1));

                        lineStartIndex = newLineIndex + 1;
                    }
                } while (tabIndex > 0 && newLineIndex > 0);

                return array;
            }
        }

		/// <summary>
		/// Reads the score name value pairs.
		/// </summary>
		/// <returns>The score name value pairs.</returns>
		/// <param name="fileName">File name.</param>
		/// <param name="ValueParser">Value parser.</param>
		/// <typeparam name="ValueType">The 1st type parameter.</typeparam>
        public static Dictionary<string, ValueType> ReadScoreNameValuePairs<ValueType>(string fileName, Func<string, ValueType> ValueParser)
        {
            using (TextReader tr = new StreamReader(fileName))
            {
                tr.ReadLine(); // Dispose of count
                return tr.ReadToEnd().Split('\n').Where(line => !string.IsNullOrWhiteSpace(line)).Select(line => 
                {
                    string[] fields = line.Split('\t');
                    return new
                            {
                                Key = fields [0],
                                Value = ValueParser(fields [1])
                            };
                })
                    .ToDictionary(x => x.Key, x => x.Value);
            }
        }


        public static void WriteSingleIndexedScore(
            ScoreType scoreType,
            string peakName, 
            string peakFile, 
            string tissue, 
            string rnaSource,
            Dictionary<string, Location> locations)
        {
            var peaks = new BedFile(peakFile, BedFile.Bed6Plus4Layout);

            TRFScorer.CreateIndex(locations, scoreType, string.Join(".", new List<string> { peakName, tissue, rnaSource }));

            string baseFileName = Path.GetFileNameWithoutExtension(peakFile);

            TRFScorer.WriteScoreNameValuePairs(
                Utilities.GetLocicoreFileName(baseFileName, Utilities.TFScoreSource(peakName, tissue, rnaSource)),
                TRFScorer.CreateScores(
                    locations,
                    peaks,
                    (items, location) => TRFScorer.GetOverlaps(
                        location, 
                        peaks.ChromosomeIndexedLocations, 
                        BedFile.IndexSize, 
                        peaks.MaxLocationSize[location.Chromosome])
                    .Sum(x => double.Parse(x.Data[BedFile.Bed6Plus4Layout.SignalValue]))), 
                null);
        }

        public static void WriteScoreEntries(string fileName, List<TRFScorer.ScoreEntry> scores, List<string> chipSets)
        {
            Func<string, List<double>, List<double>, double, string[]> tableLine = (tss, tssScores, Locicores, expression) =>
            {
                return new string[] { tss }
                    .Concat(tssScores != null ? tssScores.Select(x => x.ToString()) : new List<string>())
                    .Concat(Locicores != null ? Locicores.Select(x => x.ToString()) : new List<string>())
                    .Concat(new string[] { expression.ToString() })
                    .ToArray();
            };

            var promoterColumns = scores.First().PromoterScores != null ? 
                chipSets.Select(x => string.Format("Tss_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0])) :
                new List<string>();

            var LocusColumns = scores.First().Locicores != null ? 
                chipSets.Select(x => string.Format("Drm_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0])) :
                new List<string>();

            Tables.ToNamedTsvFile(
                fileName,
                scores.Select(x => tableLine(x.Tss, x.PromoterScores, x.Locicores, x.Exp)),
                promoterColumns.Concat(LocusColumns).Concat(new string[] { "esc" }));
        }

        /// <summary>
        /// Creates the tss Locus and combined table.
        /// </summary>
        /// <returns>The tss Locus and combined table.</returns>
        /// <param name="LocusMap">Locus map.</param>
        /// <param name="tssScores">Tss scores.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="expression">Expression.</param>
        /// <param name="outFileName">Out file name.</param>
        /// <param name="chipSets">Chip sets.</param>
        public static int CreateTssLocusAndCombinedTable(
            Dictionary<MapLink.Tss, List<MapLink.Locus>> LocusMap, 
            TRFScorer.IndexedScores tssScores, 
            TRFScorer.IndexedScores Locicores,
            Dictionary<string, double> expression,
            string outFileName, 
            List<string> chipSets)
        {
            var rawScores = LocusMap
                .Where(mapEntry => expression.ContainsKey(mapEntry.Key) && tssScores.Index.ContainsKey(mapEntry.Key))
                .Select(mapEntry => 
                {
                    var atTss = tssScores.Scores.Select(chipData => chipData[tssScores.Index[mapEntry.Key]]).ToList();
                    var atLoci = Locicores.Scores.Select(chipData => mapEntry.Value.Where(x => Locicores.Index.ContainsKey(x)).Sum(Locus => chipData[Locicores.Index[Locus]])).ToList();

                    return new
                    {
                        Tss = mapEntry.Key,
                        PromoterScores = atTss,
                        Locicores = atLoci,
                        LocusCount = mapEntry.Value.Count(x => Locicores.Index.ContainsKey(x)),
                        Exp = expression[mapEntry.Key]
                    };
                })
                .ToList();

            var expectedLocicores = rawScores[0].Locicores.Select((x, i) => rawScores.Average(y => y.Locicores[i])).ToList();

            var scores = rawScores.Select(x => 
            {
                var atLoci = x.LocusCount == 0 ? expectedLocicores : x.Locicores;
                var total = x.PromoterScores.Select((f, i) => f + atLoci[i]);

                return new 
                {
                    x.Tss,
                    x.PromoterScores,
                    Locicores = atLoci,
                    CombinedScores = total,
                    x.Exp
                };
            })
                .ToList();

            Console.WriteLine("Number of TSSes: {0}", scores.Count);

            using (TextWriter tw = Helpers.CreateStreamWriter(outFileName))
            {
                // table format
                tw.WriteLine("{0}\t{1}\t{2}\tesc", 
                    string.Join("\t", chipSets.Select(x => string.Format("Tss_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))), 
                    string.Join("\t", chipSets.Select(x => string.Format("Drm_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))),
                    string.Join("\t", chipSets.Select(x => string.Format("TssDrm_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))));
                foreach (var gene in scores)
                {
                    tw.WriteLine(
                        "{0}\t{1}\t{2}\t{3}\t{4}",
                        gene.Tss, 
                        string.Join("\t", gene.PromoterScores), 
                        string.Join("\t", gene.Locicores), 
                        string.Join("\t", gene.CombinedScores), 
                        gene.Exp);
                }
            }

            return scores.Count;
        }

        /// <summary>
        /// Creates the tss Locus and combined table.
        /// </summary>
        /// <returns>The tss Locus and combined table.</returns>
        /// <param name="LocusMap">Locus map.</param>
        /// <param name="tssScores">Tss scores.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="expression">Expression.</param>
        /// <param name="outFileName">Out file name.</param>
        /// <param name="chipSets">Chip sets.</param>
        public static int CreateLocusTable(
            TRFScorer.MissingScoreAction action,
            TssRegulatoryMap LocusMap, 
            TRFScorer.IndexedScores Locicores,
            Dictionary<string, double> expression,
            string outFileName, 
            List<string> chipSets,
            Func<IEnumerable<MapLink>, TRFScorer.IndexedScores, int, double> Locicorer)
        {
            var rawScores = LocusMap
                .Where(mapEntry => expression.ContainsKey(mapEntry.Key))
                .Select(mapEntry => 
                {
                    var atLoci = Locicores.Scores
                        .Select((chipData, index) => Locicorer(mapEntry.Value.Values, Locicores, index))
                        .ToList();

                    return new TRFScorer.ScoreEntry
                    {
                        Tss = mapEntry.Key,
                        Locicores = atLoci,
                        LocusCount = mapEntry.Value.Values.Count(link => Locicores.Index.ContainsKey(link.LocusName)),
                        Exp = expression[mapEntry.Key]
                    };
                })
                .ToList();

            List<TRFScorer.ScoreEntry> scores;

            if (action == TRFScorer.MissingScoreAction.Expected)
            {
                var expectedLocicores = rawScores[0].Locicores.Select((x, i) => rawScores.Average(y => y.Locicores[i])).ToList();

                scores = rawScores.Select(x =>
                {
                    var atLoci = x.LocusCount == 0 ? expectedLocicores : x.Locicores;

                    return new TRFScorer.ScoreEntry
                    {
                        Tss = x.Tss,
                        Locicores = atLoci,
                        Exp = x.Exp
                    };
                })
                    .ToList();
            }
            else
            {
                scores = rawScores;
            }

            Console.WriteLine("Number of TSSes: {0}", scores.Count);

            using (TextWriter tw = Helpers.CreateStreamWriter(outFileName))
            {
                // table format
                tw.WriteLine("{0}\tesc", 
                    string.Join("\t", chipSets.Select(x => string.Format("Drm_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))));

                foreach (var gene in scores)
                {
                    tw.WriteLine(
                        "{0}\t{1}\t{2}",
                        gene.Tss, 
                        string.Join("\t", gene.Locicores), 
                        gene.Exp);
                }
            }

            return scores.Count;
        }

        /// <summary>
        /// Creates table of TSS and Locus data for regression using flat indexed tss and drm scores
        /// For TSSes with no Loci linked, the expected score is filled in
        /// </summary>
        /// <returns>The count of rows in the table.</returns>
        /// <param name="LocusMap">Locus map.</param>
        /// <param name="tssIndex">Tss index.</param>
        /// <param name="tssScores">Tss scores.</param>
        /// <param name="LocusIndex">Locus index.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="expression">Expression.</param>
        /// <param name="outFileName">Out file name.</param>
        /// <param name="chipSets">Chip sets.</param>
        public static int CreateTable(
            Dictionary<MapLink.Tss, List<MapLink.Locus>> LocusMap, 
            Dictionary<string, int> tssIndex,
            List<double[]> tssScores, 
            Dictionary<string, int> LocusIndex,
            List<double[]> Locicores, 
            Dictionary<string, double> expression,
            string outFileName, 
            List<string> chipSets)
        {
            var rawScores = LocusMap
                .Where(mapEntry => expression.ContainsKey(mapEntry.Key) && tssIndex.ContainsKey(mapEntry.Key))
                .Select(mapEntry => new
                {
                    Tss = mapEntry.Key,
                    PromoterScores = tssScores.Select(chipData => chipData[tssIndex[mapEntry.Key]]),
                    Locicores = Locicores.Select(chipData => mapEntry.Value.Where(x => LocusIndex.ContainsKey(x)).Sum(Locus => chipData[LocusIndex[Locus]])).ToList(),
                    LocusCount = mapEntry.Value.Count(x => LocusIndex.ContainsKey(x)),
                    Exp = expression[mapEntry.Key]
                })
                .ToList();

            var expectedLocicores = rawScores[0].Locicores.Select((x, i) => rawScores.Average(y => y.Locicores[i])).ToList();

            var scores = rawScores.Select(x => new 
            {
                x.Tss,
                x.PromoterScores,
                Locicores = x.LocusCount == 0 ? expectedLocicores : x.Locicores,
                x.Exp
            })
                .ToList();

            //Console.WriteLine("Number of TSSes: {0}", scores.Count);

            using (TextWriter tw = Helpers.CreateStreamWriter(outFileName))
            {
                // table format
                tw.WriteLine("{0}\t{1}\tesc", 
                    string.Join("\t", chipSets.Select(x => string.Format("Tss_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))), 
                    string.Join("\t", chipSets.Select(x => string.Format("Drm_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))));
                foreach (var gene in scores)
                {
                    tw.WriteLine("{0}\t{1}\t{2}\t{3}", gene.Tss, string.Join("\t", gene.PromoterScores), string.Join("\t", gene.Locicores), gene.Exp);
                }
            }

            return scores.Count;
        }


        /// <summary>
        /// convert non finite values to NA
        /// </summary>
        /// <returns>The filter double.</returns>
        /// <param name="x">The x coordinate.</param>
        private static string NAFilterDouble(double x)
        {
            if (double.IsNaN(x)) // || double.IsInfinity(x))
            {
                return "0";
            }

            return x.ToString();
        }

        /// <summary>
        /// Creates table of TSS and DRM data for regression using flat indexed tss and drm scores
        /// </summary>
        /// <returns>The table.</returns>
        /// <param name="drmMap">Drm map.</param>
        /// <param name="tssIndex">Tss index.</param>
        /// <param name="tssScores">Tss scores.</param>
        /// <param name="drmIndex">Drm index.</param>
        /// <param name="drmScores">Drm scores.</param>
        /// <param name="expression">Expression.</param>
        /// <param name="outFileName">Out file name.</param>
        /// <param name="chipSets">Chip sets.</param>
        public static int CreateTableFilterNA(
            Dictionary<MapLink.Tss, List<MapLink.Locus>> drmMap, 
            Dictionary<string, int> tssIndex,
            List<double[]> tssScores, 
            Dictionary<string, int> drmIndex,
            List<double[]> drmScores, 
            Dictionary<string, double> expression,
            string outFileName, 
            List<string> chipSets)
        {
            var v = drmMap
                .Where(mapEntry => 
                    expression.ContainsKey(mapEntry.Key) &&
                    tssIndex.ContainsKey(mapEntry.Key))
                .Select(mapEntry => new
                {
                    Tss = mapEntry.Key,
                    TssOnlyScore = tssScores.Select(chipData => chipData[tssIndex[mapEntry.Key]]),
                    Score = drmScores.Select(chipData => mapEntry.Value.Where(drm => drmIndex.ContainsKey(drm)).Sum(drm => chipData[drmIndex[drm]])),
                    Exp = expression[mapEntry.Key]
                });

            Console.WriteLine("Number of TSSes: {0}", v.Count());

            using (TextWriter tw = Helpers.CreateStreamWriter(outFileName))
            {
                // table format
                tw.WriteLine("{0}\t{1}\tesc", 
                    string.Join("\t", chipSets.Select(x => string.Format("Tss_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))), 
                    string.Join("\t", chipSets.Select(x => string.Format("Drm_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))));
                foreach (var gene in v)
                {
                    tw.WriteLine(
                        "{0}\t{1}\t{2}\t{3}", 
                        gene.Tss, 
                        string.Join("\t", gene.TssOnlyScore.Select(x => NAFilterDouble(x))), 
                        string.Join("\t", gene.Score.Select(x => NAFilterDouble(x))), 
                        gene.Exp);
                }
            }

            return v.Count();
        }


        /// <summary>
        /// Reads the single indexed score.
        /// </summary>
        /// <returns>The single indexed score.</returns>
        /// <param name="scoreType">Score type.</param>
        /// <param name="peakName">Peak name.</param>
        /// <param name="peakFile">Peak file.</param>
        /// <param name="tissue">Cell line.</param>
        /// <param name="rnaSource">Rna source.</param>
        public static IndexedScores ReadSingleIndexedScore(
            ScoreType scoreType,
            string peakName, 
            string peakFile, 
            string tissue,
            string rnaSource)
        {
            string scoreSource = Utilities.TFScoreSource(peakName, tissue, rnaSource);
            string baseFileName = Path.GetFileNameWithoutExtension(peakFile);

            var index = TRFScorer.ReadIndex(scoreType, scoreSource); 

            return new IndexedScores
            {
                Index = index,
                Scores = new List<double[]> 
                {
                    TRFScorer.ReadScoreNameValuePairs<double>(
                        index, 
                        Utilities.GetLocicoreFileName(baseFileName, scoreSource), 
                        x => double.Parse(x))
                }    
            };
        }

		/// <summary>
		/// Finds the range start.
		/// </summary>
		/// <returns>The range start.</returns>
		/// <param name="items">Items.</param>
		/// <param name="start">Start.</param>
		/// <param name="PreStartTester">Pre start tester.</param>
		/// <typeparam name="T">The 1st type parameter.</typeparam>
        private static int FindRangeStart<T>(List<T> items, int start, Func<T, int, bool> PreStartTester)
        {
            int lower = 0;
            int upper = items.Count;
            
            while (upper > lower)
            {
                int mid = (lower + upper) / 2;
                
                if (PreStartTester(items[mid], start))
                {
                    lower = mid + 1;
                } 
                else
                {
                    upper = mid;
                }
            }

            return upper;
        }

        /// <summary>
        /// Loci the summation scorer.
        /// </summary>
        /// <returns>The summation scorer.</returns>
        /// <param name="mapEntry">Map entry.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="index">Index.</param>
        public static double LociummationScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry, 
            TRFScorer.IndexedScores Locicores, 
            int index) 
        {
            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => Locicores.Scores[index][Locicores.Index[link.LocusName]]);
        }

        public static double LocihuffledSummationScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry, 
            TRFScorer.IndexedScores Locicores, 
            int index) 
        {
            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => Locicores.Scores[index][Stats.RNG.Next(Locicores.Scores[index].Length)]);
        }

        /// <summary>
        /// Creates a weighted Locus score
        /// </summary>
        /// <returns>The weighted scorer.</returns>
        /// <param name="mapEntry">Map entry.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="index">Index.</param>
        public static double LocusCorrelationWeightedScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry,
            TRFScorer.IndexedScores Locicores,
            int index) 
        {
            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => Math.Pow(link.Correlation, 2) * Locicores.Scores[index][Locicores.Index[link.LocusName]]);
        }

        public static double[] LocihuffledCorrelationWeightedScorerSample { get; set; }

        /// <summary>
        /// Creates a weighted Locus score
        /// </summary>
        /// <returns>The weighted scorer.</returns>
        /// <param name="mapEntry">Map entry.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="index">Index.</param>
        public static double LocihuffledCorrelationWeightedScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry,
            TRFScorer.IndexedScores Locicores,
            int index) 
        {
            if (LocihuffledCorrelationWeightedScorerSample == null)
            {
                LocihuffledCorrelationWeightedScorerSample = map.Links.Select(x => x.Correlation).ToArray();
            }

            var correlations = mapEntry.Select(x => x.Correlation).ToArray();

            //return mapEntry
            //    .Where(link => Locicores.Index.ContainsKey(link.LocusName))
            //    .Sum(link => Math.Pow(correlations[Stats.RNG.Next(correlations.Length)], 2) * Locicores.Scores[index][Locicores.Index[link.LocusName]]);

            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => Math.Pow(correlations[Stats.RNG.Next(correlations.Length)], 2) * Locicores.Scores[index][Locicores.Index[link.LocusName]]);
        }


        /// <summary>
        /// Creates a weighted Locus score
        /// </summary>
        /// <returns>The weighted scorer.</returns>
        /// <param name="mapEntry">Map entry.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="index">Index.</param>
        public static double LocusConfidenceWeightedScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry,
            TRFScorer.IndexedScores Locicores,
            int index) 
        {
            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => -Math.Log10(Math.Max(link.ConfidenceScore, 1e-20)) * Locicores.Scores[index][Locicores.Index[link.LocusName]]);
        }

        public static double[] LocihuffledConfidenceWeightedScorerSample { get; set; }

        /// <summary>
        /// Creates a weighted Locus score
        /// </summary>
        /// <returns>The weighted scorer.</returns>
        /// <param name="mapEntry">Map entry.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="index">Index.</param>
        public static double LocihuffledConfidenceWeightedScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry,
            TRFScorer.IndexedScores Locicores,
            int index) 
        {
            //var correlations = mapEntry.Select(x => x.Correlation).ToArray();
            if (LocihuffledConfidenceWeightedScorerSample == null)
            {
                LocihuffledConfidenceWeightedScorerSample = map.Links.Select(x => x.Correlation).ToArray();
            }

            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => -Math.Log10(
                    Math.Max(LocihuffledConfidenceWeightedScorerSample[Stats.RNG.Next(LocihuffledConfidenceWeightedScorerSample.Length)], 1e-20)) * 
                    Locicores.Scores[index][Locicores.Index[link.LocusName]]);
        }

        /// <summary>
        /// Loci the distance weighted scorer.
        /// </summary>
        /// <returns>The distance weighted scorer.</returns>
        /// <param name="mapEntry">Map entry.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="index">Index.</param>
        public static double LocusDistanceWeightedScorer(
            TssRegulatoryMap map,
            IEnumerable<MapLink> mapEntry,
            TRFScorer.IndexedScores Locicores,
            int index)
        {
            return mapEntry
                .Where(link => Locicores.Index.ContainsKey(link.LocusName))
                .Sum(link => Math.Exp(-link.AbsLinkLength / (1000)) * Locicores.Scores[index][Locicores.Index[link.LocusName]]);
        }

        /// <summary>
        /// Selects the scorer.
        /// </summary>
        /// <returns>The scorer.</returns>
        /// <param name="scoring">Requested scoring function.</param>
        public static Func<TssRegulatoryMap, IEnumerable<MapLink>, TRFScorer.IndexedScores, int, double> SelectScorer(ScoringFunction scoring)
        {
            switch (scoring)
            {
                case ScoringFunction.Sum:
                    return LociummationScorer;

                case ScoringFunction.ShuffledSum:
                    return LocihuffledSummationScorer;

                case ScoringFunction.CorrelationWeightedSum:
                    return LocusCorrelationWeightedScorer;

                case ScoringFunction.ShuffledCorrelationWeightedSum:
                    return LocihuffledCorrelationWeightedScorer;

                case ScoringFunction.DistanceWeightedSum:
                    return LocusDistanceWeightedScorer;

                case ScoringFunction.ConfidenceWeightedSum:
                    return LocusConfidenceWeightedScorer;

                case ScoringFunction.ShuffledConfidenceWeightedSum:
                    return LocihuffledConfidenceWeightedScorer;

                default:
                    return null;
            }
        }

        /// <summary>
        /// Creates the tss Locus and combined table.
        /// </summary>
        /// <returns>The tss Locus and combined table.</returns>
        /// <param name="LocusMap">Locus map.</param>
        /// <param name="tssScores">Tss scores.</param>
        /// <param name="Locicores">Locus scores.</param>
        /// <param name="expression">Expression.</param>
        /// <param name="outFileName">Out file name.</param>
        /// <param name="chipSets">Chip sets.</param>
        public static List<ScoreEntry> ScoreMap(
            MissingScoreAction action,
            TssRegulatoryMap LocusMap, 
            TRFScorer.IndexedScores tssScores,
            TRFScorer.IndexedScores Locicores,
            Dictionary<string, double> expression,
            Func<TssRegulatoryMap, IEnumerable<MapLink>, TRFScorer.IndexedScores, int, double> Locicorer)
        {
            var rawScores = LocusMap
                .Where(mapEntry => expression.ContainsKey(mapEntry.Key))
                .Select(mapEntry => 
                {
                    var atTss = tssScores != null ? 
                        tssScores.Scores.Select(chipData => chipData[tssScores.Index[mapEntry.Key]]).ToList() :
                        null;

                    var atLoci = Locicores != null ?
                        Locicores.Scores
                            .Select((chipData, index) => Locicorer(LocusMap, mapEntry.Value.Values, Locicores, index))
                            .ToList() :
                        null;

                    return new ScoreEntry
                    {
                        Tss = mapEntry.Key,
                        PromoterScores = atTss,
                        Locicores = atLoci,
                        LocusCount = mapEntry.Value.Values.Count(link => Locicores.Index.ContainsKey(link.LocusName)),
                        Exp = expression[mapEntry.Key]
                    };
                })
                .ToList();

            List<ScoreEntry> scores;

            if (action == MissingScoreAction.Expected)
            {
                var expectedLocicores = rawScores[0].Locicores.Select((x, i) => rawScores.Average(y => y.Locicores[i])).ToList();

                scores = rawScores.Select(x =>
                {
                    var atLoci = x.LocusCount == 0 ? expectedLocicores : x.Locicores;

                    return new ScoreEntry
                    {
                        Tss = x.Tss,
                        Locicores = atLoci,
                        Exp = x.Exp
                    };
                })
                    .ToList();
            }
            else
            {
                scores = rawScores;
            }

            return scores;
        }
    }
}

