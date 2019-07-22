using System;

namespace Genomics
{
    using System.Collections.Generic;
    using System.Linq;
   
    public class KmerCollection
    {
        public KmerCollection(int kmerLength)
        {
            this.KmerLength = kmerLength;
        }

        public NucleotideAlphabet Alphabet = new DnaAmbiguousAlphabet();

        public void AddSequence(string sequence, int index)
        {
            //var reverseSequence = this.Alphabet.ReverseComplement(sequence);

            for (int i = 0; i < sequence.Length - this.KmerLength; i++)
            {
                var forwardKmer = sequence.Substring(i, this.KmerLength);
                var reverseKmer = this.Alphabet.ReverseComplement(forwardKmer);//reverseSequence.Substring(sequence.Length - this.KmerLength, this.KmerLength);

                bool foundForward = this.KmerCounts.ContainsKey(forwardKmer);
                bool foundReverse = this.KmerCounts.ContainsKey(reverseKmer);

                if (!foundForward && !foundReverse)
                {
                    this.KmerCounts.Add(forwardKmer, 1);
                    this.KmerSources.Add(forwardKmer, new List<KmerSource> { new KmerSource { Index = index, Position = i } });
                }
                else if (foundForward)
                {
                    this.KmerCounts[forwardKmer]++;
                    this.KmerSources[forwardKmer].Add(new KmerSource { Index = index, Position = i });
                }
                else
                {
                    this.KmerCounts[reverseKmer]++;
                    this.KmerSources[reverseKmer].Add(new KmerSource { Index = index, Position = -i });
                }
            }
        }
         
        public int KmerLength { get; private set; }

        public void ReportMedianSpan()
        {
            var spanSet = new HashSet<int>();

            Console.WriteLine("Kmer\tCount\tAddedCount\tCoverage");

            var counts = this.KmerCounts.Select(x => x).ToDictionary(x => x.Key, x => x.Value);


            var bestKmers = this.KmerCounts.OrderBy(x => -x.Value).ToList();

            var sequenceCount = this.KmerSources.Max(x => x.Value.Max(y => y.Index));

            int kmerCount = 0;
            while (spanSet.Count < sequenceCount && counts.Count > 0)
            {
                var bestKmer = counts.OrderBy(x => -x.Value).First();

                kmerCount++;
                var newSequences = new HashSet<int>(this.KmerSources[bestKmer.Key]
                    .Where(x => !spanSet.Contains(x.Index))
                    .Select(x => x.Index));

                var newSequenceKmers = this.KmerSources.Select(x => new 
                { 
                    Key = x.Key, 
                    Value = x.Value.Count(y => newSequences.Contains(y.Index))
                }).ToDictionary(x => x.Key, x => x.Value);

                var newCounts = counts.Select(x => new 
                {
                    Key = x.Key, 
                    Value = x.Value - (newSequenceKmers.ContainsKey(x.Key) ? newSequenceKmers[x.Key] : 0)
                }).Where(x => x.Value > 0)
                    .ToDictionary(x => x.Key, x => x.Value);


                Console.WriteLine(string.Join("\t", new string[]
                {
                    bestKmer.Key,
                    bestKmer.Value.ToString(),
                    newSequences.Count.ToString(),
                    spanSet.Count.ToString()
                }));

                spanSet.UnionWith(newSequences);

                counts = newCounts;
            }

            Console.WriteLine("Kmer count = " + kmerCount);


            /*Console.WriteLine("Sequence count: " + sequenceCount);

            for (int i = 0; i < bestKmers.Count && spanSet.Count < sequenceCount / 2; i++)
            {
                var newSequences = new HashSet<int>(this.KmerSources[bestKmers[i].Key]
                    .Where(x => !spanSet.Contains(x.Index))
                    .Select(x => x.Index));

            }*/
        }

        public Dictionary<string, int> KmerCounts
        {
            get
            {
                return this.kmerCounts;
            }
        }

        public Dictionary<string, List<KmerSource>> KmerSources
        {
            get
            {
                return this.kmerSources;
            }
        }

        private Dictionary<string, int> kmerCounts = new Dictionary<string, int>();
        private Dictionary<string, List<KmerSource>> kmerSources = new Dictionary<string, List<KmerSource>>();

        public struct KmerSource
        {
            public int Index { get; set; }
            public int Position { get; set; }
        }
    }
}

