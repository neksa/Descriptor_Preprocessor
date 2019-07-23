
//--------------------------------------------------------------------------------
// File: BedFile.cs
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
    using System.Text.RegularExpressions;
    using Shared;

    /// <summary>
    /// Basic bed file wrapper
    /// </summary>
    public class BedFile : IUnknown
    {
        /// <summary>
        /// Gets the size of the bin used for indexing chromosomal coordinates
        /// </summary>
        /// <value>The size of the index.</value>
        public static int IndexSize { get { return 100000; } }

        /// <summary>
        /// The locations of items on chromosomes.
        /// </summary>
        private Dictionary<string, Dictionary<string, Genomics.Location>> chromosomeLocations;

        /// <summary>
        /// The chromosome binned locations for use in indexed lookup
        /// </summary>
        //private Dictionary<string, Dictionary<int, Dictionary<string, Genomics.Location>>> chromosomeBinnedLocations;

        /// <summary>
        /// Locations in order by START on each chromosome
        /// </summary>
        private Dictionary<string, List<Location>> orderedLocations;

        /// <summary>
        /// Locations in order by END on each chromosome
        /// </summary>
        private Dictionary<string, List<Location>> orderedEndLocations;

        /// <summary>
        /// The indexed locations.
        /// </summary>
        private Dictionary<string, Dictionary<int, List<Location>>> indexedLocations;

        /// <summary>
        /// The size of the max location.
        /// </summary>
        private Dictionary<string, int> maxLocationSize;

        /// <summary>
        /// Writes out a Locus bed file from a list of locations
        /// </summary>
        /// <param name="locations">Locations.</param>
        /// <param name="filename">Filename.</param>
        public static void ToFileBed4Locus(List<Genomics.Location> locations, string filename)
        {
            using (TextWriter tw = Helpers.CreateStreamWriter(filename))
            {
                tw.WriteLine(string.Join("\n", locations.Select(x => string.Join("\t", new string[] { x.Chromosome, x.Start.ToString(), x.End.ToString(), x.Strand, x.Name } ))));
            }
        }

        /// <summary>
        /// Writes out a bed6 bed file from a list of locations
        /// </summary>
        /// <param name="locations">Locations.</param>
        /// <param name="filename">Filename.</param>
        public static void ToFileBed6(List<Genomics.Location> locations, string filename)
        {
            using (TextWriter tw = Helpers.CreateStreamWriter(filename))
            {
                tw.WriteLine(string.Join("\n", locations.Select(x => string.Join("\t", new string[] { x.Chromosome, x.Start.ToString(), x.End.ToString(), x.Name, x.Score.ToString(), x.Strand } ))));
            }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.BedFile"/> class.
        /// </summary>
        /// <param name='filename'>
        /// Filename.
        /// </param>
        /// <param name = "layout"></param>
        public BedFile (string filename, Layout layout, Action preParseInitCallback = null)
        {
            FileName = filename;

            FileLayout = layout;

            if (preParseInitCallback != null)
            {
                preParseInitCallback();
            }

            ParseFile();
        }

        public BedFile(Dictionary<string, Location> locations)
        {
            this.Locations = locations;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.BedFile"/> class.
        /// </summary>
        /// <param name='filename'>
        /// Filename.
        /// </param>
        /// <param name = "layout"></param>
        /// <param name = "parse"></param>
        protected BedFile (string filename, Layout layout, bool parse)
        {
            FileName = filename;

            FileLayout = layout;

            if (parse)
            {
                ParseFile();
            }
        }

        /// <summary>
        /// The bed file entry locations.
        /// </summary>
        public Dictionary<string, Genomics.Location> Locations { get; set; }

        /// <summary>
        /// Gets the locations by chromosome
        /// </summary>
        /// <value>The chromosome locations.</value>
        public Dictionary<string, Dictionary<string, Genomics.Location>> ChromosomeLocations
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.chromosomeLocations,
                    () => this.Locations
                    .ToLookup(x => x.Value.Chromosome, x => new { Name = x.Key, Value = x.Value })
                    .ToDictionary(x => x.Key, x => x.ToDictionary(y => y.Name, y => y.Value)));
            }
        }

        /// <summary>
        /// Gets the ordered locations by START position on each chromosome
        /// </summary>
        /// <value>The chromosome ordered locations.</value>
        public Dictionary<string, List<Location>> ChromosomeOrderedLocations
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.orderedLocations,
                    () => this.ChromosomeLocations.ToDictionary(x => x.Key, x => x.Value.Values.OrderBy(y => y.Start).ToList()));
            }
        }

        /// <summary>
        /// Gets the ordered locations by END position on each chromosome
        /// </summary>
        /// <value>The chromosome ordered locations.</value>
        public Dictionary<string, List<Location>> ChromosomeEndOrderedLocations
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.orderedEndLocations,
                    () => this.ChromosomeLocations.ToDictionary(x => x.Key, x => x.Value.Values.OrderBy(y => y.End).ToList()));
            }
        }

        /// <summary>
        /// Gets the chromosome locations indexed into bins by position on the chromosome
        /// </summary>
        /// <value>The chromosome indexed locations.</value>
        public Dictionary<string, Dictionary<int, List<Genomics.Location>>> ChromosomeIndexedLocations
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.indexedLocations,
                    () => this.ChromosomeLocations.ToDictionary(
                        chr => chr.Key, 
                        chr =>
                        {
                            var indexData = chr.Value.Select(loc => new 
                            {
                                Key = loc.Value.Start / IndexSize,
                                Value = loc.Value
                            });

                            return indexData.ToLookup(x => x.Key, x => x.Value).ToDictionary(x => x.Key, x => x.ToList());
                        }));
            }
        }

        /// <summary>
        /// Gets the chromosomes observed in the bed file
        /// </summary>
        /// <value>The chromosomes.</value>
        public List<string> Chromosomes 
        {
            get 
            {
                return ChromosomeLocations.Keys.ToList(); 
            }
        }

        /// <summary>
        /// Gets the ordered locations using the given field selector function.
        /// </summary>
        /// <returns>The ordered locations.</returns>
        /// <param name="selector">Selector.</param>
        public Dictionary<string, List<Genomics.Location>> GetOrderedLocations(Func<Genomics.Location, int> selector)
        {
            return this.ChromosomeLocations
                    .ToDictionary(
                        x => x.Key, 
                        x => x.Value
                            .OrderBy(y => selector(y.Value))
                            .Select(y => y.Value)
                            .ToList());
        }

        /// <summary>
        /// Gets the maximum location size. This data determines the number of index bins to search.
        /// </summary>
        /// <value>The size of the max location.</value>
        public Dictionary<string, int> MaxLocationSize
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.maxLocationSize,
                    () => this.ChromosomeLocations.ToDictionary(x => x.Key, x => x.Value.Max(loc => loc.Value.End - loc.Value.Start)));
            }
        }

        /// <summary>
        /// The names.
        /// </summary>
        public List<string> Names { get { return Locations.Keys.ToList(); } }

        /// <summary>
        /// Gets or sets the file layout.
        /// </summary>
        /// <value>The file layout.</value>
        public Layout FileLayout { get; set; }

        /// <summary>
        /// Gets or sets the name of the file.
        /// </summary>
        /// <value>The name of the file.</value>
        public string FileName { get; set; }

        /// <summary>
        /// Bed file layout structure
        /// </summary>
        public class Layout
        {
            public int Chromosome { get; set; }
            public int Start { get; set; }
            public int End { get; set; }
            public int Strand { get; set; }
            public int Score { get; set; }
            public int Name { get; set; }
            public int SignalValue { get; set; }
            public int PValue { get; set; }
            public int QValue { get; set; }
            public int Peak { get; set; }
            public int Type { get; set; }
        }

        /// <summary>
        /// The bed6 layout.
        /// </summary>
        public static Layout Bed3Layout = new Layout
        {
            Chromosome = 0,
            Start = 1,
            End = 2,
            Strand = -1,
            Score = -1,
            Name = -1,
            SignalValue = -1,
            PValue = -1,
            QValue = -1,
            Peak = -1,
            Type = -1,
        };

        /// <summary>
        /// The bed6 layout.
        /// </summary>
        public static Layout Bed6Layout = new Layout
        {
            Chromosome = 0,
            Start = 1,
            End = 2,
            Name = 3,
            Score = 4,
            Strand = 5,
            SignalValue = -1,
            PValue = -1,
            QValue = -1,
            Peak = -1,
            Type = -1,
        };

        /// <summary>
        /// The bed6 plus4 layout.
        /// </summary>
        public static Layout Bed6Plus4Layout = new Layout
        {
            Chromosome = 0,
            Start = 1,
            End = 2,
            Name = 3,
            Score = 4,
            Strand = 5,
            SignalValue = 6,
            PValue = 7,
            QValue = 8,
            Peak = 9,
            Type = -1,
        };

        /// <summary>
        /// The bed6 plus3 layout.
        /// </summary>
        public static Layout Bed6Plus3Layout = new Layout
        {
            Chromosome = 0,
            Start = 1,
            End = 2,
            Name = 3,
            Score = 4,
            Strand = 5,
            SignalValue = 6,
            PValue = 7,
            QValue = 8,
            Peak = -1,
            Type = -1,
        };

        /// <summary>
        /// The bed expression.
        /// </summary>
        public static Layout BedExpression = new Layout
        {
            Chromosome = 0,
            Start = 3,
            End = 4,
            Strand = 6,
            Name = 7,
            Score = 8,
            SignalValue = -1,
            PValue = -1,
            QValue = -1,
            Peak = -1,
            Type = -1,
        };


        /// <summary>
        /// Removes the location.
        /// </summary>
        public void RemoveLocation(string key)
        {
            this.Locations.Remove(key);
            this.Invalidate();
        }

        /// <summary>
        /// Removes the locations.
        /// </summary>
        /// <param name="keys">Keys.</param>
        public void RemoveLocations(List<string> keys)
        {
            foreach (string key in keys)
            {
                this.Locations.Remove(key);
            }

            this.Invalidate();
        }

        /// <summary>
        /// Perform parsing action on nominated file
        /// </summary>
        protected void ParseFile()
        {
            using (TextReader tr = new StreamReader(FileName))
            {
                int i = 0;

                List<Tuple<Genomics.Location, string>> data = new List<Tuple<Genomics.Location, string>>();

                string line = null;

                while ((line = tr.ReadLine()) != null)
                {
                    if (!string.IsNullOrWhiteSpace(line) && line[0] != '#')
                    {
                        ParseFields(line.Split('\t'), FileLayout, data, ref i);
                    }
                }

                var nonUnique = data.ToLookup(x => x.Item2, x => x.Item1).Where(x => x.Count() > 1).ToDictionary(x => x.Key, x => x.ToList());

                if (nonUnique.Count > 0)
                {
                    this.Locations = data.ToLookup(x => x.Item2, x => x.Item1).ToDictionary(x => x.Key, x => x.First());
                }
                else
                {
                    this.Locations = data.ToDictionary(x => x.Item2, x => x.Item1);
                }
            }
        }

        /// <summary>
        /// Parses the fields.
        /// </summary>
        /// <param name="fields">Fields.</param>
        /// <param name="layout">Layout.</param>
        /// <param name="data">Data.</param>
        /// <param name="entryCount">Entry count.</param>
        protected virtual void ParseFields(
            string[] fields, 
            Layout layout, 
            List<Tuple<Genomics.Location, string>> data, 
            ref int entryCount)
        {
            double score = 0.0;

            // Score can only be in a strictly positive location
            if (layout.Score > 0)
            {
                score = double.Parse(fields [layout.Score]);
            }


            string name = fields[layout.Chromosome] + ":" + fields[layout.Start] + "-" + fields[layout.End];
            if (layout.Name != -1)
            {
                name = this.ParseName(fields[layout.Name], ref entryCount);
            }

            data.Add(new Tuple<Genomics.Location, string>(
                new Genomics.Location
                {
                    Name = name,
                    Chromosome = fields[layout.Chromosome], 
                    Start = int.Parse(fields[layout.Start]),
                    End = int.Parse(fields[layout.End]),
                    Strand = layout.Strand != -1 ? fields[layout.Strand] : "+",
                    Score = score,
                    Data = fields,
                    AlternateName = layout.Name != -1 ? fields[layout.Name] : name
                },
                name));
        }

        /// <summary>
        /// Parses the name field
        /// </summary>
        /// <returns>The name.</returns>
        /// <param name="name">Name fields of the bed file</param>
        /// <param name="i">Number of unnamed bed entries</param>
        protected virtual string ParseName(string name, ref int i)
        {
            if (name != ".")
            {
                return name;
            }
            else
            {
                return string.Format("{0}", i++);
            }
        }

        /// <summary>
        /// Invalidate this instance.
        /// </summary>
        private void Invalidate()
        {
            this.chromosomeLocations = null;
            this.orderedLocations = null;
            this.orderedEndLocations = null;
            this.indexedLocations = null;
            this.maxLocationSize = null;
        }
    }
}

