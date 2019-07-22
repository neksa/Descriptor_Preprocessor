//--------------------------------------------------------------------------------
// <copyright file="Tables.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Data
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.IO;
    using System.Text.RegularExpressions;
    using Shared;
    using Genomics;

    /// <summary>
    /// Manages creation of data tables suitable for R.
    /// </summary>
    public static class Tables
    {

        /// <summary>
        /// Dumps data to column file.
        /// </summary>
        /// <returns>The to column file.</returns>
        /// <param name="data">Data.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public static string ToNsvFile<T>(IEnumerable<T> data)
        {
            string dataFile = Path.GetTempFileName();
            return ToNamedNsvFile(dataFile, data);
        }

        public static string ToNamedNsvFile<T>(string fileName, IEnumerable<T> data)
        {
            using (TextWriter tw = Helpers.CreateStreamWriter(fileName))
            {
                tw.WriteLine(string.Join("\n", data));
            }

            return fileName;
        }



        /// <summary>
        /// Dumps data to a table file.
        /// </summary>
        /// <returns>The to column file.</returns>
        /// <param name="data">Data.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <param name = "labels"></param>
        public static string ToTsvFile<T>(IEnumerable<IEnumerable<T>> data, IEnumerable<string> labels)
        {
            string dataFile = Path.GetTempFileName();
            //var d = data.ToArray();
            using (TextWriter tw = Helpers.CreateStreamWriter(dataFile))
            {
                if (labels != null)
                {
                    tw.WriteLine(string.Join("\t", labels));
                }

                foreach (var line in data.Select(x => string.Join("\t", x)))
                {
                    tw.WriteLine(line);
                }
            }

            return dataFile;
        }

        /// <summary>
        /// Dumps data to a table file at the specified file name.
        /// </summary>
        /// <returns>The to column file.</returns>
        /// <param name = "filename"></param>
        /// <param name="data">Data.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        /// <param name = "labels"></param>
        public static string ToNamedTsvFile<T>(string filename, IEnumerable<IEnumerable<T>> data, IEnumerable<string> labels)
        {
            string dataFile = filename;
            //var d = data.ToArray();
            using (TextWriter tw = Helpers.CreateStreamWriter(dataFile))
            {
                if (labels != null)
                {
                    tw.WriteLine(string.Join("\t", labels));
                }

                foreach (var line in data.Select(x => string.Join("\t", x)))
                {
                    tw.WriteLine(line);
                }
            }

            return dataFile;
        }

        /// <summary>
        /// Appends to named tsv file.
        /// </summary>
        /// <returns>The to named tsv file.</returns>
        /// <param name="filename">Filename.</param>
        /// <param name="data">Data.</param>
        /// <param name="labels">Labels.</param>
        public static string AppendToNamedTsvFile<T>(string filename, IEnumerable<IEnumerable<T>> data, IEnumerable<string> labels)
        {
            string dataFile = filename;
            //var d = data.ToArray();
            using (TextWriter tw = Helpers.AppendStreamWriter(dataFile))
            {
                if (labels != null)
                {
                    tw.WriteLine(string.Join("\t", labels));
                }

                foreach (var line in data.Select(x => string.Join("\t", x)))
                {
                    tw.WriteLine(line);
                }
            }

            return dataFile;
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
        /*public static int CreateTssTable(
            Dictionary<string, List<string>> LocusMap, 
            IndexedScores tssScores, 
            Dictionary<string, double> expression,
            string outFileName, 
            List<string> chipSets)
        {
            var rawScores = LocusMap
                .Where(mapEntry => expression.ContainsKey(mapEntry.Key) && tssScores.Index.ContainsKey(mapEntry.Key))
                .Select(mapEntry => 
                {
                    var atTss = tssScores.Scores.Select(chipData => chipData[tssScores.Index[mapEntry.Key]]).ToList();

                    return new
                    {
                        Tss = mapEntry.Key,
                        PromoterScores = atTss,
                        Exp = expression[mapEntry.Key]
                    };
                })
                .ToList();

            var scores = rawScores.Select(x => 
            {
                return new 
                {
                    x.Tss,
                    x.PromoterScores,
                    x.Exp
                };
            })
                .ToList();

            Console.WriteLine("Number of TSSes: {0}", scores.Count);

            using (TextWriter tw = Helpers.CreateStreamWriter(outFileName))
            {
                // table format
                tw.WriteLine("{0}\tesc", 
                    string.Join("\t", chipSets.Select(x => string.Format("Tss_{0}", x.Trim().Split('/').Last().Split('.').First().Split('-')[0]))), 
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
        }*/




    }
}

