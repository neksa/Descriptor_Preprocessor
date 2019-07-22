//--------------------------------------------------------------------------------
// File: CuffDiffFile.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using Shared;
    using System.Text.RegularExpressions;

    /// <summary>
    /// Class for encapsulating data coming from a cuffdiff run
    /// </summary>
    public class CuffDiffData
    {
        public string TestId { get; set; }
        public string GeneId { get; set; }
        public string Gene { get; set; }
        public Genomics.Location Locus { get; set; }
        public string Status { get; set; }
        public double Fpkm1 { get; set; }
        public double Fpkm2 { get; set; }
        public double FoldChange { get; set; }
        public double PValue { get; set; }
        public double QValue { get; set; }
        public bool Significant { get; set; }
    }

    /// <summary>
    /// Object representation of a cuff diff output file
    /// </summary>
    public class CuffDiffFile
    {
        /// <summary>
        /// The filename.
        /// </summary>
        private readonly string filename;

        /// <summary>
        /// The transcript data.
        /// </summary>
        private Dictionary<string, CuffDiffData> transcriptData;

        /// <summary>
        /// Significantly changed genes
        /// </summary>
        private List<string> significantGenes;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.CuffDiffData"/> class.
        /// </summary>
        public CuffDiffFile(string filename)
        {
            this.filename = filename;
        }

        /// <summary>
        /// Significantly changed genes
        /// </summary>
        public List<string> SignificantlyChangedGenes
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.significantGenes,
                    () => this.TranscriptData.Values.Where(x => x.Significant).Select(x => x.Gene).ToList());
            }
        }

        /// <summary>
        /// Gets the transcript data from the differential expression experiment
        /// </summary>
        /// <value>The transcript data.</value>
        public Dictionary<string, CuffDiffData> TranscriptData
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcriptData,
                    () =>
                    {
                        using (TextReader tr = new StreamReader(this.filename))
                        {
                            var headers = tr.ReadLine().Trim().Split('\t')
                                .Select((x, i) => new { Field = x, Index = i })
                                .ToDictionary(x => x.Field, x => x.Index);

                            return tr.ReadToEnd().Split('\n').Where(line => !string.IsNullOrWhiteSpace(line))
                                .Select(line => 
                                {
                                    var fields = line.Split('\t');

                                    if (fields.Length != headers.Count)
                                    {
                                        throw new Exception(string.Format("Invalid number of fields ({0}) in line:\n\t{1}", fields.Length, line));
                                    }

                                    var locusFields = fields[headers["locus"]].Split(new char[] { ':', '-' } );

                                    Func<string, double> parseDouble = (string arg) => 
                                    {
                                        if (arg == "inf")
                                        {
                                            return double.PositiveInfinity;
                                        }
                                        else if (arg == "-inf")
                                        {
                                            return double.NegativeInfinity;
                                        }
                                        else
                                        {
                                            return double.Parse(arg);
                                        }
                                    };

                                    Regex ensembleFormat = new Regex("^ENS[TG]");

                                    Func<string, string> RemoveEnsemblSuffix = (id) => 
                                    {
                                        if (ensembleFormat.IsMatch(id))
                                        {
                                            return id.Split('.')[0];
                                        }

                                        return id;
                                    };

                                    return new CuffDiffData
                                    {
                                        TestId = RemoveEnsemblSuffix(fields[headers["test_id"]]),
                                        GeneId = RemoveEnsemblSuffix(fields[headers["gene_id"]]),
                                        Gene = fields[headers["gene"]],
                                        Locus = new Genomics.Location
                                        {
                                            Name = RemoveEnsemblSuffix(fields[headers["gene_id"]]),
                                            Chromosome = locusFields[0],
                                            Start = int.Parse(locusFields[1]),
                                            End = int.Parse(locusFields[2]),
                                        },
                                        Status = fields[headers["status"]],
                                        Fpkm1 = parseDouble(fields[headers["value_1"]]),
                                        Fpkm2 = parseDouble(fields[headers["value_2"]]),
                                        FoldChange = parseDouble(fields[headers["log2(fold_change)"]]),
                                        PValue = parseDouble(fields[headers["p_value"]]),
                                        QValue = parseDouble(fields[headers["q_value"]]),
                                        Significant = fields[headers["significant"]] == "yes"
                                    };
                                })
                                .ToLookup(x => x.GeneId, x => x)
                                .ToDictionary(x => x.Key, x => x.OrderBy(y => Math.Abs(y.FoldChange)).Last());
                        }
                    });
            }
        }
    }
}

