//--------------------------------------------------------------------------------
// File: IExpressionData.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Genomics;
    using Shared;
    using Tools;
    using Data;

    /// <summary>
    /// Expression data set representation
    /// </summary>
    abstract public class IExpressionData
    {
        /// <summary>
        /// Factory method for implementations of this interface
        /// </summary>
        /// <returns>The expression data.</returns>
        /// <param name="filename">Filename.</param>
        /// <param name="rnaType">Rna type.</param>
        /// <param name="filetype">Filetype.</param>
        /// <param name="annotation">Annotation.</param>
        public static IExpressionData LoadExpressionData(string filename, string rnaType, string filetype, IAnnotation annotation)
        {
            Console.WriteLine("\tLoading expression file " + filename + " ... ");
            IExpressionData data = null;
            if (filetype == "gtf" || filetype == ".gtf")
            {
                data = IUnknown.QueryInterface<IExpressionData>(new GtfExpressionFile(GtfExpressionFile.ExpressionTypeFromString(rnaType), filename, annotation));
            } 
            else
            {
                data = IUnknown.QueryInterface<IExpressionData>(new BedExpressionFile(filename, BedFile.BedExpression, annotation));
            }

            Console.WriteLine("\ttranscript count: " + data.Transcripts.Count);
            Console.WriteLine("\tgene count: " + data.Genes.Count);

            return data;
        }

        //- Abstract Methods -------------------------------------------------------------------------------

        /// <summary>
        /// All transcripts in the expression data set.
        /// </summary>
        /// <value>The transcripts.</value>
        abstract public Dictionary<string, Location> Transcripts { get; }

        //- Base Implementation -------------------------------------------------------------------------------

        /// <summary>
        /// Cache of expression data
        /// </summary>
        protected Dictionary<string, double> transcriptExpression;

        /// <summary>
        /// Cache of expression data zscores
        /// </summary>
        protected Dictionary<string, double> transcriptExpressionZScore;

        /// <summary>
        /// Cache of expression data by gene
        /// </summary>
        protected Dictionary<string, double> geneExpression;
        protected Dictionary<string, double> log2GeneExpression;

        /// <summary>
        /// Cache of expression zscore data by gene
        /// </summary>
        protected Dictionary<string, double> geneExpressionZScore;
        protected Dictionary<string, double> log2GeneExpressionZScore;

        /// <summary>
        /// The genes.
        /// </summary>
        protected Dictionary<string, Location> genes;

        /// <summary>
        /// The gene expression mean.
        /// </summary>
        private double geneExpressionMean = double.NaN;
        private double log2GeneExpressionMean = double.NaN;

        /// <summary>
        /// The gene expression std dev.
        /// </summary>
        private double geneExpressionStdDev = double.NaN;
        private double log2GeneExpressionStdDev = double.NaN;

        /// <summary>
        /// Genes the expression from transcripts.
        /// </summary>
        /// <returns>The expression from transcripts.</returns>
        /// <param name="transcriptSubset">Transcript subset.</param>
        public Dictionary<string, double> GeneExpressionFromTranscripts(HashSet<string> transcriptSubset)
        {
            var geneSet = new HashSet<string>(this.Transcripts
                .Where(x => transcriptSubset.Contains(x.Value.Name))
                .Select(x => x.Value.AlternateName));

            return this.GeneExpression
                    .Where(x => geneSet.Contains(x.Key))
                    .ToDictionary(x => x.Key, x => x.Value);
        }

        /// <summary>
        /// Gets the transcript expression.
        /// </summary>
        /// <value>The transcript expression.</value>
        public Dictionary<string, double> TranscriptExpression
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcriptExpression,
                    () => this.transcriptExpression = this.Transcripts.ToDictionary(x => x.Key, x => x.Value.Score));
            }
        }

        /// <summary>
        /// Gets the transcript expression Z score.
        /// </summary>
        /// <value>The transcript expression Z score.</value>
        public Dictionary<string, double> TranscriptExpressionZScore
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcriptExpressionZScore,
                    () => 
                {
                    var mean = this.TranscriptExpression.Values.Average();
                    var sd = Stats.StdDev(this.TranscriptExpression.Values);

                    return this.TranscriptExpression.ToDictionary(x => x.Key, x => (x.Value - mean) / sd);
                });
            }
        }

        /// <summary>
        /// Gets the gene expression.
        /// </summary>
        /// <value>The gene expression.</value>
        public Dictionary<string, double> GeneExpression
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.geneExpression,
                    () => this.geneExpression = this.Transcripts
                        .ToLookup(x => x.Value.AlternateName, x => x.Value)
                        .ToDictionary(x => x.Key, x => x.Sum(y => y.Score)));
            }
        }

        /// <summary>
        /// Gets the log2 gene expression.
        /// </summary>
        /// <value>The log2 gene expression.</value>
        public Dictionary<string, double> Log2GeneExpression
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.log2GeneExpression,
                    () => this.log2GeneExpression = this.Transcripts
                        .ToLookup(x => x.Value.AlternateName, x => x.Value)
                        .ToDictionary(
                            x => x.Key, 
                        x => Math.Log(x.Sum(y => y.Score + 1)) / Math.Log(2)
                    ));
            }
        }

        /// <summary>
        /// Gets the gene expression mean.
        /// </summary>
        /// <value>The gene expression mean.</value>
        public double GeneExpressionMean
        {
            get
            {
                if (double.IsNaN(this.geneExpressionMean))
                {
                    this.geneExpressionMean = this.GeneExpression.Values.Average();
                }

                return this.geneExpressionMean;
            }
        }

        /// <summary>
        /// Gets the log2 gene expression mean.
        /// </summary>
        /// <value>The log2 gene expression mean.</value>
        public double Log2GeneExpressionMean
        {
            get
            {
                if (double.IsNaN(this.log2GeneExpressionMean))
                {
                    this.log2GeneExpressionMean = this.Log2GeneExpression.Values.Average();
                }

                return this.log2GeneExpressionMean;
            }
        }

        /// <summary>
        /// Gets the gene expression std dev.
        /// </summary>
        /// <value>The gene expression std dev.</value>
        public double GeneExpressionStdDev
        {
            get
            {
                if (double.IsNaN(this.geneExpressionStdDev))
                {
                    this.geneExpressionStdDev = Stats.StdDev(this.GeneExpression.Values);
                }

                return this.geneExpressionStdDev;
            }
        }

        /// <summary>
        /// Gets the log2 gene expression std dev.
        /// </summary>
        /// <value>The log2 gene expression std dev.</value>
        public double Log2GeneExpressionStdDev
        {
            get
            {
                if (double.IsNaN(this.log2GeneExpressionStdDev))
                {
                    this.log2GeneExpressionStdDev = Stats.StdDev(this.Log2GeneExpression.Values);
                }

                return this.log2GeneExpressionStdDev;
            }
        }

        /// <summary>
        /// Gets the gene expression Z score.
        /// </summary>
        /// <value>The gene expression Z score.</value>
        public Dictionary<string, double> GeneExpressionZScore
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.geneExpressionZScore,
                    () => this.GeneExpression
                        .ToDictionary(x => x.Key, x => (x.Value - this.GeneExpressionMean) / this.GeneExpressionStdDev));
            }
        }

        /// <summary>
        /// Gets the log2 gene expression Z score.
        /// </summary>
        /// <value>The log2 gene expression Z score.</value>
        public Dictionary<string, double> Log2GeneExpressionZScore
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.log2GeneExpressionZScore,
                    () => this.Log2GeneExpression
                        .ToDictionary(x => x.Key, x => (x.Value - this.Log2GeneExpressionMean) / this.Log2GeneExpressionStdDev));
            }
        }

        /// <summary>
        /// Gets the genes.
        /// </summary>
        /// <value>The genes.</value>
        public Dictionary<string, Genomics.Location> Genes
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.genes,
                    () => this.Transcripts
                        .ToLookup(x => x.Value.AlternateName, x => x.Value)
                        .ToDictionary(x => x.Key, x => x.First()));
            }
        }
    }
}

