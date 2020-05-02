//--------------------------------------------------------------------------------
// <copyright file="CorrelationMapBuilder.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Data;
    using Genomics;
    using Shared;
    using Tools;

    /// <summary>
    /// Correlation map builder.
    /// </summary>
    public class CorrelationMapBuilder : MapBuilder<CorrelationMapBuilder.HistoneExpressionLinkData, CorrelationMapBuilder.LinkCorrelation>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.CorrelationMapBuilder"/> class.
        /// </summary>
        /// <param name="xmlFile">Xml configuration file.</param>
        /// <param name = "omittedTissues">CSV list of cell lines to omit from mapping.</param>
        public CorrelationMapBuilder(string xmlFile, string omittedTissues)
            : base(
                xmlFile,
                omittedTissues != null ? omittedTissues.Split(',').ToArray() : new string[0])
        {
        }

        /// <summary>
        /// Gets the name of the mapping method name
        /// </summary>
        /// <value>The name of the method.</value>
        public override string MethodName
        {
            get
            {
                return "Correlation";
            }
        }

        /// <summary>
        /// Queries the interface.
        /// </summary>
        /// <returns>The interface.</returns>
        /// <param name="t">Type to query.</param>
        protected override object QueryInterface(Type t)
        {
            return t == typeof(IMapBuilder) ? new CorrelationMapBuilderProxy(this) as IMapBuilder : base.QueryInterface(t);
        }

        /// <summary>
        /// Construct the data used to establish a confidence score for the link
        /// </summary>
        /// <returns>The link data.</returns>
        /// <param name="tss">Transcription start site name.</param>
        /// <param name="locus">Locus name.</param>
        /// <param name="expressionByTissue">Expression by cell line.</param>
        /// <param name="histoneByTissue">Histone by cell line.</param>
        protected override HistoneExpressionLinkData BuildLinkData(
            string tss, 
            string locus, 
            Dictionary<string, double> expressionByTissue, 
            Dictionary<string, double> histoneByTissue)
        {
            if (histoneByTissue.Count <= 3)
            {
                return null;
            }

            Stats.Point2[] vector = new Stats.Point2[histoneByTissue.Count];
            int i = 0;
            foreach (var data in histoneByTissue)
            {
                vector[i].x = Math.Log10(data.Value + 1);
                vector[i].y = Math.Log10(expressionByTissue[data.Key] + 1);
                i++;
            }

            if (vector.Length >= this.MinFeatureCount && 
                vector.Count(point => point.x > 0) > 0 && 
                vector.Count(point => point.y >= Math.Log10(this.MinimaxExpressionValue + 1)) > 0 && 
                IsTwoFoldChange(vector.Select(point => point.y)))
            {
                double correlation = Stats.Pearson(vector);
                double zscore = Stats.FisherTransformZScore(correlation, vector.Length);
                return new HistoneExpressionLinkData
                {
                    Tss = tss,
                    Locus = locus,
                    Correlation = correlation,
                    ZScore = zscore
                    ////Vector = vector
                };
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// Determines whether specified link data is valid
        /// </summary>
        /// <returns>true if the link data meets mapping requirements.</returns>
        /// <param name="data">Link data.</param>
        protected override bool IsValidLinkData(HistoneExpressionLinkData data)
        {
            return data != null && !double.IsNaN(data.ZScore);
            ////return data.Vector.Count() >= this.MinFeatureCount && 
            ////    data.Vector.Count(point => point.x > 0) > 0 && 
            ////    data.Vector.Count(point => point.y >= Math.Log10(this.MinimaxExpressionValue + 1)) > 0 && 
            ////    IsTwoFoldChange(data.Vector.Select(point => point.y));
        }

        /// <summary>
        /// Determines whether the specified expressionData is valid
        /// </summary>
        /// <returns>true if the expression data fits the mapping requirements.</returns>
        /// <param name="expressionData">Expression data.</param>
        protected override bool IsValidExpressionData(MapBuilderData.TissueExpressionData expressionData)
        {
            return expressionData.Expression
                .Count >= this.MinFeatureCount && expressionData.Expression.Count(y => y.Value >= this.MinimaxExpressionValue) > 0;
        }

        /// <summary>
        /// Gets the link score from the given data.
        /// </summary>
        /// <returns>The link scores.</returns>
        /// <param name="corrVectors">Link data.</param>
        /// <param name="count">Total number of links.</param>
        protected override LinkCorrelation[] GetCorrelations(IEnumerable<HistoneExpressionLinkData> corrVectors, int count)
        {
            LinkCorrelation[] data = new LinkCorrelation[count];

            var zscores = corrVectors.Select(x => x.ZScore);
            double[] logpvalues = Calculations.ZscoreToLogPvalue(zscores);

            int i = 0;
            foreach (var corr in corrVectors)
            {
                var corrData = (HistoneExpressionLinkData)corr;
                double lp = logpvalues[i];
                data[i] = new LinkCorrelation
                {
                    Tss = corrData.Tss,
                    Locus = corrData.Locus,
                    Correlation = corrData.Correlation,
                    Pvalue = Math.Exp(lp),
                    LogPvalue = lp
                };
                i++;
            }

            return data;
        }

        /// <summary>
        /// Formats the score into a bed file line.
        /// </summary>
        /// <returns>The score bed file line.</returns>
        /// <param name="scoreData">Score data.</param>
        /// <param name="location">Location of the target gene.</param>
        /// <param name="distance">Locus-TSS distance.</param>
        protected override string FormatScoreTsvLine(LinkCorrelation scoreData, Location location, int distance)
        {
            return string.Format(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:F15}\t{8:" + 
	          (scoreData.Pvalue == 0 ? "" : "0.000e+000") + "}\t{9}\t{10}\t{11}",
                location.Chromosome,
                location.DirectionalStart,
                location.DirectionalStart,
                scoreData.Tss,
                "NA", ////scoreData.CorrectedPvalue,
                location.Strand,
                scoreData.Locus,
                scoreData.Correlation,
                //scoreData.Pvalue == 0 ? 2e-322 : scoreData.Pvalue,
                scoreData.Pvalue,
                distance,
                this.HistoneName,
                location.AlternateName);
        }

        /// <summary>
        /// Determines if is two fold change the specified data.
        /// </summary>
        /// <returns><c>true</c> if is two fold change the specified data; otherwise, <c>false</c>.</returns>
        /// <param name="data">Expression data across cell lines.</param>
        protected static bool IsTwoFoldChange(IEnumerable<double> data)
        {
            var dataArray = data.ToArray();
            double mean = dataArray.Average();
            double min = dataArray.Min();
            double max = dataArray.Max();

            if (Math.Abs(mean) < double.Epsilon)
            {
                return false;
            }

            if (max - mean >= 2 * mean)
            {
                return true;
            }

            if (Math.Abs(mean - min) < double.Epsilon)
            {
                return false;
            }

            return mean - min >= mean / 3;
        }

        /// <summary>
        /// Histone expression link data.
        /// </summary>
        public class HistoneExpressionLinkData : MapBuilderData.BaseLinkScore
        {
            /// <summary>
            /// Gets or sets the correlation.
            /// </summary>
            /// <value>The correlation.</value>
            public double Correlation { get; set; }

            /// <summary>
            /// Gets or sets the Z score.
            /// </summary>
            /// <value>The Z score.</value>
            public double ZScore { get; set; }

            ////public IEnumerable<Stats.Point2> Vector { get; set; }
        }

        /// <summary>
        /// Link correlation.
        /// </summary>
        public class LinkCorrelation : MapBuilderData.BaseLinkScore
        {
            /// <summary>
            /// Gets or sets the correlation.
            /// </summary>
            /// <value>The correlation.</value>
            public double Correlation { get; set; }

            /// <summary>
            /// Gets or sets the log pvalue.
            /// </summary>
            /// <value>The log pvalue.</value>
            public double LogPvalue { get; set; }

            /// <summary>
            /// Gets or sets the pvalue.
            /// </summary>
            /// <value>The pvalue.</value>
            public double Pvalue { get; set; }

            ////public double CorrectedPvalue { get; set; }
        }

        /// <summary>
        /// Executes the correlation map builder task.
        /// </summary>
        public class Executor : IAnalysisExecutor<CorrelationMapBuilder, CorrelationMapBuilder.Executor.Arguments>
        {
            /// <summary>
            /// Command line arguments.
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The cell line to be mapped.
                /// </summary>
                Tissue,

                /// <summary>
                /// XML configuration file name.
                /// </summary>
                Config,

                /// <summary>
                /// The cell lines omitted from mapping.
                /// </summary>
                OmittedTissues,

                /// <summary>
                /// The name of the histone to use for mapping.
                /// </summary>
                HistoneName,

                /// <summary>
                /// The RNA source type to use for mapping.
                /// </summary>
                RnaSource,

                /// <summary>
                /// The stage in the mapping pipeline.
                /// </summary>
                Stage,

                /// <summary>
                /// The locus filter range.
                /// </summary>
                LocusFilterRange,

                /// <summary>
                /// The name of the Locus file.
                /// </summary>
                LocusFileName,

                /// <summary>
                /// The use genes.
                /// </summary>
                UseGenes,

                /// <summary>
                /// The name of the output file.
                /// </summary>
                MapFileName,

                /// <summary>
                /// The name of the output directory.
                /// </summary>
                OutDir,
            }

            /// <summary>
            /// Stage of the map building process.
            /// </summary>
            private enum Stage
            {
                /// <summary>
                /// Preprocess the expression and histone data only.
                /// </summary>
                Preprocess,

                /// <summary>
                /// Build the map only using the pre-processed data.
                /// </summary>
                BuildMap,

                /// <summary>
                /// Run the full pipeline in one execution.
                /// </summary>
                RunPipeline,

            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Map builder that calculates a confidence score that a locus regulates a TSS within the specified distance range.\n" +
                        "Confidence is based on the cross-tissue Pearson Correlation Coefficient between histone marks at loci and expression at TSSes";
                }
            }

            /// <summary>
            /// Gets the command line option descriptions.
            /// </summary>
            /// <value>The options data.</value>
            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { Arguments.Tissue, "Cell line whose loci will be mapped; can be \"None\"" },
                        { Arguments.HistoneName, "Histone mark to use for mapping" },
                        { Arguments.Config, "XML configuration file name" },
                        { Arguments.OmittedTissues, "CSV of cell lines whose expression and histone data should be omitted from the mapping process" },
                        { Arguments.RnaSource, "RNA-seq source whose TSSes will be used as putative locus targets" },
                        { Arguments.Stage, "Stage in the mapping process: {" + string.Join(", ", Enum.GetNames(typeof(Stage))) + "}" },
                        { Arguments.LocusFilterRange, "Optional parameter to remove loci within the specified distance of the transcript body (TSS to TTS)" },
                        { Arguments.LocusFileName, "Optional manual override of loci filename" },
                        { Arguments.UseGenes, "Optional flag to map genes instead of transcripts" },
                        { Arguments.MapFileName, "Optional parameter to specify file to contain the map" },
                        { Arguments.OutDir, "Optional parameter to specify the output directory" },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis for which this class is the factor with the given command line arguments.
            /// </summary>
            /// <param name="commandArgs">Command arguments.</param>
            public override void Execute(Args commandArgs)
            {
                IMapBuilder builder = IUnknown.QueryInterface<IMapBuilder>(new CorrelationMapBuilder(
                    commandArgs.StringEnumArgs[Arguments.Config],
                    commandArgs.StringEnumArgs.GetOptional(Arguments.OmittedTissues)));

                this.ReflectStringArgs(builder, new Arguments[]
                {
                    Arguments.Tissue,
                    Arguments.HistoneName,
                    Arguments.RnaSource,
                });

                this.ReflectOptionalStringArgs(builder, new Arguments[]
                {
                    Arguments.LocusFileName,
                    Arguments.MapFileName,
                });

                this.ReflectOptionalString(builder, Arguments.OutDir, ".");

                this.ReflectFlag(builder, Arguments.UseGenes);
                this.ReflectOptionalInt(builder, Arguments.LocusFilterRange, 0);

                Stage stage = (Stage)Enum.Parse(typeof(Stage), commandArgs.StringEnumArgs[Arguments.Stage]);
                switch (stage)
                {
                    case Stage.Preprocess:
                        builder.PreprocessExpression();
                        break;

                    case Stage.BuildMap:
                        builder.BuildMap();
                        break;

                    case Stage.RunPipeline:
                        builder.PreprocessExpression();
                        builder.BuildMap();
                        break;

                    default: 
                        throw new Exception("Invalid stage argument: " + commandArgs.StringArgs["Stage"] + "\nValid values: Preprocess, BuildMap");
                }
            }
        }

        /// <summary>
        /// Correlation map builder proxy.
        /// </summary>
        private class CorrelationMapBuilderProxy : IMapBuilder
        {
            /// <summary>
            /// The data.
            /// </summary>
            private CorrelationMapBuilder instance;

            /// <summary>
            /// Initializes a new instance of the <see cref="CorrelationMapBuilderProxy"/> class.
            /// </summary>
            /// <param name="instance">Instance to wrap in a proxy.</param>
            public CorrelationMapBuilderProxy(CorrelationMapBuilder instance)
            {
                this.instance = instance;
            }

            public string LocusFileName
            {
                get
                {
                    return this.instance.LocusFileName;
                }

                set
                {
                    this.instance.LocusFileName = value;
                }
            }

            public string MapFileName { get { return this.instance.MapFileName; } set { this.instance.MapFileName = value; } }

            public string OutDir { get { return this.instance.OutDir; } set { this.instance.OutDir = value; } }


            /// <summary>
            /// Gets or sets the cell line.
            /// </summary>
            /// <value>The cell line.</value>
            public string Tissue 
            {
                get
                {
                    return this.instance.Tissue;
                }

                set
                {
                    this.instance.Tissue = value;
                }
            }

            public string[] TissueSources 
            {
                get
                {
                    return this.instance.TissueSources;
                }

                set
                {
                    this.instance.TissueSources = value;
                }
            }

            /// <summary>
            /// Gets or sets the name of the histone.
            /// </summary>
            /// <value>The name of the histone.</value>
            public string HistoneName
            {
                get
                {
                    return this.instance.HistoneName;
                }

                set
                {
                    this.instance.HistoneName = value;
                }
            }

            /// <summary>
            /// Gets or sets the rna source.
            /// </summary>
            /// <value>The rna source.</value>
            public string RnaSource
            {
                get
                {
                    return this.instance.RnaSource;
                }

                set
                {
                    this.instance.RnaSource = value;
                }
            }

            /// <summary>
            /// Gets or sets a value indicating whether this
            /// <see cref="Analyses.CorrelationMapBuilder+CorrelationMapBuilderProxy"/> use genes.
            /// </summary>
            /// <value><c>true</c> if use genes; otherwise, <c>false</c>.</value>
            public bool UseGenes
            {
                get
                {
                    return this.instance.UseGenes;
                }

                set 
                {
                    this.instance.UseGenes = value;
                }
            }

            /// <summary>
            /// Gets or sets a value indicating the distance to remove loci withing a minimum distance of an annotated transcript
            /// </summary>
            /// <value>Range to remove around the promoter region</value>
            public int LocusFilterRange
            {
                get
                {
                    return this.instance.LocusFilterRange;
                }

                set
                {
                    this.instance.LocusFilterRange = value;
                }
            }

            /// <summary>
            /// Preprocesses the expression.
            /// </summary>
            public void PreprocessExpression()
            {
                this.instance.PreprocessExpression();
            }

            /// <summary>
            /// Builds the map.
            /// </summary>
            public void BuildMap()
            {
                this.instance.BuildMap();
            }
        }
    }
}
