
namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Data;
    using Genomics;
    using Shared;
    using Tools;

    public class CorrelationMapEligibleGenes : CorrelationMapBuilder
    {
        public CorrelationMapEligibleGenes(string xmlFile, string omittedTissues)
            : base(xmlFile, omittedTissues)
        {
        }

        public string AnnotationFileName { get; set; }

        public void Execute()
        {
            var validTss = this.TranscriptExpression
                .Select(x => new
                {
                    Values = x.Value.Values,
                    TissueExpressionData = new Analyses.MapBuilderData.TissueExpressionData
                    {   
                        Tss = x.Key,
                        Location = this.TranscriptLocations[x.Key],
                        Expression = x.Value,
                    }
                })
                .Where(x => this.IsValidExpressionData(x.TissueExpressionData) && CorrelationMapBuilder.IsTwoFoldChange(x.Values))
                .ToList();

            string stem = string.Join("_", new string[]
            {
                "TssSet",
                this.Tissue,
                this.RnaSource,
                this.HistoneName,
            }
                .Concat(this.OmittedTissues != null ? this.OmittedTissues : new string[] {}));

            string file = string.Format("../temp/results/TssSets/{0}.nsv", stem);

            Tables.ToNamedNsvFile(file, validTss.Select(x => x.TissueExpressionData.Tss));

            if (this.AnnotationFileName != null)
            {
                var annotation = new GtfExpressionFile(
                                     GtfExpressionFile.ExpressionTypeFromString(this.RnaSource), 
                                     this.AnnotationFileName);

                var tssToGene = annotation.Transcripts.ToDictionary(x => x.Key, x => x.Value.AlternateName);

                string geneStem = string.Join("_", new string[]
                {
                    "GeneSet",
                    this.Tissue,
                    this.RnaSource,
                    this.HistoneName,
                }
                    .Concat(this.OmittedTissues != null ? this.OmittedTissues : new string[] { }));

                string geneFile = string.Format("../temp/results/TssSets/{0}.nsv", geneStem);

                Tables.ToNamedNsvFile(geneFile, validTss
                    .Where(x => tssToGene.ContainsKey(x.TissueExpressionData.Tss))
                    .ToLookup(x => tssToGene[x.TissueExpressionData.Tss], x => x)
                    .Select(x => x.Key));
            }


        }

        public new class Executor : IAnalysisExecutor<CorrelationMapEligibleGenes, CorrelationMapEligibleGenes.Executor.Arguments>
        {
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
                /// The use genes.
                /// </summary>
                UseGenes,

                /// <summary>
                /// The name of the annotation file.
                /// </summary>
                AnnotationFileName,
            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Simple filter function to get the CorrelationMapBuilder elligible TSS set";
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
                        { Arguments.Tissue, "Cell line whose Loci will be mapped; can be \"None\"" },
                        { Arguments.HistoneName, "Histone mark to use for mapping" },
                        { Arguments.Config, "XML configuration file name" },
                        { Arguments.OmittedTissues, "CSV of cell lines whose expression and histone data should be omitted from the mapping process" },
                        { Arguments.RnaSource, "RNA-seq source whose TSSes will be used as putative Locus targets" },
                        { Arguments.AnnotationFileName, "Optional annotation file for converting TSSs to Genes" },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis for which this class is the factor with the given command line arguments.
            /// </summary>
            /// <param name="commandArgs">Command arguments.</param>
            public override void Execute(Args commandArgs)
            {
                Console.WriteLine(commandArgs.StringEnumArgs[Arguments.Tissue]);

                var builder = new CorrelationMapEligibleGenes(
                    commandArgs.StringEnumArgs[Arguments.Config],
                    commandArgs.StringEnumArgs.GetOptional(Arguments.OmittedTissues));

                builder.Tissue = commandArgs.StringEnumArgs[Arguments.Tissue];
                builder.HistoneName = commandArgs.StringEnumArgs[Arguments.HistoneName];
                builder.RnaSource = commandArgs.StringEnumArgs[Arguments.RnaSource];
                builder.UseGenes = commandArgs.Flags.Contains(Arguments.UseGenes.ToString());

                builder.AnnotationFileName = this.GetOptionalStringArg(commandArgs, Arguments.AnnotationFileName);

                builder.Execute();
            }
        }
    }
}

