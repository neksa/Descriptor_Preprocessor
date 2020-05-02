//--------------------------------------------------------------------------------
// File: NullMapBuilder.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Shared;
    using Genomics;
    using System.IO;

    /// <summary>
	/// Constructs a map of all possible Locus-TSS links
	/// </summary>
    public class NullMapBuilder : MapBuilder<MapBuilderData.BaseLinkScore, MapBuilderData.BaseLinkScore>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.NullMapBuilder"/> class.
        /// </summary>
        /// <param name="xmlFile">Xml file.</param>
		public NullMapBuilder(string xmlFile)
			: base(xmlFile)
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
				return "Null";
			}
		}

		/// <summary>
		/// Queries the interface.
		/// </summary>
		/// <returns>The interface.</returns>
		/// <param name="t">T.</param>
		override protected object QueryInterface(Type t)
		{
			return t == typeof(IMapBuilder) ? new NullMapBuilderProxy(this) as IMapBuilder : base.QueryInterface(t);
		}

		/// <summary>
		/// Construct the data used to establish a confidence score for the link
		/// </summary>
		/// <returns>The link data.</returns>
		/// <param name="tss">Tss.</param>
		/// <param name="Locus">Locus.</param>
		/// <param name="expressionByTissue">Expression by cell line.</param>
		/// <param name="histoneByTissue">Histone by cell line.</param>
        protected override MapBuilderData.BaseLinkScore BuildLinkData(string tss, string Locus, Dictionary<string, double> expressionByTissue, Dictionary<string, double> histoneByTissue)
		{
            return new MapBuilderData.BaseLinkScore { Tss = tss, Locus = Locus, Score = histoneByTissue.First().Value };
		}

		/// <summary>
		/// Determines whether specified link data is valid
		/// </summary>
		/// <returns>true</returns>
		/// <c>false</c>
		/// <param name="data">Data.</param>
        protected override bool IsValidLinkData(MapBuilderData.BaseLinkScore data)
		{
			return true;
		}

		/// <summary>
		/// Determines whether the specified expressionData is valid
		/// </summary>
		/// <returns>true</returns>
		/// <c>false</c>
		/// <param name="expressionData">Expression data.</param>
        protected override bool IsValidExpressionData(MapBuilderData.TissueExpressionData expressionData)
		{
			return true;
		}

		/// <summary>
		/// Gets the link score from the given data.
		/// </summary>
		/// <returns>The link scores.</returns>
		/// <param name="corrVectors">Link data.</param>
		/// <param name="count">Total number of links.</param>
        protected override MapBuilderData.BaseLinkScore[] GetCorrelations(System.Collections.Generic.IEnumerable<MapBuilderData.BaseLinkScore> corrVectors, int count)
		{
			return corrVectors.ToArray();
		}

        public bool GeneMap { get; set; }

        public void FilterMap(MapLinkFilter.LinkType linkType)
        {
            var mapLinkFilter = new MapLinkFilter
            {
                LinkTypeFilter = linkType,
            };

            var map = TssRegulatoryMap.LoadMap(this.MapFileName, mapLinkFilter);

            var outputMapFile = this.MapFileName.Replace(".bed", "." + linkType.ToString() + ".bed");

            var tssSet = this.GeneMap ? 
                this.ExpressionFiles.First().Value.Genes :
                this.ExpressionFiles.First().Value.Transcripts;

            WriteMap(map, tssSet, this.HistoneName, outputMapFile);
        }

        public static void WriteMap(
            TssRegulatoryMap map, 
            Dictionary<string, Location> tssSet, 
            string histoneName, 
            string outputMapFile)
        {
            using (TextWriter tw = new StreamWriter(outputMapFile))
            {
                foreach (var link in map.Links)
                {
                    var tss = tssSet[link.TranscriptName];

                    tw.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", 
                        tss.Chromosome, 
                        tss.DirectionalStart, 
                        tss.DirectionalStart, 
                        link.TranscriptName, 
                        1.0, 
                        tss.Strand,
                        link.LocusName,
                        link.Correlation, 
                        link.ConfidenceScore, 
                        link.LinkLength, 
                        histoneName,
                        tss.AlternateName);
                }
            }
        }

		/// <summary>
		/// Formats the score into a bed file line.
		/// </summary>
		/// <returns>The score bed file line.</returns>
		/// <param name="scoreData">Score data.</param>
		/// <param name="location">Location.</param>
		/// <param name="distance">Distance.</param>
        protected override string FormatScoreTsvLine(MapBuilderData.BaseLinkScore scoreData, Genomics.Location location, int distance)
		{
			return string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}", 
				location.Chromosome, 
				location.DirectionalStart, 
				location.DirectionalStart, 
				scoreData.Tss, 
				1.0, 
				location.Strand, 
				scoreData.Locus, 
				0.0, 
				1.0, 
				distance, 
                "None",
				location.AlternateName);
		}

        public class Executor : IAnalysisExecutor<NullMapBuilder, NullMapBuilder.Executor.Arguments>
        {
            public enum Arguments
            {
                Tissue,
                TissueSources,
                Config,
                RnaSource,
                Stage,
                RemoveGeneBodyTranscripts,
                MapType,
                LocusFileName,
                MapFileName,
	 	OutDir,

                /// <summary>
                /// Flag indicating the map file name maps genes
                /// </summary>
                GeneMap,
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

                /// <summary>
                /// Filter a map
                /// </summary>
                FilterMap,
            }

            /// <summary>
            /// The options.
            /// </summary>
            /// <value>The options.</value>
            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { Arguments.Tissue, "Cell line whose Loci will be mapped" },
                        { Arguments.TissueSources, "Optional CSV of cell lines whose expression files should be used as TSS targets" },
                        { Arguments.Config, "XML configuration file name" },
                        { Arguments.RnaSource, "RNA-seq source whose TSSes will be used as putative Locus targets" },
                        { Arguments.Stage, "Stage in the mapping process: {Preprocess, BuildMap}" },
                        { Arguments.RemoveGeneBodyTranscripts, "Optional parameter to remove Loci within the specified distance of the transcript body (TSS to TTS)" },
                        { Arguments.MapType, "Parameter specifying type of connectivity {" + string.Join(", ", Enum.GetNames(typeof(MapLinkFilter.LinkType))) + "}" },
                        { Arguments.LocusFileName, "Optional manual override of Locus filename" },
                        { Arguments.MapFileName, "Optional manual override of map filename" },
                        { Arguments.OutDir, "Optional manual override of output directory" },
                        { Arguments.GeneMap, "Flag indicating the map file name maps genes" },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis
            /// </summary>
            /// <param name = "commandArgs"></param>
            override public void Execute(Args commandArgs)
            {
                Console.WriteLine(commandArgs.StringEnumArgs[Arguments.Tissue]);

                var nullMapBuilder = new NullMapBuilder(commandArgs.StringEnumArgs[Arguments.Config]);
                IMapBuilder builder = IUnknown.QueryInterface<IMapBuilder>(nullMapBuilder);

                if (commandArgs.StringEnumArgs.ContainsKey(Arguments.LocusFileName))
                {
                    builder.LocusFileName = commandArgs.StringEnumArgs[Arguments.LocusFileName];
                }

                if (commandArgs.StringEnumArgs.ContainsKey(Arguments.MapFileName))
                {
                    builder.MapFileName = commandArgs.StringEnumArgs[Arguments.MapFileName];
                    nullMapBuilder.GeneMap = commandArgs.Flags.Contains(Arguments.GeneMap.ToString());
                }

                builder.OutDir = commandArgs.StringEnumArgs.ContainsKey(Arguments.OutDir) ?
		    commandArgs.StringEnumArgs[Arguments.OutDir] :
                    builder.OutDir = ".";

                builder.Tissue = commandArgs.StringEnumArgs[Arguments.Tissue];
                builder.RnaSource = commandArgs.StringEnumArgs[Arguments.RnaSource];
                builder.LocusFilterRange = commandArgs.IntArgs.ContainsKey(Arguments.RemoveGeneBodyTranscripts.ToString()) ? 
                    commandArgs.IntArgs[Arguments.RemoveGeneBodyTranscripts.ToString()] :
                    0;

                if (commandArgs.StringEnumArgs.ContainsKey(Arguments.TissueSources))
                {
                    builder.TissueSources = commandArgs.StringEnumArgs[Arguments.TissueSources].Split('\t');
                }

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

                    case Stage.FilterMap:
                        nullMapBuilder.FilterMap((MapLinkFilter.LinkType)Enum.Parse(typeof(MapLinkFilter.LinkType), commandArgs.StringEnumArgs[Arguments.MapType]));
                        break;

                    default: 
                        throw new Exception("Invalid stage argument: " + commandArgs.StringEnumArgs[Arguments.Stage] + "\nValid values: Preprocess, BuildMap");
                }


            }

            /// <summary>
            /// Gets the description.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Map builder that maps all Loci to TSSes within the specified maximum range. Used to build a background model for all possible regulatory links";
                }
            }
        }

		private class NullMapBuilderProxy : IMapBuilder
		{
            private readonly NullMapBuilder instance;

			public NullMapBuilderProxy(NullMapBuilder data)
			{
				this.instance = data;
			}

            public string Tissue { get { return this.instance.Tissue; } set { this.instance.Tissue = value; } }
            public string[] TissueSources { get { return this.instance.TissueSources; } set { this.instance.TissueSources = value; } }
            public string LocusFileName { get { return this.instance.LocusFileName; } set { this.instance.LocusFileName = value; } }
            public string MapFileName { get { return this.instance.MapFileName; } set { this.instance.MapFileName = value; } }
            public string OutDir { get { return this.instance.OutDir; } set { this.instance.OutDir = value; } }
	    public string HistoneName { get { return this.instance.HistoneName; } set { this.instance.HistoneName = value; } }
	    public string RnaSource { get { return this.instance.RnaSource; } set { this.instance.RnaSource = value; } }
            public bool UseGenes { get { return this.instance.UseGenes; } set { this.instance.UseGenes = value; } }
            public int LocusFilterRange { get { return this.instance.LocusFilterRange; } set { this.instance.LocusFilterRange = value; } }
			public void PreprocessExpression() { this.instance.PreprocessExpression(); }
            public void BuildMap() { this.instance.BuildMap(); }
		}
    }
}

