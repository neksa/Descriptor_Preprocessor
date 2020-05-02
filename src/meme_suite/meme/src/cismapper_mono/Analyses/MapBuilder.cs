//--------------------------------------------------------------------------------
// File: MapBuilder.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using System.IO;
    using Genomics;
    using Shared;

	namespace MapBuilderData
	{
		/// <summary>
		/// Base link score.
		/// </summary>
		public class BaseLinkScore
		{
			public string Tss { get; set; }
			public string Locus { get; set; }
            		public double Score { get; set; }
		}

		/// <summary>
		/// Cell line expression data.
		/// </summary>
		public struct TissueExpressionData 
		{
			public string Tss { get; set; }
			public Genomics.Location Location { get; set; }
			public Dictionary<string, double> Expression { get; set; }
		}
	}

	/// <summary>
	/// Map building interface
	/// </summary>
	public interface IMapBuilder
	{
        /// <summary>
        /// Gets or sets the cell line.
        /// </summary>
        /// <value>The cell line.</value>
	string Tissue { get; set; }

        /// <summary>
        /// Gets or sets the cell line sources used for mapping expression data and TSS positions.
        /// </summary>
        /// <value>The cell line sources.</value>
        string[] TissueSources { get; set; }

        /// <summary>
        /// Gets or sets the name of the histone.
        /// </summary>
        /// <value>The name of the histone.</value>
	string HistoneName { get; set; }

        /// <summary>
        /// Gets or sets the name of the locus file.
        /// </summary>
        /// <value>The name of the locus file.</value>
        string LocusFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the map file.
        /// </summary>
        /// <value>The name of the map file.</value>
        string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the name of the output directory.
        /// </summary>
        /// <value>The name of the output directory.</value>
        string OutDir { get; set; }

        /// <summary>
        /// Gets or sets the rna source.
        /// </summary>
        /// <value>The rna source.</value>
	string RnaSource { get; set; }

        /// <summary>
        /// Gets or sets a value indicating the distance to remove loci withing a minimum distance of an annotated transcript
        /// </summary>
        /// <value><c>true</c> if filter proximal loci; otherwise, <c>false</c>.</value>
        int LocusFilterRange { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Analyses.IMapBuilder"/> use genes.
        /// </summary>
        /// <value><c>true</c> if use genes; otherwise, <c>false</c>.</value>
        bool UseGenes { get; set; }

        /// <summary>
        /// Preprocesses the expression.
        /// </summary>
	void PreprocessExpression();

        /// <summary>
        /// Builds the map.
        /// </summary>
        void BuildMap();

	}

    /// <summary>
    /// Core methods of a map builder.
    /// </summary>
	abstract public class MapBuilder<TLinkData, TLinkScore> : BaseData where TLinkData : MapBuilderData.BaseLinkScore where TLinkScore : MapBuilderData.BaseLinkScore
    {
	// Private Members ===============================================================================================//

        /// <summary>
        /// The cell line expression files.
        /// </summary>
        private Dictionary<string, string> expressionFilenames;

        /// <summary>
        /// The histone filenames.
        /// </summary>
	private Dictionary<string, string> histoneFilenames;

        /// <summary>
        /// The cell line sources.
        /// </summary>
        private string[] tissueSources;

        /// <summary>
        /// The minimum expression test.
        /// </summary>
	private Func<double, bool> minExpressionTest;

        /// <summary>
        /// The annotation set.
        /// </summary>
	private IAnnotation annotationSet;

        /// <summary>
        /// The transcript set.
        /// </summary>
	private HashSet<string> transcriptSet;

        /// <summary>
        /// The bad locus set.
        /// </summary>
        private HashSet<string> badLocusSet;

        /// <summary>
        /// The expression files.
        /// </summary>
	private Dictionary<string, IExpressionData> expressionFiles;

        /// <summary>
        /// The tss locations.
        /// </summary>
	private Dictionary<string, Genomics.Location> transcriptLocations;

        /// <summary>
        /// The tss expression.
        /// </summary>
	private Dictionary<string, Dictionary<string, double>> transcriptTissueExpression;

        /// <summary>
        /// The transcript chromosomes.
        /// </summary>
        private Dictionary<string, string> transcriptChromosomes;

        /// <summary>
        /// The chromosome transcripts.
        /// </summary>
        private ILookup<string, string> chromosomeTranscripts;

        /// <summary>
        /// The transcript locus map.
        /// </summary>
        private Dictionary<string, List<Genomics.Location>> transcriptLocusMap;

        /// <summary>
        /// The histone locus data.
        /// </summary>
        private Dictionary<string, Dictionary<string, double>> locusData;

	/// <summary>
	/// The processed locus filename.
	/// </summary>
        private string processedLociFilename;

	/// <summary>
	/// The processed transcript location filename.
	/// </summary>
	private string processedTranscriptLocationFilename;

	/// <summary>
	/// The processed expression filename.
	/// </summary>
	private string processedExpressionFilename;

	/// <summary>
	/// The processed histone filename.
	/// </summary>
	private string processedHistoneFilename;

        /// <summary>
        /// The transcript types.
        /// </summary>
        private string[] transcriptTypes;

	/// <summary>
	/// The omitted cell lines.
	/// </summary>
	private string[] omittedTissues;

	/// <summary>
	/// The current stage.
	/// </summary>
	private Stage currentStage;

	// Enumerations ===============================================================================================//

	private enum Stage
	{
		Preprocess,
		BuildMap
	}

	// Constructors ===============================================================================================//

        const string OutputRoot = "../temp/results/MapGeneration/";

	/// <summary>
	/// Initializes a new instance of the <see cref="Genomics.MapBuilder`2"/> class.
	/// </summary>
	/// <param name="xmlFile">Xml file.</param>
	protected MapBuilder(string xmlFile, string[] omittedTissues)
            : base(OutputRoot, OutputRoot, OutputRoot, xmlFile)
        {
			this.omittedTissues = omittedTissues;
	}

	/// <summary>
	/// Initializes a new instance of the <see cref="Genomics.MapBuilder`2"/> class.
	/// </summary>
	/// <param name="xmlFile">Xml file.</param>
	protected MapBuilder(string xmlFile)
    : base(OutputRoot, OutputRoot, OutputRoot, xmlFile)
	{
		this.omittedTissues = new string[0];
	}

	// Abstract Members ===============================================================================================//

	/// <summary>
	/// Gets the name of the mapping method name
	/// </summary>
	/// <value>The name of the method.</value>
	abstract public string MethodName { get; }

	/// <summary>
	/// Construct the data used to establish a confidence score for the link
	/// </summary>
	/// <returns>The link data.</returns>
	/// <param name="tss">Tss.</param>
	/// <param name="locus">Locus name.</param>
	/// <param name="expressionByTissue">Expression by cell line.</param>
	/// <param name="histoneByTissue">Histone by cell line.</param>
	abstract protected TLinkData BuildLinkData(
		string tss, 
		string locus, 
		Dictionary<string, double> expressionByTissue,
		Dictionary<string, double> histoneByTissue);

	/// <summary>
	/// Determines whether specified link data is valid
	/// </summary>
	/// <returns><c>true</c> if specified link data is valid; otherwise, <c>false</c>.</returns>
	/// <param name="data">Data.</param>
	abstract protected bool IsValidLinkData(TLinkData data);

	/// <summary>
	/// Determines whether the specified expressionData is valid
	/// </summary>
	/// <returns><c>true</c> if the specified expressionData is valid; otherwise, <c>false</c>.</returns>
	/// <param name="expressionData">Expression data.</param>
	abstract protected bool IsValidExpressionData(MapBuilderData.TissueExpressionData expressionData);

	/// <summary>
	/// Gets the link score from the given data.
	/// </summary>
	/// <returns>The link scores.</returns>
	/// <param name="corrVectors">Link data.</param>
	/// <param name="count">Total number of links.</param>
	abstract protected TLinkScore[] GetCorrelations(IEnumerable<TLinkData> corrVectors, int count);

	/// <summary>
	/// Formats the score into a bed file line.
	/// </summary>
	/// <returns>The score bed file line.</returns>
	/// <param name="scoreData">Score data.</param>
	/// <param name="location">Location.</param>
	/// <param name="distance">Distance.</param>
	abstract protected string FormatScoreTsvLine(TLinkScore scoreData, Genomics.Location location, int distance);

	// Public Methods ===============================================================================================//

        /// <summary>
        /// Preprocesses the expression.
        /// </summary>
	public void PreprocessExpression()
	{
		this.currentStage = Stage.Preprocess;

		BedFile.ToFileBed6(this.LocusFile.Locations.Values.Where(x => !this.BadLocusSet.Contains(x.Name)).ToList(), this.ProcessedLociFilename);

		BedFile.ToFileBed6(this.TranscriptLocations.Values.ToList(), this.ProcessedTranscriptLocationFilename);

		using (TextWriter tw = Helpers.CreateStreamWriter(this.ProcessedExpressionFilename))
		{
			tw.WriteLine("TSS_ID\t" + string.Join("\t", this.TissueSources));
			tw.WriteLine(string.Join(
				"\n", 
				this.TranscriptExpression
					.Select(x => x.Key + "\t" + string.Join(
						"\t", 
					this.TissueSources.Select(c => x.Value.ContainsKey(c) ? x.Value[c].ToString() : "NaN")))));
		}
		using (TextWriter tw2 = Helpers.CreateStreamWriter(this.ProcessedHistoneFilename))
		{
			tw2.WriteLine("RE_Locus\t" + string.Join("\t", this.TissueSources));
			tw2.WriteLine(string.Join(
				"\n", 
				this.LocusData
					.Select(x => x.Key + "\t" + string.Join(
						"\t", 
					this.TissueSources.Select(c => x.Value.ContainsKey(c) ? x.Value[c].ToString() : "NaN")))));
		
		}
	}
            
        private string mapFileName;
        public string MapFileName
        {
            get 
            {
                return Helpers.CheckInit(
                    ref this.mapFileName,
                    () => string.Format(
                        "../temp/results/MapGeneration/{0}{1}{2}.{3}.{4}.{5}.{6}{7}.bed", 
                        this.OmittedTissues.Length > 0 ? "Sans" + string.Join("_", this.OmittedTissues) + "." : "Full.",
                        string.IsNullOrEmpty(this.Tissue) ? string.Empty : string.Format("{0}.", this.Tissue), 
                        this.MethodName,
                        this.RnaSource, 
                        this.HistoneName, 
                        this.TranscriptTypesString,
                        this.ConfigData.StringValues["MaximumLinkDistance"],
                        this.ConfigData.StringValues.ContainsKey("Batch") ? ".Batch_" + this.ConfigData.StringValues["Batch"] : string.Empty));
            }

            set
            {
                this.mapFileName = value;
            }
        }

	private string outDir;
        public string OutDir 
        {
            get 
	    { 
              return this.outDir; 
            }
	    set
	    {
		this.outDir = value;
	    } 
        }
	

	/// <summary>
	/// Builds a map for a given cell line's loci. 
	/// Optionally omits source line from data.
	/// </summary>
        public void BuildMap()
	{
	    this.currentStage = Stage.BuildMap;
	    this.LocusFileName = this.ProcessedLociFilename;

	    Console.WriteLine("Building drm-tss pairs using histone " + this.HistoneName + ", expression type " + this.RnaSource + ", and Tissues:");
            Console.WriteLine("\t" + string.Join("\n\t", this.TissueSources));

            using (TextWriter tw = Helpers.CreateStreamWriter(this.MapFileName, true)) // append to map file
            {

                foreach (var chr in this.ChromosomeTranscripts)
                {
                    //GC.Collect();
		    Console.WriteLine("Chromosome {0}\n", chr.Key);

                    var corrVectors = BuildLinkData(chr.Key);

                    int count = corrVectors.Count();

                    if (count > 0)
                    {
                        var cordata = GetCorrelations(corrVectors, count);

                        Console.WriteLine("########");
                        Console.WriteLine(
                            "\t### {0} tss-drm pairs with valid data; {1} tsses; {2} drms", 
                            cordata.Length, 
                            cordata.ToLookup(x => x.Tss).Count, cordata.ToLookup(x => x.Locus).Count);

                        foreach (var entry in cordata)
                        {
                            var tss = this.TranscriptLocations[entry.Tss];
                            var locus = this.LocusFile.Locations[entry.Locus];
                            int distToLocusStart = tss.DirectionalStart - this.LocusFile.Locations[entry.Locus].Start;
                            int distToLocusEnd = tss.DirectionalStart - this.LocusFile.Locations[entry.Locus].End;
                            int distance = Math.Sign(distToLocusStart) * Math.Min(Math.Abs(distToLocusStart), Math.Abs(distToLocusEnd));

                            // If the sign of the distance to each end of the locus is different,
                            // then the TSS is in the MIDDLE of the locus, so the distance is zero.
                            if (Math.Sign(distToLocusStart) != Math.Sign(distToLocusEnd))
                            {
                                distance = 0;
                            }

                            tw.WriteLine(FormatScoreTsvLine(entry, tss, distance));
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Gets or sets a value indicating the distance to remove loci withing a minimum distance of an annotated transcript
        /// </summary>
        /// <value><c>true</c> if filter proximal loci; otherwise, <c>false</c>.</value>
        public int LocusFilterRange { get; set; }

	// Protected Properties ===============================================================================================//

        /// <summary>
        /// Gets a value indicating whether to omit expression from the mapped cell line in mapping
        /// </summary>
        /// <value><c>true</c> if omit loci tissue source; otherwise, <c>false</c>.</value>
		protected string[] OmittedTissues
        {
            get
            {
				return Helpers.CheckInit(
					ref this.omittedTissues,
					() => 
					{
						if (this.ConfigData.CsvValues.ContainsKey("OmitLocusSource"))
						{
							return this.ConfigData.CsvValues["OmitLocusSource"].ToArray();
						}

						return  new string[0];
					});
            }
        }

        /// <summary>
        /// Gets the expression filenames.
        /// </summary>
        /// <value>The expression filenames.</value>
        protected Dictionary<string, string> ExpressionFilenames
        {
            get
            {
				if (this.expressionFilenames == null)
				{
					this.expressionFilenames = this.TissueSources
                        .Select(cl => new
                        {
                            Tissue = cl,
							Filename = string.Format("{0}/{1}/{2}.{3}", ExpressionRoot, cl, this.RnaSource, this.ConfigData.StringValues["ExpressionFileType"])
                        })
                        .Where(x => File.Exists(x.Filename))
						.ToDictionary(x => x.Tissue, x => x.Filename);
				}

				return this.expressionFilenames;
            }
        }

        /// <summary>
        /// Gets the minimum feature count.
        /// </summary>
        /// <value>The minimum feature count.</value>
        protected int MinFeatureCount
        {
            get
            {
                return this.ConfigData.IntValues ["MinFeatureCount"];
            }
        }

        /// <summary>
        /// Gets the cell line sources.
        /// </summary>
        /// <value>The cell line sources.</value>
        public string[] TissueSources
        {
            get
            {
				return this.tissueSources ?? (this.tissueSources = this.ConfigData.CsvValues ["Tissues"]
					.Where(cl => !string.IsNullOrEmpty(this.Tissue) && !this.OmittedTissues.Contains(cl))
					.ToArray());
            }

            set
            {
                this.tissueSources = value;
            }
        }

        /// <summary>
        /// Gets the minimax expression value.
        /// </summary>
        /// <value>The minimax expression value.</value>
        protected double MinimaxExpressionValue
        {
            get
			{
				return double.Parse(this.ConfigData.StringValues ["MinimaxExpression"]);
            }
        }

        /// <summary>
        /// Gets the histone root.
        /// </summary>
        /// <value>The histone root.</value>
        protected string HistoneRoot 
        {
            get
            {
                return this.ConfigData.StringValues ["HistoneRoot"];
            }
        }

        /// <summary>
        /// Gets the expression root.
        /// </summary>
        /// <value>The expression root.</value>
        protected string ExpressionRoot
        {
            get
            {
                return this.ConfigData.StringValues ["ExpressionRoot"];
            }
        }

        /// <summary>
        /// Gets the minimum expression test.
        /// </summary>
        /// <value>The minimum expression test.</value>
	protected Func<double, bool> MinExpressionTest
	{
	    get
		{
		    if (this.minExpressionTest == null)
		    {
		       if (this.ConfigData.StringValues.ContainsKey("MinExpression"))
	    	       {
			   double minExpression = double.Parse(this.ConfigData.StringValues ["MinExpression"]);
			   this.minExpressionTest = x => x > minExpression;
	    	       } else
	               {
			   this.minExpressionTest = x => true;
			}
		    }

		    return this.minExpressionTest;
	        }
	}

        /// <summary>
        /// Gets the annotation set.
        /// </summary>
        /// <value>The annotation set.</value>
	protected IAnnotation AnnotationSet
        {
	    get
	    {
		return Helpers.CheckInit(
		    ref this.annotationSet, () =>
		        {   
			    Console.WriteLine("Getting annotation set");
			    IAnnotation annotationSet = null;
			    switch (ConfigData.StringValues ["Annotation/Type"])
			    {
				case "Gencode":
				    annotationSet = this.annotationSet = IUnknown.QueryInterface<IAnnotation>(this.UseGenes ?
					GtfGencodeFile.LoadGenes(ConfigData.StringValues ["Annotation/FileName"]) :
					GtfGencodeFile.LoadTranscripts(ConfigData.StringValues ["Annotation/FileName"])
				    );
				    break;
				case "RefSeq":
				    annotationSet = this.annotationSet = IUnknown.QueryInterface<IAnnotation>(this.UseGenes ?
				    GtfGencodeFile.LoadGenes(ConfigData.StringValues ["Annotation/FileName"]) : 
				    GtfGencodeFile.LoadTranscripts(ConfigData.StringValues ["Annotation/FileName"]));
				    break;
				default:
				    throw new Exception("Invalid annotation type: " + ConfigData.StringValues ["Annotation/Type"]);
			     }
			     annotationSet.ValidTranscriptTypes = this.TranscriptTypes;
			     return annotationSet;
			 }
                );
	    }
        }
        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Analyses.MapBuilder`2"/> uses genes.
        /// </summary>
        /// <value><c>true</c> if use genes; otherwise, <c>false</c>.</value>
        public bool UseGenes { get; set; }

        /// <summary>
        /// Gets the transcript set.
        /// </summary>
        /// <value>The transcript set.</value>
        protected HashSet<string> TranscriptSet
	{
		get
		{
			return this.transcriptSet ?? (this.transcriptSet = this.UseGenes ? 
			this.AnnotationSet.ValidGenes : 
			this.AnnotationSet.ValidTranscripts);
		}
	}

        /// <summary>
        /// Gets the transcript types.
        /// </summary>
        /// <value>The transcript types.</value>
        protected string[] TranscriptTypes
        {
            get
            {
                return this.transcriptTypes ?? (this.transcriptTypes = this.ConfigData.CsvValues["TranscriptTypes"].ToArray());
            }
        }

        /// <summary>
        /// String giving the transcript type set used for map generation delimited by '-'
        /// </summary>
        protected string TranscriptTypesString
        {
            get
            {
                return string.Join("-", this.TranscriptTypes.OrderBy(x => x));
            }
        }

        /// <summary>
        /// Gets the bad locus set.
        /// </summary>
        /// <value>The bad locus set.</value>
		protected HashSet<string> BadLocusSet
		{
			get
			{
                return Helpers.CheckInit(
                    ref this.badLocusSet,
                    () =>
                    {
                        if (this.LocusFilterRange == 0)
                        {
                            Console.WriteLine("Omitting transcript filtering");
                            return new HashSet<string>();
                        }

                        Console.WriteLine("Getting loci within " + this.LocusFilterRange + " of any annotated gene bodies");
                        return new HashSet<string>(this.LocusFile.ChromosomeOrderedLocations
                            .Select(x => TRFScorer.IntersectLeft(this.LocusFile, this.AnnotationSet, x.Key, this.LocusFilterRange))
                            .SelectMany(x => x.Select(y => y.Name)));
                    });
			}
		}

		/// <summary>
		/// Gets the processed locus filename.
		/// </summary>
		/// <value>The processed locus filename.</value>
		protected string ProcessedLociFilename
		{
			get
			{
				return Helpers.CheckInit(
					ref this.processedLociFilename,
                    () => string.Format(this.OutDir + "/" + ConfigData.StringValues["ProcessedFileFormat/Loci"], this.Tissue, this.RnaSource, this.HistoneName, this.TranscriptTypesString));
			}
		}

        /// <summary>
        /// Gets the processed transcript location filename.
        /// </summary>
        /// <value>The processed transcript location filename.</value>
		protected string ProcessedTranscriptLocationFilename
		{
			get
			{
				return Helpers.CheckInit(
					ref this.processedTranscriptLocationFilename,
                    () => string.Format(this.OutDir + "/" + ConfigData.StringValues["ProcessedFileFormat/TranscriptLocation"], this.Tissue, this.RnaSource, this.HistoneName, this.TranscriptTypesString));
			}
		}

        /// <summary>
        /// Gets the processed expression filename.
        /// </summary>
        /// <value>The processed expression filename.</value>
	protected string ProcessedExpressionFilename
	{
	  get
	  {
	    return Helpers.CheckInit(
	      ref this.processedExpressionFilename,
              () => string.Format(this.OutDir + "/" + ConfigData.StringValues["ProcessedFileFormat/TranscriptExpression"], this.Tissue, this.RnaSource, this.HistoneName, this.TranscriptTypesString));
 	  }
	}

        /// <summary>
        /// Gets the processed histone filename.
        /// </summary>
        /// <value>The processed histone filename.</value>
	protected string ProcessedHistoneFilename
	{
	  get
	  {
	    return Helpers.CheckInit(
	      ref this.processedHistoneFilename,
              () => string.Format(this.OutDir + "/" + ConfigData.StringValues["ProcessedFileFormat/Histone"], this.Tissue, this.RnaSource, this.HistoneName, this.TranscriptTypesString));
	  }
        }

        /// <summary>

        /// <summary>
        /// Gets the maximum link distance.
        /// </summary>
        /// <value>The maximum link distance.</value>
        protected int MaximumLinkDistance
        {
            get
            {
                return this.ConfigData.IntValues["MaximumLinkDistance"];
            }
        }

        /// <summary>
        /// Gets the expression files.
        /// </summary>
        /// <value>The expression files.</value>
	protected Dictionary<string, IExpressionData> ExpressionFiles
	{
		get
		{
			return Helpers.CheckInit(
				ref this.expressionFiles,
				() => this.ExpressionFilenames.ToDictionary(x => x.Key, x =>
				  {
				    IExpressionData data = IExpressionData.LoadExpressionData(
				      Path.Combine(this.ConfigData.StringValues["ExpressionRoot"], x.Key, this.RnaSource + "." + this.ConfigData.StringValues ["ExpressionFileType"]), 
				      this.RnaSource, 
				      this.ConfigData.StringValues ["ExpressionFileType"], 
				      this.AnnotationSet
				    );
				    
				    Console.WriteLine("\ttranscript count: " + data.Transcripts.Count);

				    return data;
				  }
				 ));
		}
	}
 
        /// <summary>
        /// Gets the histone filenames.
        /// </summary>
        /// <value>The histone filenames.</value>
		protected Dictionary<string, string> HistoneFilenames
		{
			get
			{
				return Helpers.CheckInit(
					ref this.histoneFilenames,
					() =>
					{
						var hf = this.TissueSources
							.Select(cl => new 
								{
									Tissue = cl,
									HistoneFilename = GetHistoneModFilesForTissue(cl, this.HistoneRoot)
								})
							.Where(x => x.HistoneFilename != null)
							.ToDictionary(x => x.Tissue, x => x.HistoneFilename);

						if (hf.Count < this.MinFeatureCount)
						{
							string message = "Error: insufficient number of cell lines for histone modification " + this.HistoneName;
							Console.Error.WriteLine(message);
							throw new Exception(message);
						}

						return hf;
					});
			}
		}

		/// <summary>
		/// Gets the transcript locus map.
		/// </summary>
		/// <value>The transcript locus map.</value>
        protected Dictionary<string, List<Genomics.Location>> TranscriptLocusMap
        {
            get
            {
				return Helpers.CheckInit(
					ref this.transcriptLocusMap,
					() => this.TranscriptLocations
                        .Where(tss => this.LocusFile.ChromosomeOrderedLocations.ContainsKey(tss.Value.Chromosome))
                            .Select(tss => new 
                            {
                                Tss = tss.Key,
                                Loci= TRFScorer.GetOverlaps(
                                    new Location
                                    {
                                        Chromosome = tss.Value.Chromosome,
                                        Start = tss.Value.DirectionalStart - this.MaximumLinkDistance,
                                        End = tss.Value.DirectionalStart + this.MaximumLinkDistance, 
                                    },
                                    this.LocusFile.ChromosomeIndexedLocations,
                                    BedFile.IndexSize,
                                    this.LocusFile.MaxLocationSize[tss.Value.Chromosome])
                        /*Loci = TRFScorer.EvaluateRange(
								this.LocusFile.ChromosomeOrderedLocations [tss.Value.Chromosome].Where(locus => this.currentStage == Stage.Preprocess ? !this.BadLocusSet.Contains(locus.Name) : true).ToList(),
                                tss.Value.DirectionalStart - this.MaximumLinkDistance,
                                tss.Value.DirectionalStart + this.MaximumLinkDistance, 
                                    new List<Genomics.Location>(),
                                    (item, start) => item.End < start,
                                    (item, end) => item.Start <= end,
                                    (items, length, location) => 
                                    {
                                    items.Add(location);
                                    return items; 
                                })*/
							}).ToDictionary(x => x.Tss, x => x.Loci));
            }
        }

		/// <summary>
		/// Gets the histone locus data.
		/// </summary>
		/// <value>The histone locus data.</value>
        protected Dictionary<string, Dictionary<string, double>> LocusData
        {
            get
            {
                if (this.locusData != null)
                {
                    return this.locusData;
                }
                else
                {
                    if (this.UseHistoneData)
                    {
                        this.locusData = this.HistoneFilenames
                            .Select(cl =>
                        {
                            Console.WriteLine("Loading histone data from " + cl.Value);
                            BedFile peaks = new BedFile(
                                                    cl.Value, 
                                                    cl.Value.Contains("broad") ? BedFile.Bed6Plus3Layout : BedFile.Bed6Plus4Layout);

                            return this.LocusFile.Locations
    										.Where(locus => this.currentStage != Stage.Preprocess || !this.BadLocusSet.Contains(locus.Key))
                                            .Where(x => peaks.ChromosomeOrderedLocations.ContainsKey(x.Value.Chromosome))
                                            .Select(x => new
                                                {
                                                    Tissue = cl.Key,
                                                    Locus = x.Key,
                                                    Value = TRFScorer.EvaluateRange<double>(
                                peaks.ChromosomeOrderedLocations[x.Value.Chromosome],
                                x.Value.Start,
                                x.Value.End,
                                0,
                                (item, start) => item.End < start,
                                (item, end) => item.Start <= end,
                                (aggregateScore, length, peak) =>
                                {
                                    double d = 0.0;
                                    if (double.TryParse(peak.Data[peaks.FileLayout.SignalValue], out d))
                                    {
                                        // Get max peak value in range
                                        return d > aggregateScore ? d : aggregateScore;
                                    }
                                    return aggregateScore;
                                })
                                                });
                        })
                            .SelectMany(x => x)
                            .ToLookup(x => x.Locus, x => x)
                            .ToDictionary(x => x.Key, x => x.ToDictionary(y => y.Tissue, y => y.Value));
                    }
                    else
                    {
                        this.locusData = this.LocusFile.Locations
                            .Where(locus => this.currentStage != Stage.Preprocess || !this.BadLocusSet.Contains(locus.Key))
                            .ToDictionary(x => x.Key, x => this.TissueSources.ToDictionary(y => y, y =>
                            {
                                double value;
                                if (x.Value.Data.Length >= 7 && double.TryParse(x.Value.Data[6], out value))
                                {
                                    return value;
                                }
                                else
                                {
                                    return 0.0;
                                }
                            }));
                    }

                    return this.locusData;
                }
            }
        }

	/// <summary>
        /// Gets the tss locations.
        /// </summary>
        /// <value>The tss locations.</value>
	protected Dictionary<string, Genomics.Location> TranscriptLocations
	{
		get
		{
			if (this.transcriptLocations == null)
			{
				this.InitExpression();
			}

			return this.transcriptLocations;
		}
	}

        /// <summary>
        /// Gets the tss expression.
        /// </summary>
        /// <value>The tss expression.</value>
	protected Dictionary<string, Dictionary<string, double>> TranscriptExpression
	{
		get
		{
			if (this.transcriptTissueExpression == null)
			{
				this.InitExpression();
			}

			return this.transcriptTissueExpression;
		}
	}

        /// <summary>
        /// Gets the transcript chromosomes.
        /// </summary>
        /// <value>The transcript chromosomes.</value>
        protected Dictionary<string, string> TranscriptChromosomes
        {
            get
            {
                if (this.transcriptChromosomes == null)
                {
                    this.InitExpression();
                }

                return this.transcriptChromosomes;
            }
        }

        /// <summary>
        /// Gets the chromosome transcripts.
        /// </summary>
        /// <value>The chromosome transcripts.</value>
        protected ILookup<string, string> ChromosomeTranscripts
        {
            get
            {
                if (this.chromosomeTranscripts == null)
                {
                    this.InitExpression();
                }

                return this.chromosomeTranscripts;
            }
        }

	// Private Methods ===============================================================================================//

	/// <summary>
        /// Initializes the expression data backing properties TssLocations and TssExpression
        /// </summary>
	private void InitExpression()
	{

	    if (this.currentStage == Stage.Preprocess)
	    {
	        var data = this.ExpressionFiles
	            .Select(x => x.Value.Transcripts
		    .Where(loc => this.TranscriptSet.Contains(loc.Value.Name) && this.MinExpressionTest(loc.Value.Score))
		    .Select(loc => new 
		    {
			Tissue = x.Key,
			Location = loc.Value
		    }))
		    .SelectMany(x => x)
		    .ToLookup(x => x.Location.Name, x => x).Where(x => x.Any())
		    .Select(x =>
		    {
			var chromosomes = x.Select(y => y.Location.Chromosome).ToArray();
			var starts = x.Select(y => y.Location.DirectionalStart).ToArray();

			if ((chromosomes.Count(y => y == chromosomes[0]) != chromosomes.Length) || (starts.Count(y => y == starts[0]) != starts.Length))
			{
			    Console.Error.WriteLine("Error: inconsistent transcript data:\t" + x.Key + "\t" + chromosomes[0] + "\t" + string.Join(",", chromosomes) + "\t" + starts.First() + "\t" + string.Join(",", starts));
			}
			return new MapBuilderData.TissueExpressionData
			{
			    Tss = x.Key,
			    Location = x.First().Location,
			    Expression = x.ToLookup(y => y.Tissue, y => y.Location)
			    .ToDictionary(y => y.Key, y => y.Max(z => z.Score))
			};
		    })
		    .Where(IsValidExpressionData)
		    .ToArray();

		    this.transcriptLocations = data.ToDictionary(x => x.Tss, x => x.Location);
		    this.transcriptChromosomes = data.ToDictionary(x => x.Tss, x => x.Location.Chromosome);
		    this.chromosomeTranscripts = this.transcriptChromosomes.ToLookup(x => x.Value, x => x.Key);
		    this.transcriptTissueExpression = data.ToDictionary(x => x.Tss, x => x.Expression);
	    }
	    else
	    {
	        this.transcriptLocations = new BedFile(this.ProcessedTranscriptLocationFilename, BedFile.Bed3Layout).Locations;
                this.transcriptChromosomes = this.transcriptLocations.ToDictionary(x => x.Key, x => x.Value.Chromosome);
                this.chromosomeTranscripts = this.transcriptChromosomes.ToLookup(x => x.Value, x => x.Key);

		using (TextReader tr = new StreamReader(this.ProcessedExpressionFilename))
		{
		   var Tissues = tr.ReadLine().Split('\t').ToArray();
		   this.transcriptTissueExpression = tr.ReadToEnd().Split('\n').Where(line => !string.IsNullOrEmpty(line))
		       .Select(line =>
		       {
		           var fields = line.Split('\t');
			   return new {
			      Transcript = fields[0],
			      Expression = Tissues.Select((x, i) => new 
			      {
				      Tissue = x,
				      ExpressionString = fields[i + 1]
			      })
			      .Where(x => x.ExpressionString != "NaN")
			      .ToDictionary(x => x.Tissue, x => double.Parse(x.ExpressionString))
		            };
			})
			.ToDictionary(x => x.Transcript, x => x.Expression);
		  }
	    }
	}

        /// <summary>
        /// Builds the correlation vectors.
        /// </summary>
        /// <returns>The correlation vectors.</returns>
        private List<TLinkData> BuildLinkData(string chromosome)
        {
            Console.Write("Building correlation vectors ");

            List<TLinkData> data = new List<TLinkData>();

            int i = 0;
            Console.WriteLine(chromosome);
            foreach (var tss in this.TranscriptLocusMap.Where(x => this.TranscriptChromosomes[x.Key] == chromosome)) //.Where(tss => this.TranscriptExpression.ContainsKey(tss.Key)))
            {
                if (this.TranscriptExpression.ContainsKey(tss.Key))
                {
                    Dictionary<string, double> tssExpression = this.TranscriptExpression[tss.Key];
                    foreach (var locus in tss.Value.Where(locus => this.LocusData.ContainsKey(locus.Name)))
                    {
                        var linkData = BuildLinkData(
                            tss.Key, 
                            locus.Name, 
                            tssExpression, 
                            this.LocusData[locus.Name]
                                .Where(tissue => tssExpression.ContainsKey(tissue.Key))
                                .ToDictionary(x => x.Key, x => x.Value));

                        if (IsValidLinkData(linkData))
                        {
                            data.Add(linkData);
                            i++;
                            if (i % 100000 == 0)
                            {
                                Console.Write(".");
                            }
                        }
                    }
                }
            }
            Console.WriteLine();

            return data;

            /*return this.TranscriptLocusMap.Where(tss => this.TranscriptExpression.ContainsKey(tss.Key)).Select(tss => 
            {
                Dictionary<string, double> tssExpression = this.TranscriptExpression[tss.Key];
                return tss.Value
					.Where(drm => this.HistoneLocusData.ContainsKey(drm.Name))
					.Select(drm => BuildLinkData(
						tss.Key, 
						drm.Name, 
						tssExpression, 
						this.HistoneLocusData[drm.Name].Where(data => tssExpression.ContainsKey(data.Key)).ToDictionary(x => x.Key, x => x.Value)))
					.Where(IsValidLinkData);
			}).SelectMany(x => x)
                    .ToList();*/
        }

        /// <summary>
        /// Writes the map.
        /// </summary>
        /// <param name="cordata">Cordata.</param>
        /// <param name = "histoneNameOverride"></param>
        /*private void WriteMap(IEnumerable<IEnumerable<TLinkScore>> cordata, string histoneNameOverride)
        {

                int i = 0;
                foreach (var chrcordata in cordata)
                {
                    foreach (var entry in chrcordata)//for (int i = 0; i < cordata.Count; i++)
                    {
                        var location = this.TranscriptLocations[entry.Tss];
                        int drmStartDist = location.DirectionalStart - this.LocusFile.Locations[entry.Locus].Start;
                        int drmEndDist = location.DirectionalStart - this.LocusFile.Locations[entry.Locus].End;
                        int distance = drmStartDist;
                        if (Math.Abs(drmEndDist) < Math.Abs(drmStartDist))
                        {
                            distance = drmEndDist;
                        }
                        //allTssDrmPairs.Add(string.Format("{0}#{1}", entry.Tss, entry.Drm));
                        tw.WriteLine(FormatScoreTsvLine(entry, location, distance));

                        i++;
                    }
                }
            }
        }*/

		/// <summary>
		/// Gets the histone modification files for a cell line.
		/// </summary>
		/// <returns>The histone mod files for cell line.</returns>
		/// <param name="tissue">Cell line.</param>
		/// <param name="histoneRoot">Histone root.</param>
		private string GetHistoneModFilesForTissue(string tissue, string histoneRoot)
		{
			string histoneDirectory = string.Format("{0}/{1}", histoneRoot, tissue);
			string[] histoneFiles = Directory.GetFiles(histoneDirectory, "*" + this.HistoneName + "*Peak");

			return histoneFiles.Length > 0 ? histoneFiles[0] : null;

		}
	}
}

