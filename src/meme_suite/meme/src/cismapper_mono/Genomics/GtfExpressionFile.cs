//--------------------------------------------------------------------------------
// File: GtfExpressionFile.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Text.RegularExpressions;
    using System.Linq;
    using Shared;

    /// <summary>
    /// Wrapper for a sequence collection such as a genome
    /// </summary>
	public class GtfExpressionFile : BedFile//, IExpressionData
    {
        public enum ExpressionType
        {
            LongPap,
            LongPapMouse,
            LongPam,
            Short,
            Cage
        };

        private ExpressionType expressionType;

        /// <summary>
        /// Expression level of the last parsed line
        /// </summary>
        private double expressionLevel;

        /// <summary>
        /// Expression field names
        /// </summary>
        private static Dictionary<ExpressionType, Regex> ExpressionFieldFormat = new Dictionary<ExpressionType, Regex>
        {
            { ExpressionType.LongPap, new Regex("[RF]PKM[12]?") },
            { ExpressionType.LongPapMouse, new Regex("[RF]PKM[12]?") },
            { ExpressionType.LongPam, new Regex("[RF]PKM[12]?") },
            { ExpressionType.Short, new Regex("rpm[12]") },
            { ExpressionType.Cage, new Regex("rpm[12]") },
        };

        /// <summary>
        /// Transcript field names
        /// </summary>
        private static Dictionary<ExpressionType, string> TranscriptFieldFormat = new Dictionary<ExpressionType, string>
        {
            { ExpressionType.LongPap, "transcript_id" },
            { ExpressionType.LongPapMouse, "transcript_id" },
            { ExpressionType.LongPam, "transcript_id" },
            { ExpressionType.Short, "transcript_id" },
            { ExpressionType.Cage, "trlist" },
        };

        /// <summary>
        /// Entry types for expression
        /// </summary>
        private static Dictionary<ExpressionType, string> EntryTypeFieldFormat = new Dictionary<ExpressionType, string>
        {
            { ExpressionType.LongPap, "transcript" },
            { ExpressionType.LongPapMouse, "start_codon" },
            { ExpressionType.LongPam, "transcript" },
            { ExpressionType.Short, "exon" },
            { ExpressionType.Cage, "TSS" },
        };

		/// <summary>
		/// Tests if a transcript name should be included in the expression set
		/// </summary>
		private Func<string, bool> isValidTranscriptName;

		/// <summary>
		/// Initializes a new instance of the <see cref="Genomics.GtfExpressionFile"/> class.
		/// </summary>
		/// <param name="type">Type.</param>
		/// <param name="filename">Filename.</param>
		public GtfExpressionFile(ExpressionType type, string filename, IAnnotation annotation)
            : base(filename, new Layout 
                    {
                        Chromosome = 0,
                        Start = 3,
                        End = 4,
                        Strand = 6,
                        Name = 8
                    },
                false)
        {
            this.expressionType = type;
	    this.Annotation = annotation;

            ParseFile();
        }

	/// <summary>
	/// Initializes a new instance of the <see cref="Genomics.GtfExpressionFile"/> class.
	/// </summary>
	/// <param name="type">Type.</param>
	/// <param name="filename">Filename.</param>
	public GtfExpressionFile(ExpressionType type, string filename)
		: base(filename, new Layout 
			{
				Chromosome = 0,
				Start = 3,
				End = 4,
				Strand = 6,
				Name = 8
			},
			false)
	{
		this.expressionType = type;

		ParseFile();
	}

	/// <summary>
	/// Queries the interface.
	/// </summary>
	/// <returns>The interface.</returns>
	/// <param name="t">T.</param>
	override protected object QueryInterface(Type t)
	{
		if (t == typeof(IExpressionData))
		{
			return new ExpressionDataProxy(this) as IExpressionData;
		}

		return base.QueryInterface(t);
	}

	/// <summary>
	/// Gets a value indicating whether this instance is valid transcript name.
	/// </summary>
	/// <value><c>true</c> if this instance is valid transcript name; otherwise, <c>false</c>.</value>
	protected Func<string, bool> IsValidTranscriptName
	{
		get
		{
			return Helpers.CheckInit(
				ref this.isValidTranscriptName,
				() =>
				{
					if (this.Annotation == null)
					{
						return transcriptName => true;
					}
					else
					{
						return transcriptName => this.Annotation.ValidTranscripts.Contains(transcriptName);
					}
				});
		}
	}

	/// <summary>
	/// Gets or sets the annotation.
	/// </summary>
	/// <value>The annotation.</value>
		protected IAnnotation Annotation { get; set; }

        /// <summary>
        /// Gets the type of the expression.
        /// </summary>
        /// <returns>The expression type.</returns>
        /// <param name="type">Type.</param>
        public static ExpressionType ExpressionTypeFromString(string type)
        {
            return (GtfExpressionFile.ExpressionType)Enum.Parse(
                typeof(GtfExpressionFile.ExpressionType),
                type);
        }

        /// <summary>
        /// Gets the transcripts.
        /// </summary>
        /// <value>The transcripts.</value>
        public Dictionary<string, Location> Transcripts
        {
            get
            {
                return this.Locations;
            }
        }

        private Dictionary<string, Location> genes;

        public Dictionary<string, Location> Genes
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.genes,
                    () => this.Transcripts
                        .ToLookup(x => x.Value.AlternateName, x => x.Value)
                        .ToDictionary(x => x.Key, x => 
			  {
			      var gene = x.First();
			      gene.Score = x.Sum(y => y.Score);
			      return gene;
			  }));
            }
        }

        /// <summary>
        /// Extracts transcript name and FPKM/RPM values and creates bed locations with score
        /// When multiple transcripts are present, adds all transcripts
        /// </summary>
        /// <param name="fields">Fields.</param>
        /// <param name="layout">Layout.</param>
        /// <param name="data">Data.</param>
        /// <param name="entryCount">Entry count.</param>
        protected override void ParseFields(string[] fields, Layout layout, List<Tuple<Genomics.Location, string>> data, ref int entryCount)
        {
            //Console.WriteLine(EntryTypeFieldFormat [expressionType] + "\t" + expressionType + "\t" + fields [2]);

            if (EntryTypeFieldFormat [expressionType] == fields [2])
            {
		string geneName = "";
		List<string> transcriptNames = ParseNameData(fields [layout.Name], ref entryCount, out geneName);

                var startSite = new Genomics.Location
		    {
			Start = int.Parse(fields [layout.Start]),
			End = int.Parse(fields [layout.End])
		    };

                foreach (string transcriptName in transcriptNames)
                {
                    string name = transcriptName;

                    var location = new Genomics.Location
                    {
                        Name = name,
                        Chromosome = fields [layout.Chromosome], 
                        Start = int.Parse(fields [layout.Start]),
                        End = int.Parse(fields [layout.End]),
                        Strand = fields [layout.Strand],
                        Score = expressionLevel,
			Data = fields,
                        AlternateName = geneName
                    };

		    // Short expression only measures exons, so use each exon
		    // start as it's own 'transcript' and 'gene'
                    if (expressionType == ExpressionType.Short)
                    {
			name = transcriptName + "." + location.DirectionalStart;

                        location = new Genomics.Location
                        {
                            Name = name,
                            Chromosome = fields [layout.Chromosome], 
                            Start = int.Parse(fields [layout.Start]),
                            End = int.Parse(fields [layout.End]),
                            Strand = fields [layout.Strand],
                            Score = expressionLevel,
			    Data = fields,
                            AlternateName = geneName
                        };
                    }

                    if (this.IsValidTranscriptName(name)) {
			 data.Add(new Tuple<Genomics.Location, string>(location, name));
		    }
                }
            }
        }

        /// <summary>
        /// Extracts the transcript ids and rpm/fpkm values
        /// </summary>
        /// <returns>The name data.</returns>
        /// <param name="name">Name.</param>
        /// <param name="i">The index.</param>
		protected List<string> ParseNameData(string name, ref int i, out string geneName)
        {
            var data = name.Split(';').Where(x => !string.IsNullOrWhiteSpace(x))
                .Select(x =>
                {
                    var fields = x.Trim().Split(' ');

                    return new 
                    {
                        Key = fields [0],
                        Value = fields [1].Split('"') [1].Split(',').Select(y => y.Trim()).Where(y => !string.IsNullOrEmpty(y)).ToList()
                    };
                })
                .ToDictionary(x => x.Key, x => x.Value);

            //string genename = data ["gene_id"].First().Split('.') [0];
            List<string> transcriptNames = data[TranscriptFieldFormat[expressionType]].Select(x => x.Split('.') [0]).ToList();

            if (this.expressionType == ExpressionType.Cage && transcriptNames.Count > 1)
            {
                transcriptNames.RemoveRange(1, transcriptNames.Count - 1);
            }

	    geneName = data["gene_id"].Select(x => x.Split('.') [0]).First();

            var expressionData = data
                .Where(x => ExpressionFieldFormat[expressionType].IsMatch(x.Key))
                .ToList();

            // Get the average of the two expression replicates
            expressionLevel = expressionData.Count > 0 ? expressionData.Average(x => double.Parse(x.Value.First())) : double.NaN;

            //if (expressionLevel > 0)
            //{
            //    Console.WriteLine("Name field {0}; expression level : {1}; transcript {2}", name, expressionLevel, string.Join(",", transcriptNames));
            //}

            return transcriptNames;
        }

		private class ExpressionDataProxy : IExpressionData
		{
			private GtfExpressionFile data;

			public ExpressionDataProxy(GtfExpressionFile data)
			{
				this.data = data;
			}

            override public Dictionary<string, Location> Transcripts { get { return this.data.Transcripts; } }
		}
    }
}

