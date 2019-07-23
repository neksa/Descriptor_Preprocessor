//--------------------------------------------------------------------------------
// File: GtfGencodeFile.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;

    public class GtfGencodeGeneFile : GtfGencodeFile
    {
        override protected string Feature
        {
            get
            {
                return "gene";
            }
        }

        public GtfGencodeGeneFile(string filename)
            : base(filename)
        {
        }
    }

    public class GtfGencodeTranscriptFile : GtfGencodeFile
    {
        override protected string Feature
        {
            get
            {
                return "transcript";
            }
        }

        public GtfGencodeTranscriptFile(string filename)
            : base(filename)
        {
        }
    }

    /// <summary>
    /// Gtf gencode file.
    /// </summary>
	public abstract class GtfGencodeFile : BedFile
    {
        public enum TranscriptType
        {
            protein_coding = 0x01,
            processed_transcript = 0x02,
        };

	/// <summary>
	/// The gencode files.
	/// </summary>
	private static readonly Dictionary<string, GtfGencodeFile> GencodeFiles = new Dictionary<string, GtfGencodeFile>();

        /// <summary>
        /// Cache of valid transcript sets. Key is filename.TranscriptType.TranscriptType...
        /// </summary>
        private static readonly Dictionary<string, HashSet<string>> staticValidTranscripts = new Dictionary<string, HashSet<string>>();

        /// <summary>
        /// The static valid genes.
        /// </summary>
        private static readonly Dictionary<string, HashSet<string>> staticValidGenes = new Dictionary<string, HashSet<string>>();

        /// <summary>
        /// Cache of  ordered gene locations.
        /// </summary>
        private static readonly Dictionary<string, Dictionary<string, List<Genomics.Location>>> staticOrderedTranscriptLocations = new Dictionary<string, Dictionary<string, List<Genomics.Location>>>();

        /// <summary>
        /// The valid types backint the ValidGeneTypes property
        /// </summary>
	private TranscriptType[] validTranscriptTypes;

        /// <summary>
        /// Filename of this instance
        /// </summary>
        private readonly string filename;

        /// <summary>
        /// The feature type to use
        /// </summary>
        abstract protected string Feature { get; }

	/// <summary>
        /// Load the specified filename and feature.
        /// </summary>
        /// <param name="filename">Filename.</param>
        /// <param name="feature">Feature.</param>
        public static GtfGencodeFile LoadGenes(string filename)
	{
            if (!GencodeFiles.ContainsKey(filename + "_gene"))
	    {
                GencodeFiles.Add(filename, new GtfGencodeGeneFile(filename));
	    }

	    return GencodeFiles [filename];
	}

        public static GtfGencodeFile LoadTranscripts(string filename)
        {
            if (!GencodeFiles.ContainsKey(filename + "_transcript"))
            {
                GencodeFiles.Add(filename, new GtfGencodeTranscriptFile(filename));
            }

            return GencodeFiles [filename];
        }

	/// <summary>
	/// Initializes a new instance of the <see cref="Genomics.GtfGencodeFile"/> class.
	/// </summary>
	/// <param name="filename">Filename.</param>
	protected GtfGencodeFile(string filename)
            : base (filename, new Layout 
	     {
		  Chromosome = 0,
		  Start = 3,
		  End = 4,
		  Strand = 6,
		  Name = 8
             }
          )
	  {
	      this.filename = filename;
	  }

	/// <summary>
	/// Queries the interface.
	/// </summary>
	/// <returns>The interface.</returns>
	/// <param name="t">T.</param>
	override protected object QueryInterface(Type t)
	{
		return t == typeof(IAnnotation) ? new AnnotationProxy(this) as IAnnotation : base.QueryInterface(t);
	}

        /// <summary>
        /// Gets or sets the valid gene types.
        /// </summary>
        /// <value>The valid gene types.</value>
        public string[] ValidTranscriptTypes
        {
            get
            {
				return this.validTranscriptTypes.Select(x => x.ToString()).ToArray();
            }

            set
            {
				if (value.Count(x => !string.IsNullOrWhiteSpace(x)) == 0)
				{
					this.validTranscriptTypes = new TranscriptType[]{ };
				}
				else
				{
					this.validTranscriptTypes = value.Select(x => (TranscriptType)Enum.Parse(typeof(TranscriptType), x)).ToArray();
				}
            }
        }

        private string InitValidFeature(string feature, Dictionary<string, HashSet<string>> validFeatures)
        {
            string validGeneKey = GenerateKey(this.validTranscriptTypes);

            if (!validFeatures.ContainsKey(validGeneKey))
            {
                var featureSet = this.Locations.Where(x => x.Value.Data [2] == feature);

                Func<string, bool> TestType = l => this.ValidTranscriptTypes.Contains(l);

                // All transcript types are valid
                if (this.validTranscriptTypes.Length == 0)
                {
                    TestType = l => true;
                }

                validFeatures.Add(validGeneKey, new HashSet<string>(featureSet.Where(x => TestType(x.Value.AdditionalFields ["transcript_type"])).Select(x => x.Value.Name)));
            }

            return validGeneKey;
        }

        /// <summary>
        /// Gets the valid transcripts.
        /// </summary>
        /// <returns>The valid transcripts.</returns>
        public HashSet<string> ValidTranscripts
        {
            get
            {
                string validTranscriptKey = this.InitValidFeature("transcript", staticValidTranscripts);

                return staticValidTranscripts [validTranscriptKey];
            }
        }

        /// <summary>
        /// Gets the valid genes.
        /// </summary>
        /// <value>The valid genes.</value>
        public HashSet<string> ValidGenes
        {
            get
            {
                string validGeneKey = this.InitValidFeature("gene", staticValidGenes);

                return staticValidGenes [validGeneKey];
            }
        }

        /// <summary>
        /// Gets the ordered gene locations by chromosome and then position
        /// </summary>
        /// <returns>The ordered gene locations.</returns>
        public Dictionary<string, List<Genomics.Location>> OrderedTranscriptLocations
        {
            get
            {
		var key = GenerateKey(this.validTranscriptTypes);

                if (!staticOrderedTranscriptLocations.ContainsKey(key))
                {
                    staticOrderedTranscriptLocations.Add(key, this.ChromosomeOrderedLocations.ToDictionary(x => x.Key, x => x.Value.Where(y => this.ValidTranscripts.Contains(y.Name)).ToList()));
                }

                return staticOrderedTranscriptLocations [key];
            }
        }

        /// <summary>
        /// Gets maximum feature size
        /// </summary>
        /// <value>The size of the max feature.</value>
        public int MaxFeatureSize
        {
            get
            {
                return this.MaxLocationSize.Max(x => x.Value);
            }
        }

        /// <summary>
        /// Gets the indexed transcript locations.
        /// </summary>
        /// <value>The indexed transcript locations.</value>
        public Dictionary<string, Dictionary<int, List<Genomics.Location>>> IndexedTranscriptLocations
        {
            get
            {
                return this.ChromosomeIndexedLocations;
            }
        }

	/// <summary>
	/// Parses the fields.
	/// </summary>
	/// <param name="fields">Fields.</param>
	/// <param name="layout">Layout.</param>
	/// <param name="data">Data.</param>
	/// <param name="entryCount">Entry count.</param>
        protected override void ParseFields(string[] fields, Layout layout, List<Tuple<Genomics.Location, string>> data, ref int entryCount) {
            if (fields[0][0] == '#') { return; }

            if (fields [2] != this.Feature) { return; }

  	    var nameData = ParseNameData(fields [layout.Name]);

            var startSite = new Genomics.Location
	    {
		Start = int.Parse(fields [layout.Start]),
		End = int.Parse(fields [layout.End])
	    };

	    linecount++;
	    if (linecount % 10000 == 0)
	    {
		Console.WriteLine(linecount);
	    }

	    var name = nameData["transcript_id"];
	    string geneName = nameData["gene_id"];

            var location = new Genomics.Location
	    {
		Name = nameData["transcript_id"],
		Chromosome = fields [layout.Chromosome], 
		Start = int.Parse(fields [layout.Start]),
		End = int.Parse(fields [layout.End]),
		Strand = fields [layout.Strand],
		Data = fields,
		AdditionalFields = nameData,
		AlternateName = geneName + "." + startSite.DirectionalStart
	    };

	    data.Add(new Tuple<Genomics.Location, string>(
		location,
		nameData["transcript_id"])
	    );
	}
        
        private int linecount = 0;

        /// <summary>
        /// Parses the name data producing a dictionary of all named data element lists
        /// </summary>
        /// <returns>The name data.</returns>
        /// <param name="name">Name.</param>
	protected Dictionary<string, string> ParseNameData(string name)
        {
            return  name.Split(';').Where(x => !string.IsNullOrWhiteSpace(x))
	        .Select(x =>
		    {
			var fields = x.Trim().Split(' ');
			var value = (fields[0] == "level" ? fields[1] : fields[1].Replace("\"", "")).Trim();

			if (fields[0] == "gene_id" || fields[0] == "transcript_id")
			{
			    value = value.Split('.')[0].Trim();
			}

			return new 
			{
			    Key = fields[0],
			    Value = value
			};
		      }
                )

		.ToLookup(x => x.Key, x => x.Value)
		.ToDictionary(x => x.Key, x => x.First()); // Only unused fields 'ont' and 'tag' are multiply used per transcript
        }

        /// <summary>
        /// Generates keys for cached data
        /// </summary>
        /// <returns>The key.</returns>
        /// <param name="types">Types.</param>
        string GenerateKey(TranscriptType[] types) {
	    if (types == null)
	    {
		throw new Exception("Transcript types not set in annotation implementation " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
	    }
            return this.filename + "." + string.Join(".", types.OrderBy(x => x).Select(x => x.ToString()));
        }

	private class AnnotationProxy : IAnnotation
	{
	    private readonly GtfGencodeFile data;
	    public AnnotationProxy(GtfGencodeFile data)
	    {
		this.data = data;
	    }

	    override public string[] ValidTranscriptTypes { get { return this.data.ValidTranscriptTypes; } set { this.data.ValidTranscriptTypes = value; } }
            override public HashSet<string> ValidGenes { get { return this.data.ValidGenes; } }
	    override public HashSet<string> ValidTranscripts { get { return this.data.ValidTranscripts; } }
	    override public Dictionary<string, List<Genomics.Location>> OrderedTranscriptLocations { get { return this.data.OrderedTranscriptLocations; } }
	    override public Dictionary<string, Dictionary<int, List<Genomics.Location>>> IndexedTranscriptLocations { get { return this.data.IndexedTranscriptLocations; } }
	    override public int MaxFeatureSize { get { return this.data.MaxFeatureSize; } }
	    // Analysis disable AccessToStaticMemberViaDerivedType
	    // Analysis disable MemberHidesStaticFromOuterClass
	    override public int IndexSize { get { return GtfGencodeFile.IndexSize; } }
	    // Analysis restore MemberHidesStaticFromOuterClass
	    // Analysis restore AccessToStaticMemberViaDerivedType
	}
    }
}

