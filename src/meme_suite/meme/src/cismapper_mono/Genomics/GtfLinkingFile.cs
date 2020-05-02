
namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Shared;

    public class GtfLinkingFile : BedFile
    {
        /// <summary>
        /// The gene transcripts.
        /// </summary>
        private Dictionary<string, List<string>> geneTranscripts;

        /// <summary>
        /// The transcript exons.
        /// </summary>
        private Dictionary<string, Dictionary<int, string>> transcriptExons;

        /// <summary>
        /// The transcript junctions in sequential order.
        /// </summary>
        private Dictionary<string, List<string>> transcriptJunctions;

        /// <summary>
        /// The junction transcripts.
        /// </summary>
        private Dictionary<string, List<string>> junctionTranscripts;

        /// <summary>
        /// Exon locations keyed by chr:start:end:strand
        /// </summary>
        private List<string> exons;

        /// <summary>
        /// Junction locations keyed by chr:start-end
        /// </summary>
        private Dictionary<string, Location> junctionLocations;
       
        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.GtfLinkingFile"/> class.
        /// </summary>
        /// <param name="filename">Filename.</param>
        public GtfLinkingFile(string filename)
            : base(filename, new Layout
            {
                Chromosome = 0,
                Start = 3,
                End = 4,
                Strand = 6,
                Name = 8,
                Type = 2
            })
        {
        }

        /// <summary>
        /// Gets the gene transcripts.
        /// </summary>
        /// <value>The gene transcripts.</value>
        public Dictionary<string, List<string>> GeneTranscripts
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.geneTranscripts,
                    () => new Dictionary<string, List<string>>());
            }
        }

        /// <summary>
        /// Gets the transcript exons.
        /// </summary>
        /// <value>The transcript exons.</value>
        public Dictionary<string, Dictionary<int, string>> TranscriptExons
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcriptExons,
                    () => new Dictionary<string, Dictionary<int, string>>());
            }
        }

        /// <summary>
        /// Gets the exons.
        /// </summary>
        /// <value>The exons.</value>
        public List<string> Exons
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.exons,
                    () => this.TranscriptExons
                        .Select(x => x.Value.Values)
                        .SelectMany(x => x)
                        .ToLookup(x => x, x => x)
                        .Select(x => x.Key)
                        .ToList());
            }
        }

        /// <summary>
        /// Gets the transcript junctions.
        /// </summary>
        /// <value>The transcript junctions.</value>
        public Dictionary<string, List<string>> TranscriptJunctions
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcriptJunctions,
                    () => this.TranscriptExons
                    .Select(x => new 
                    {
                        Transcript = x.Key,
                        Exons = x.Value.OrderBy(y => y.Key).ToList()
                    })
                    .Select(x => new
                    {
                        Transcript = x.Transcript,
                        Junctions = x.Exons
                            .Where((y, i) => i + 1 < x.Exons.Count)
                            .Select((y, i) => this.JunctionFromExons(y.Value, x.Exons[i + 1].Value))
                            .ToList()
                    })
                    .ToDictionary(x => x.Transcript, x => x.Junctions));
            }
        }

        /// <summary>
        /// Gets the junction transcripts.
        /// </summary>
        /// <value>The junction transcripts.</value>
        public Dictionary<string, List<string>> JunctionTranscripts
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.junctionTranscripts,
                    () => this.TranscriptJunctions.Select(x => x.Value.Select(y => new
                    {
                        Transcript = x.Key,
                        Junction = y,
                    }))
                    .SelectMany(x => x)
                    .ToLookup(x => x.Junction, x => x.Transcript)
                    .ToDictionary(x => x.Key, x => x.ToList()));
            }
        }

        /// <summary>
        /// Gets the gene junctions.
        /// </summary>
        /// <value>The gene junctions.</value>
        public Dictionary<string, List<string>> GeneJunctions
        {
            get
            {
                return this.GeneTranscripts.ToDictionary(
                    x => x.Key,
                    x => x.Value
                        .Select(y => this.TranscriptJunctions[y])
                        .SelectMany(y => y)
                        .Distinct()
                        .ToList());
            }
        }

        public Dictionary<string, Location> JunctionLocations
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.junctionLocations,
                    () => this.TranscriptJunctions
                        .Select(x => x.Value.Select(y => new { Junction = y, Score = this.Locations[x.Key].Score }))
                    .SelectMany(x => x)
                    .ToLookup(x => x.Junction, x => x.Score)
                    .ToDictionary(x => x.Key, x => new Location
                    {
                        Name = x.Key,
                        Score = x.Sum()
                    }));
            }
        }

        /// <summary>
        /// Parses the fields.
        /// </summary>
        /// <param name="fields">Fields.</param>
        /// <param name="layout">Layout.</param>
        /// <param name="data">Data.</param>
        /// <param name="entryCount">Entry count.</param>
        protected override void ParseFields(string[] fields, Layout layout, List<Tuple<Genomics.Location, string>> data, ref int entryCount)
        {
            if (fields[0][0] == '#')
            {
                return;
            }

            var nameData = ParseNameData(fields [layout.Name]);

            switch (fields[layout.Type])
            {
                case "transcript":
                    {
                        string transcriptName = nameData["transcript_id"];
                        string geneName = nameData["gene_id"];
                        var trLocation = new Location
                        {
                            Name = transcriptName,
                            AlternateName = geneName,
                            Start = int.Parse(fields[layout.Start]),
                            End = int.Parse(fields[layout.End]),
                            Strand = fields[layout.Strand],
                            Score = double.Parse(nameData["FPKM"])
                        };

                        data.Add(new Tuple<Location, string>(trLocation, transcriptName));
                        entryCount++;

                        this.AddTranscriptToGene(transcriptName, geneName);
                    }
                    break;
                case "exon":
                    {
                        string exonName = ExonName(fields[layout.Chromosome], fields[layout.Start], fields[layout.End], fields[layout.Strand]);
                        string transcriptName = nameData["transcript_id"];

                        var exonLocation = new Location
                        {
                            Name          = exonName, 
                            AlternateName = transcriptName,
                            Chromosome    = fields[layout.Chromosome],
                            Start         = int.Parse(fields[layout.Start]),
                            End           = int.Parse(fields[layout.End]),
                            Strand        = fields[layout.Strand],
                        };

                        data.Add(new Tuple<Location, string>(exonLocation, exonName));
                        entryCount++;

                        this.AddExonToTranscript(transcriptName, exonName, int.Parse(nameData["exon_number"]));
                    }
                    break;
            }
        }

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

                    if (fields[0] == "gene_id" || fields[0] == "transcript_id" || fields[0] == "gene_name")
                    {
                        value = value.Split('.')[0].Trim();
                    }

                    return new 
                    {
                        Key = fields[0],
                        Value = value
                    };
                })
                .ToLookup(x => x.Key, x => x.Value)
                .ToDictionary(x => x.Key, x => x.First()); // Only unused fields 'ont' and 'tag' are multiply used per transcript
        }

        /// <summary>
        /// Creates a name for a junction based on the exons
        /// </summary>
        /// <returns>The from exons.</returns>
        /// <param name="exonName1">Exon name1.</param>
        /// <param name="exonName2">Exon name2.</param>
        private string JunctionFromExons(string exonName1, string exonName2)
        {
            return JunctionName(this.Locations[exonName1], this.Locations[exonName2]);
        }

        /// <summary>
        /// Junctions the name.
        /// </summary>
        /// <returns>The name.</returns>
        /// <param name="exon1">Exon1.</param>
        /// <param name="exon2">Exon2.</param>
        public static string JunctionName(Location exon1, Location exon2)
        {
            if (exon1.Chromosome != exon2.Chromosome || exon1.Strand != exon2.Strand)
            {
                throw new Exception("Junction crossing chromosomes or strands requested");
            }

            if (exon1.Overlaps(exon2) || exon2.Overlaps(exon1))
            {
                throw new Exception("Overlapping exons");
            }

            if (exon1.Start > exon2.Start)
            {
                var v = exon1;
                exon1 = exon2;
                exon2 = v;
            }

            //if (exon1.Strand == "-")
            //{
            //    return string.Format("{0}:{1}:{2}:{3}", exon1.Chromosome, exon2.End, exon1.Start, exon1.Strand);
            //}
            //else
            //{
            return JunctionName(exon1.Chromosome, Math.Max(exon1.Start, exon1.End), Math.Min(exon2.Start, exon2.End) - 1, exon1.Strand);
        }

        /// <summary>
        /// Creates a junction name
        /// </summary>
        /// <returns>The name.</returns>
        /// <param name="chromosome">Chromosome.</param>
        /// <param name="start">Start.</param>
        /// <param name="end">End.</param>
        /// <param name="strand">Strand.</param>
        public static string JunctionName(string chromosome, int start, int end, string strand)
        {
            return string.Format("{0}:{1}:{2}:{3}", chromosome, start, end, strand);
        }

        /// <summary>
        /// Creates a name for an exon
        /// </summary>
        /// <returns>The name.</returns>
        /// <param name="chr">Chr.</param>
        /// <param name="start">Start.</param>
        /// <param name="end">End.</param>
        /// <param name="strand">Strand.</param>
        public static string ExonName(string chr, string start, string end, string strand)
        {
            return string.Format("{0}:{1}:{2}:{3}", chr, start, end, strand);
        }

        /// <summary>
        /// Adds the transcript to the gene-to-transcript map
        /// </summary>
        /// <param name="transcriptName">Transcript name.</param>
        /// <param name="geneName">Gene name.</param>
        private void AddTranscriptToGene(string transcriptName, string geneName)
        {
            if (!this.GeneTranscripts.ContainsKey(geneName))
            {
                this.GeneTranscripts.Add(geneName, new List<string>());
            }

            this.GeneTranscripts[geneName].Add(transcriptName);
        }

        /// <summary>
        /// Adds the exon to the transcript-to-exon map
        /// </summary>
        /// <param name="transcriptName">Transcript name.</param>
        /// <param name="exonName">Exon name.</param>
        /// <param name="number">Number.</param>
        private void AddExonToTranscript(string transcriptName, string exonName, int number)
        {
            if (!this.TranscriptExons.ContainsKey(transcriptName))
            {
                this.TranscriptExons.Add(transcriptName, new Dictionary<int, string>());
            }

            this.TranscriptExons[transcriptName].Add(number, exonName);
        }
    }
}

