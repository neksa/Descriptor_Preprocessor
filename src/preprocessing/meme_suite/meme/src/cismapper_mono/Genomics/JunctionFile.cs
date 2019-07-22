
namespace Genomics
{
    using System;
    using Shared;
    using System.Collections.Generic;
    using System.Linq;

    public class JunctionFile : BedFile
    {
        private Dictionary<string, List<string>> geneJunctions;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.JunctionFile"/> class.
        /// </summary>
        /// <param name="fileName">File name.</param>
        public JunctionFile(string fileName)
            : base(fileName, new Layout
            {
                Chromosome = 0,
                Start = 10,
                End = 11,
                Name = 3,
                Score = 4,
                Strand = 5,
            })
        {
        }

        /// <summary>
        /// Gets the gene junctions.
        /// </summary>
        /// <value>The gene junctions.</value>
        public Dictionary<string, List<string>> GeneJunctions
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.geneJunctions,
                    () => this.Locations
                        .ToLookup(x => x.Value.AlternateName, x => x.Value.Name)
                        .ToDictionary(x => x.Key, x => x.ToList()));
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
            if (fields[0].Length >= 5 && fields[0].Substring(0, 5) == "track")
            {
                return;
            }

            string geneName = fields[layout.Name].Split('_')[0];

            int readStart = int.Parse(fields[1]);

            string[] startOffsetData = fields[layout.Start].Split(',');
            string[] endOffsetData = fields[layout.End].Split(',');

            // Ensure we don't add TSS or TTS as junctions
            if (startOffsetData.Length > 1 && endOffsetData.Length > 1)
            {
                int start = readStart + int.Parse(startOffsetData[0]);
                int end   = readStart + int.Parse(endOffsetData[1]);

                string junctionName = GtfLinkingFile.JunctionName(
                    fields[layout.Chromosome],
                    start,
                    end,
                    fields[layout.Strand]);

                var location = new Location
                {
                    Name = junctionName,
                    AlternateName = geneName,
                    Chromosome = fields[layout.Chromosome],
                    Start = start,
                    End = end,
                    Strand = fields[layout.Strand],
                    Score = double.Parse(fields[layout.Score])
                };

                data.Add(new Tuple<Location, string>(location, junctionName));
                entryCount++;
            }
        }
    }
}

