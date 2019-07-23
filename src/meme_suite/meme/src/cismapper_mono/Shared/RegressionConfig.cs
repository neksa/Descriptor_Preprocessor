//--------------------------------------------------------------------------------
// File: Config.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Shared
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using System.Text.RegularExpressions;
    using System.Xml.Linq;

	/// <summary>
	/// Experiment XML configuration file wrapper
	/// </summary>
    public class RegressionConfig : ConfigBase
    {
        public enum DataSourceType
        {
            Histone,
            TRF,
            AllChIP,
            SingleTF,
        }

		/// <summary>
		/// The data sources.
		/// </summary>
        public static Dictionary<DataSourceType, string> DataSources = new Dictionary<DataSourceType, string>
        {
            { DataSourceType.Histone, "HistoneListFileName" },
            { DataSourceType.AllChIP, "AllChIPListFileName" },
            { DataSourceType.TRF, "TRFListFileName" },
            { DataSourceType.SingleTF, "TFSingeChipListFileName" },
        };

        public static DataSourceType DataSourceTypeFromString(string dataSource)
        {
            return (DataSourceType)Enum.Parse(typeof(DataSourceType), dataSource);
        }

        /// <summary>
        /// The promoter range.
        /// </summary>
        private List<int> promoterRange;

		/// <summary>
        /// Initializes a new instance of the <see cref="Shared.RegressionConfig"/> class.
		/// </summary>
		/// <param name="xmlFile">Xml file.</param>
        public RegressionConfig(string xmlFile) : base(xmlFile)
        {
        }

		/// <summary>
		/// Gets the chip data source file list.
		/// </summary>
		/// <returns>The chip data source file list.</returns>
		/// <param name="dataSource">Data source.</param>
        public string GetChipDataSourceFileList(DataSourceType dataSource)
        {
            return Document.Descendants(DataSources[dataSource]).First().Value.ToString();
        }

		/// <summary>
		/// Gets the chip data source file list.
		/// </summary>
		/// <returns>The chip data source file list.</returns>
		/// <param name="dataSource">Data source.</param>
		/// <param name="tissue">Cell line.</param>
        public string GetChipDataSourceFileList(DataSourceType dataSource, string tissue)
        {
            return this.TissueTag.Replace(Document.Descendants(DataSources[dataSource]).First().Value.ToString(), tissue);
        }

        /// <summary>
        /// Gets the name of the tss file.
        /// </summary>
        /// <value>The name of the tss file.</value>
        public string TssFileName
        {
            get {  return this.StringValues ["TssFileName"]; }
        }

        /// <summary>
        /// Gets the name of the Locus file.
        /// </summary>
        /// <value>The name of the Locus file.</value>
        public string LocusFileName
        {
            get {  return this.StringValues["LocusFileName"]; }
        }

        /// <summary>
        /// Gets the Locus file format.
        /// </summary>
        /// <value>The Locus file format.</value>
        public string LocusFileFormat
        {
            get { return this.StringValues ["LocusFileFormat"]; }
        }

        /// <summary>
        /// Gets the Locus file name for cell line.
        /// </summary>
        /// <returns>The Locus file name for cell line.</returns>
        /// <param name="tissue">Cell line.</param>
        public string GetLocusFileNameForParam(string param)
        {
            return string.Format(this.LocusFileFormat, param);
        }

        /// <summary>
        /// Gets the map file format.
        /// </summary>
        /// <value>The map file format.</value>
        public string MapFileFormat
        {
            get
            {
                return this.StringValues.ContainsKey("MapFileFormat") ?
                    this.StringValues ["MapFileFormat"] :
                    null;
            }
        }

        /// <summary>
        /// Gets the expression file format.
        /// </summary>
        /// <value>The expression file format.</value>
        public string ExpressionFileFormat
        {
            get {  return this.StringValues ["ExpressionFileFormat"]; }
        }

        /// <summary>
        /// Gets the name of the expression file.
        /// </summary>
        /// <returns>The expression file name.</returns>
        /// <param name="rnaSource">Rna source.</param>
        /// <param name="tissue">Cell line.</param>
        public string GetExpressionFileName(string rnaSource, string tissue)
        {
            return this.RnaSourceTag.Replace(this.TissueTag.Replace(ExpressionFileFormat, tissue), rnaSource);
        }

        /// <summary>
        /// Gets the out file format.
        /// </summary>
        /// <value>The out file format.</value>
        public string OutFileFormat
        {
            get { return this.StringValues["OutputFileFormat"]; }
        }

        /// <summary>
        /// Gets the genome wide expression file.
        /// </summary>
        /// <returns>The genome wide expression file.</returns>
        /// <param name="tissue">Cell line.</param>
        public string GetGenomeWideExpressionFile(string tissue)
        {
            return string.Format(this.StringValues["GenomeWideRegressionFile"], tissue);
        }

        /// <summary>
        /// Gets the d0 decay halflife for exponential decay weight
        /// </summary>
        /// <value>The d0.</value>
        public double D0
        {
            get { return double.Parse(this.StringValues["D0"]); }
        }

        /// <summary>
        /// Maximimum range of the promoter region. 
        /// If the range has two values, then it is the upstream and downstream ranges, respectively
        /// </summary>
        /// <value>The range.</value>
        public List<int> PromoterRange
        {
            get 
            {
                return Helpers.CheckInit(
                    ref this.promoterRange,
                    () =>
                    {
                        var range = this.CsvIntValues["PromoterRange"];
                        if (range.Count > 2)
                        {
                            throw new Exception("Invalid promoter range: " + string.Join(", ", range));
                        }

                        return range;
                    });
            }
        }

        /// <summary>
        /// Gets the rna-seq type source count.
        /// </summary>
        /// <value>The source count.</value>
        public int SourceCount
        {
            get { return this.CsvValues["RnaSources"].Count(); }
        }
            
        /// <summary>
        /// Gets the cell lines.
        /// </summary>
        /// <value>The cell lines.</value>
        public string[] Tissues
        {
            get
            {
                if (this.Tissues == null && this.CsvValues.ContainsKey("Tissues"))
                {
                    this.tissues = this.CsvValues ["Tissues"].ToArray();
                }
                return Tissues;
            }
        }
         
        /// <summary>
        /// Gets the rna sources from the CSV delimited field
        /// </summary>
        /// <value>The rna sources.</value>
        public string[] RnaSources
        {
            get
            {
                if (rnaSources == null)
                {
                    var rnaSourceData = Document.Descendants("RnaSources").FirstOrDefault();
                    rnaSources = rnaSourceData == null ? null : rnaSourceData.Value.Split(',').ToArray();
                }
                return rnaSources;
            }
        }

        /// <summary>
        /// Gets the rna source. If multiple sources are listed, only the first is given
        /// </summary>
        /// <value>The rna source.</value>
        public string RnaSource
        {
            get
            {
                return RnaSources [0];
            }
        }

		/// <summary>
		/// Gets the thresholds.
		/// </summary>
		/// <value>The thresholds.</value>
        public string[] Thresholds
        {
            get
            {
                if (thresholds == null)
                {
                    var thresholdData = Document.Descendants("SourceThreshholds").FirstOrDefault();
                    thresholds = thresholdData == null ? null : thresholdData.Value.Split(',').ToArray();
                }
                return thresholds;
            }
        }

		/// <summary>
		/// Gets the cell line tag.
		/// </summary>
		/// <value>The cell line tag.</value>
        public Regex TissueTag
        {
            get
            { 
                if (tissueTag == null)
                {
					var v = Document.Descendants("TissueTag").FirstOrDefault();
					if (!string.IsNullOrEmpty((string)v))
					{
						tissueTag = new Regex(v.Value); 
					}
					else
					{
						return null;
					}
                }
                return tissueTag;
            }
        }

		/// <summary>
		/// Gets the histone tag.
		/// </summary>
		/// <value>The histone tag.</value>
        public Regex HistoneTag
        {
            get
            {
                return new Regex(Document.Descendants("HistoneTag").First().Value);
            }
        }

		/// <summary>
		/// Gets the threshold tag.
		/// </summary>
		/// <value>The threshold tag.</value>
        public Regex ThresholdTag
        {
            get
            { 
                if (thresholdTag == null)
                {
                    thresholdTag = new Regex(Document.Descendants("ThreshholdTag").First().Value); 
                }
                return thresholdTag;
            }
        }

		/// <summary>
		/// Gets the rna source tag.
		/// </summary>
		/// <value>The rna source tag.</value>
        public Regex RnaSourceTag
        {
            get
            { 
                if (rnaSourceTag == null)
                {
					var v = Document.Descendants("RnaTag").FirstOrDefault();
					if (!string.IsNullOrEmpty((string)v))
					{
						rnaSourceTag = new Regex(v.Value); 
					}
					else
					{
						return null;
					}
                }
                return rnaSourceTag;
            }
        }

        /// <summary>
        /// Gets the map file names for cell line and all rna sources
        /// </summary>
        /// <returns>The map file names.</returns>
        /// <param name="tissue">Cell line.</param>
        public string[][] GetMapFileNamesWithThreshold(string tissue)
        {
            return RnaSources.Select(rnaSource => 
                       Thresholds.Select(threshhold =>
                          TissueTag.Replace(
                              RnaSourceTag.Replace(
                                  ThresholdTag.Replace(MapFileFormat, threshhold), 
                                  rnaSource),  
                              tissue))
                           .ToArray())
                       .ToArray();
        }

        /// <summary>
        /// Gets the map file names for cell line a single rna sources
        /// without any threshold selection
        /// </summary>
        /// <returns>The map file names.</returns>
        /// <param name="tissue">Cell line.</param>
        public string GetMapFileName(string tissue, string rnaSource, string histoneName)
        {
            return HistoneTag.Replace(TissueTag.Replace(RnaSourceTag.Replace(StringValues["MapFileName"], rnaSource), tissue), histoneName);
        }

		/// <summary>
		/// Gets the map file name for transcript types.
		/// </summary>
		/// <returns>The map file name for transcript types.</returns>
		/// <param name="tissue">Cell line.</param>
		/// <param name="rnaSource">Rna source.</param>
		/// <param name="histoneName">Histone name.</param>
		/// <param name="transcriptTypeString">Transcript type string.</param>
		public string GetMapFileNameForTranscriptTypes(string tissue, string rnaSource, string histoneName, string transcriptTypeString)
		{
			string mapFileName = this.GetMapFileName(tissue, rnaSource, histoneName);
			if (string.IsNullOrEmpty(transcriptTypeString))
			{
				if (this.StringValues.ContainsKey("OmitLociource"))
				{
					mapFileName = Path.Combine(Path.GetDirectoryName(mapFileName), "Full." + Path.GetFileName(mapFileName));
				}
				return mapFileName;
			}
				
			var transcriptTypes = transcriptTypeString.Split(',');
			return string.Format(Path.Combine(
				Path.GetDirectoryName(mapFileName), "{0}{1}.Correlation.{2}.{3}.{4}.bed"), 
				this.StringValues.ContainsKey("OmitLociource") ? string.Empty : "Full.",
				tissue,
				rnaSource, 
				histoneName,
				string.Join("-", transcriptTypes.OrderBy(x => x)));
		}

		/// <summary>
		/// Gets the map file name for omitted cell lines and transcript types.
		/// </summary>
		/// <returns>The map file name for omitted cell lines and transcript types.</returns>
		/// <param name="tissue">Cell line.</param>
		/// <param name="otherTissue">Other cell line.</param>
		/// <param name="rnaSource">Rna source.</param>
		/// <param name="histoneName">Histone name.</param>
		/// <param name="transcriptTypeString">Transcript type string.</param>
		public string GetMapFileNameForOmittedTissuesAndTranscriptTypes(string tissue, string otherTissue, string rnaSource, string histoneName, string transcriptTypeString)
		{
			string mapFileName = this.GetMapFileName(tissue, rnaSource, histoneName);
			if (string.IsNullOrEmpty(transcriptTypeString))
			{
				if (this.StringValues.ContainsKey("OmitLociource"))
				{
					mapFileName = Path.Combine(Path.GetDirectoryName(mapFileName), "Full." + Path.GetFileName(mapFileName));
				}
				return mapFileName;
			}

			var transcriptTypes = transcriptTypeString.Split(',');
			return string.Format(Path.Combine(
                Path.GetDirectoryName(mapFileName), "{0}{1}.Correlation.{2}.{3}.{4}.{5}bed"), 
				"Sans" + tissue + (otherTissue != null ? "_" + otherTissue : string.Empty) + ".",
				tissue,
				rnaSource, 
				histoneName,
                string.Join("-", transcriptTypes.OrderBy(x => x)),
                this.StringValues.ContainsKey("MaxRange") ? this.StringValues["MaxRange"] + "." : string.Empty);
		}

        /// <summary>
        /// Gets the table file names.
        /// </summary>
        /// <returns>The table file names.</returns>
        /// <param name="tissue">Cell line.</param>
        public string[][] GetTableFileNames(string tissue)
        {
            return RnaSources.Select(source =>
                       Thresholds.Select(threshold => 
                           tissueTag.Replace(ThresholdTag.Replace(RnaSourceTag.Replace(
                               OutFileFormat, 
                               source),
                           threshold), 
                       tissue)).ToArray()).ToArray();
        }

        //private XDocument Document;

        private string[] tissues;
        private string[] rnaSources;
        private string[] thresholds;

        private Regex tissueTag;
        private Regex thresholdTag;
        private Regex rnaSourceTag;
    }
}

