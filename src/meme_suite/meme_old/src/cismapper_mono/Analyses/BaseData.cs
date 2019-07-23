//--------------------------------------------------------------------------------
// File: BaseData.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using Genomics;
    using Shared;

    /// <summary>
    /// Base data for map experiments using RNA measurements, histone measurements, and loci with regulatory potential in a particular cell line.
    /// </summary>
	public class BaseData : BaseAnalysis
    {
        /// <summary>
        /// The cell line.
        /// </summary>
        private string tissue;

        /// <summary>
        /// The rna source.
        /// </summary>
        private string rnaSource;

        /// <summary>
        /// The locus file.
        /// </summary>
        private BedFile locusFile;

        /// <summary>
        /// The name of the bed file.
        /// </summary>
        private string locusFileName;

        /// <summary>
        /// Histone used to build the map
        /// </summary>
        private String histoneName;
            
        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.BaseData"/> class.
        /// </summary>
        /// <param name="xmlFile">Xml file.</param>
        public BaseData(string outputRoot, string plotRoot, string scriptRoot, string xmlFile)
            : base(outputRoot, plotRoot, scriptRoot)
        {
            ConfigData = new RegressionConfig(xmlFile);
        }

		/// <summary>
		/// The config data.
		/// </summary>
		public RegressionConfig ConfigData { get; private set; }

        public bool UseHistoneData
        {
            get
            {
                return !string.IsNullOrEmpty(this.histoneName);
            }
        }

        /// <summary>
        /// Gets or sets the name of the histone used to build the map
        /// </summary>
        /// <value>The name of the histone.</value>
        public string HistoneName
        {
            get
            {
                return !this.UseHistoneData ? "None" : this.histoneName;
            }

            set
            {
                this.histoneName = value;
            }
        }

        /// <summary>
        /// Gets or sets the cell line.
        /// </summary>
        /// <value>The cell line.</value>
        public string Tissue
        {
            get
            {
                return this.tissue;
            }

            set
            {
                this.tissue = value;
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
                return rnaSource;
            }

            set
            {
                rnaSource = value;
            }
        }

        /// <summary>
        /// Gets or sets the name of the Locus file.
        /// </summary>
        /// <value>The name of the Locus file.</value>
		public string LocusFileName
        {
            get
            {
                if (this.locusFileName == null)
                {
                    return ConfigData.GetLocusFileNameForParam(this.LocusFileNameParam);
                }

                return this.locusFileName;
            }

            set
            {
                this.locusFileName = value;
            }
        }

        /// <summary>
        /// Gets the Locus file name parameter.
        /// </summary>
        /// <value>The Locus file name parameter.</value>
        public virtual string LocusFileNameParam
        {
            get
            {  
                return this.Tissue;
            }
        }

        /// <summary>
        /// Locus data 
        /// </summary>
        /// <value>The Locus file.</value>
		public BedFile LocusFile
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.locusFile,
                    () => new BedFile(this.LocusFileName, BedFile.Bed3Layout));
            }
        }
    }
}

    
