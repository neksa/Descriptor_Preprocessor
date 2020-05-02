// <copyright file="LocusFileElement.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using Genomics;
    using Shared;

    /// <summary>
    /// Locus file element.
    /// </summary>
    public class LocusFileElement : AnalysisElement<BedFile>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.LocusFileElement"/> class.
        /// </summary>
        /// <param name="propertyPrefix">Property prefix.</param>
        public LocusFileElement(string propertyPrefix)
            : base(propertyPrefix)
        {
        }

        /// <summary>
        /// Gets or sets the name of the Locus file.
        /// </summary>
        /// <value>The name of the Locus file.</value>
        public string LocusFileName { get; set; }

        /// <summary>
        /// Gets the element.
        /// </summary>
        /// <value>The element.</value>
        public override BedFile Element
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.element,
                    () => new BedFile(this.LocusFileName, BedFile.Bed3Layout));
            }
        }

        /// <summary>
        /// Register the specified analysis.
        /// </summary>
        /// <param name="analysis">Analysis.</param>
        public override void Register(BaseAnalysis analysis)
        {
            analysis.ElementArgRegistry.Add("LocusFileName", (val) => this.LocusFileName = val);
        }
    }
}

