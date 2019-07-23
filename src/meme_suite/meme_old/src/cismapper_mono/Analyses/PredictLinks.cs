//--------------------------------------------------------------------------------
// <copyright file="PredictLinks.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2015. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Data;

    /// <summary>
    /// Predict links.
    /// </summary>
    public class PredictLinks : PredictionBase
    {
        /// <summary>
        /// Predict this instance.
        /// </summary>
        public override void Predict()
        {
            Tables.ToNamedTsvFile(
                this.OutputFile,
                this.Map.Links
                    .OrderBy(x => x.ConfidenceScore)
                    .ThenBy(x => x.TssName)
		    .ThenBy(x => (x.Correlation > 0 ? 0 : 1))
 		    .ThenBy(x => x.LocusName.Chr)
		    .ThenBy(x => x.LocusStart)
		    .ThenBy(x => x.LocusEnd)
                    .Select(x => new string[]
                    {
                        x.TssName, 
                        x.LocusName, 
                        x.LinkLength.ToString(),
                        x.Correlation < 0 ? "-" : "+", 
			x.ConfidenceScore.ToString(x.ConfidenceScore == 0 ? "" : "0.000e+000"),
                    }),
                new string[] { UseGenes ? "Gene" : "Tss", "Locus", "Distance", "Correlation_Sign", "Score", });
        }

        /// <summary>
        /// Executes PredictLinks.
        /// </summary>
        public class Executor : PredictionBase.ExecutorBase<PredictLinks, PredictionBase.BasicArgs>
        {
            /// <summary>
            /// Gets the description of the task.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Predicts regulatory links from a map using a link score threshold.";
                }
            }

            /// <summary>
            /// Gets the command line option descriptions.
            /// </summary>
            /// <value>The options data.</value>
            protected override Dictionary<BasicArgs, string> OptionsData
            {
                get
                {
                    return this.BasicOptionsData;
                }
            }
        }
    }
}
