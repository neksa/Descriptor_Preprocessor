//--------------------------------------------------------------------------------
// <copyright file="PredictElements.cs" 
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
    using Shared;

    /// <summary>
    /// Predict elements.
    /// </summary>
    public class PredictElements : PredictionBase
    {
        /// <summary>
        /// Sets the element arguments.
        /// </summary>
        /// <value>The element arguments.</value>
        public enum ElementArgs
        {
            /// <summary>
            /// The output file.
            /// </summary>
            OutputFile,

            /// <summary>
            /// The name of the map file.
            /// </summary>
            MapFileName,

            /// <summary>
            /// The use genes.
            /// </summary>
            UseGenes,

            /// <summary>
            /// The threshold.
            /// </summary>
            Threshold,

            /// <summary>
            /// The type of the threshold.
            /// </summary>
            ThresholdType,

            /// <summary>
            /// The target gene(s)
            /// </summary>
            Targets,
        }


        /// <summary>
        /// Predict this instance.
        /// </summary>
        public override void Predict()
        {
            var targets = this.Targets != null ? 
                this.Targets
                    .Where(x => this.Map.ContainsKey(x))
                    .OrderBy(x => this.Map[x].Values.Min(y => y.ConfidenceScore))
                    .ToArray() :
                this.Map.Keys
                    .OrderBy(x => this.Map[x].Values.Min(y => y.ConfidenceScore))
                    .Select(x => (string)x)
                    .ToArray();

            var elementLists = targets
                .Select(t => this.Map[t].Values
                  //.OrderBy(x => x.Correlation > 0 ? 0 : 1)
                  .OrderBy(x => x.ConfidenceScore)
                  .ThenBy(x => x.LocusName.Chr)
                  .ThenBy(x => x.LocusStart)
                  .ThenBy(x => x.LocusEnd)
                  .ToList());
            
            Tables.ToNamedTsvFile(
                this.OutputFile,
                elementLists.Select(x => x.Take(this.ThresholdType == ThresholdTypes.Rank ? 
                    Math.Min((int)this.Threshold, x.Count) :
                    x.Count))
                .SelectMany(x => x)
                .Select(x => new string[]
                {
                    x.TssName,
                    x.LocusName,
		    x.Correlation < 0 ? "-" : "+",
		    x.ConfidenceScore.ToString(x.ConfidenceScore == 0 ? "" : "0.000e+000"),
                }),
                new string[] { UseGenes ? "Gene" : "Tss", "Locus", "Correlation_Sign", "Score", });
        }

        /// <summary>
        /// Gets or sets the targets.
        /// </summary>
        /// <value>The targets.</value>
        public string[] Targets { get; set; }

        /// <summary>
        /// Executes PredictLinks.
        /// </summary>
        public class Executor : PredictionBase.ExecutorBase<PredictElements, ElementArgs>
        {
            /// <summary>
            /// Gets the description of the task.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Predicts regulatory elements most likely to regulate a given list of targets.";
                }
            }

            /// <summary>
            /// Gets the command line option descriptions.
            /// </summary>
            /// <value>The options data.</value>
            protected override Dictionary<ElementArgs, string> OptionsData
            {
                get
                {
                    return this.BasicOptionsData.Select(x => new Tuple<ElementArgs, string>(
                                Helpers.ConvertEnum<BasicArgs, ElementArgs>(x.Key),
                                x.Value))
                        .Concat(new List<Tuple<ElementArgs, string>>
                        {
                            new Tuple<ElementArgs, string>(
                                ElementArgs.Targets, 
                                "TSS/Gene targets for which to predict regulating elements, listed in CSV format")
                        })
                        .ToDictionary(x => x.Item1, x => x.Item2);
                }
            }

            /// <summary>
            /// Reflects the arguments.
            /// </summary>
            public override void ReflectArgs(PredictElements predictor)
            {
                if (this.CommandArgs.StringEnumArgs.ContainsKey(ElementArgs.Targets))
                {
                    predictor.Targets = this.CommandArgs.StringEnumArgs[ElementArgs.Targets].Split(',');
                }
            }
        }
    }
}

