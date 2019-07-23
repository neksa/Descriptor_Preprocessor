//--------------------------------------------------------------------------------
// <copyright file="PredictionBase.cs" 
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
    using Genomics;
    using Shared;

    /// <summary>
    /// Prediction base.
    /// </summary>
    public abstract class PredictionBase
    {
        /// <summary>
        /// The filter.
        /// </summary>
        private MapLinkFilter filter;

        /// <summary>
        /// The map.
        /// </summary>
        private TssRegulatoryMap map;

        /// <summary>
        /// Threshold types.
        /// </summary>
        public enum ThresholdTypes
        {
            /// <summary>
            /// The score.
            /// </summary>
            Score,

            /// <summary>
            /// The rank.
            /// </summary>
            Rank,
        }

        /// <summary>
        /// Basic arguments.
        /// </summary>
        public enum BasicArgs
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
        }

        /// <summary>
        /// Gets or sets the name of the map file.
        /// </summary>
        /// <value>The name of the map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the output file.
        /// </summary>
        /// <value>The output file.</value>
        public string OutputFile { get; set; }

        /// <summary>
        /// Gets or sets the threshold.
        /// </summary>
        /// <value>The threshold.</value>
        public double Threshold { get; set; }

        /// <summary>
        /// Gets or sets the type of the threshold.
        /// </summary>
        /// <value>The type of the threshold.</value>
        public ThresholdTypes ThresholdType { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Analyses.PredictionBase"/> use genes.
        /// </summary>
        /// <value><c>true</c> if use genes; otherwise, <c>false</c>.</value>
        public bool UseGenes { get; set; }

        /// <summary>
        /// Gets the filter.
        /// </summary>
        /// <value>The filter.</value>
        public MapLinkFilter Filter
        {
            get
            {
                if (this.filter == null)
                {
                    var filter = new MapLinkFilter();
                    switch (this.ThresholdType)
                    {
                        case ThresholdTypes.Score:
                            if (this.Threshold >= 0)
                            {
                                filter.ConfidenceThreshold = this.Threshold;
                            }
                            break;
                    }

                    this.filter = filter;
                }

                return this.filter;
            }
        }

        /// <summary>
        /// Gets the map.
        /// </summary>
        /// <value>The map.</value>
        public TssRegulatoryMap Map 
        {
            get
            {
                if (this.map == null)
                {
                    var map = TssRegulatoryMap.LoadMap(
                                  this.MapFileName,
                                  this.Filter);

                    if (this.UseGenes)
                    {
                        map = map.ConvertToGenes();
                    }

                    this.map = map;
                }

                return this.map;   
            }
        }

        /// <summary>
        /// Predict this instance.
        /// </summary>
        public abstract void Predict();

        /// <summary>
        /// Executor base.
        /// </summary>
        /// <typeparam name="TPredictor">Type of predictor</typeparam>
        /// <typeparam name="TPredictionArgs">Enum of predictor arguments</typeparam>
        public abstract class ExecutorBase<TPredictor, TPredictionArgs> : IAnalysisExecutor<TPredictor, TPredictionArgs> where TPredictor : PredictionBase
        {
            /// <summary>
            /// Required arguments.
            /// </summary>
            public enum RequiredArgs
            {
                /// <summary>
                /// The output file.
                /// </summary>
                OutputFile,

                /// <summary>
                /// The name of the map file.
                /// </summary>
                MapFileName,
            }

            /// <summary>
            /// Optional arguments.
            /// </summary>
            public enum OptionalArgs
            {
                /// <summary>
                /// The threshold.
                /// </summary>
                Threshold,

                /// <summary>
                /// The type of the threshold.
                /// </summary>
                ThresholdType,

                /// <summary>
                /// The use genes.
                /// </summary>
                UseGenes,
            }

            /// <summary>
            /// Gets the name of the class implementing the executed task.
            /// </summary>
            /// <value>The name.</value>
            public override string Name
            {
                get
                {
                    return typeof(TPredictor).Name;
                }
            }

            /// <summary>
            /// Gets the basic options data.
            /// </summary>
            protected Dictionary<BasicArgs, string> BasicOptionsData
            {
                get
                {
                    return new Dictionary<BasicArgs, string>
                    {
                        { BasicArgs.MapFileName, "Name of map to use for predictions" },
                        { BasicArgs.OutputFile, "Output file of predictions" },
                        { BasicArgs.UseGenes, "Flag to indicate if gene predictions should be made if underlying map is TSS-based" },
                        { BasicArgs.Threshold, "Optional threshold to limit predictions" },
                        { BasicArgs.ThresholdType, "Optiona parameter to indicate the type of threshold (" +
                            string.Join(", ", Enum.GetNames(typeof(ThresholdTypes))) + "); " +
                            ThresholdTypes.Score + " is default"
                        },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis for which this class is the factor with the given command line arguments.
            /// </summary>
            /// <param name="commandArgs">Command arguments.</param>
            public override void Execute(Args commandArgs)
            {
                var predictor = (TPredictor)System.Activator.CreateInstance(typeof(TPredictor));

                predictor.MapFileName = this.CommandArgs.StringArgs[BasicArgs.MapFileName.ToString()];
                predictor.OutputFile = this.CommandArgs.StringArgs[BasicArgs.OutputFile.ToString()];
                predictor.UseGenes = this.CommandArgs.Flags.Contains(BasicArgs.UseGenes.ToString());

                if (this.CommandArgs.StringArgs.ContainsKey(OptionalArgs.Threshold.ToString()))
                {
                    predictor.Threshold = double.Parse(this.CommandArgs.StringArgs[OptionalArgs.Threshold.ToString()]);
                }
                else
                {
                    predictor.Threshold = -1;
                }

                if (this.CommandArgs.StringArgs.ContainsKey(OptionalArgs.ThresholdType.ToString()))
                {
                    predictor.ThresholdType = (ThresholdTypes)Enum.Parse(
                        typeof(ThresholdTypes), 
                        this.CommandArgs.StringArgs[OptionalArgs.ThresholdType.ToString()]);
                }
                else
                {
                    predictor.ThresholdType = ThresholdTypes.Score;
                }

                this.ReflectArgs(predictor);

                predictor.Predict();
            }

            /// <summary>
            /// Reflects the arguments.
            /// </summary>
            public virtual void ReflectArgs(TPredictor predictor)
            {
            }
        }
    }
}