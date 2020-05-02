
namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Analyses;
    using Shared;

    public class MapStats : BaseAnalysis
    {
        /// <summary>
        /// Stat.
        /// </summary>
        public enum Stat
        {
            /// <summary>
            /// The link count.
            /// </summary>
            LinkCount,

            /// <summary>
            /// The tss count.
            /// </summary>
            TssCount,
        }

        /// <summary>
        /// The stat list.
        /// </summary>
        private List<Stat> statList;

        /// <summary>
        /// The stat calculate.
        /// </summary>
        private Dictionary<Stat, Func<string>> statCalc;

        /// <summary>
        /// The stat set.
        /// </summary>
        private HashSet<string> statSet = new HashSet<string>(Enum.GetNames(typeof(Stat)));

        /// <summary>
        /// The stat map.
        /// </summary>
        private Dictionary<string, Stat> statMap = Enum.GetValues(typeof(Stat))
            .Cast<Stat>()
            .ToDictionary(x => x.ToString(), x => x);

        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.MapStats"/> class.
        /// </summary>
        public MapStats()
            : base("", "", "")
        {
            this.MapFile.Register(this);
        }

        /// <summary>
        /// The map file.
        /// </summary>
        public MapFileElement MapFile = new MapFileElement(string.Empty);

        /// <summary>
        /// Gets or sets the stats.
        /// </summary>
        /// <value>The stats.</value>
        public string Stats { get; set; }

        /// <summary>
        /// Gets the stat list.
        /// </summary>
        /// <value>The stat list.</value>
        public List<Stat> StatList
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.statList,
                    () => this.Stats.Split(',').Select(x => (Stat)Enum.Parse(typeof(Stat), x)).ToList());
            }
        }

        /// <summary>
        /// Gets the stat calculators.
        /// </summary>
        /// <value>The stat calculators.</value>
        public Dictionary<Stat, Func<string>> StatCalculators
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.statCalc,
                    () => this.GetType().GetMethods()
                        .Where(x => this.statSet.Contains(x.Name))
                        .ToDictionary(x => this.statMap[x.Name], x => 
                        {
                            Func<string> f = () => (string)x.Invoke(this, null); 
                            return f;
                        }));
            }
        }

        /// <summary>
        /// Execute this instance.
        /// </summary>
        public void Execute()
        {
            Console.WriteLine(string.Join("\n", this.StatList.Select(x => string.Format("{0}\t{1}", x, this.StatCalculators[x]()))));
        }

        /// <summary>
        /// Tsses the count.
        /// </summary>
        /// <returns>The count.</returns>
        public string TssCount()
        {
            return this.MapFile.Element.Keys.Count.ToString();
        }

        /// <summary>
        /// Executor.
        /// </summary>
        public class Executor : IAnalysisExecutor<MapStats, MapStats.Executor.Arguments>
        {
            /// <summary>
            /// Arguments.
            /// </summary>
            public enum Arguments
            {
                /// <summary>
                /// The stat.
                /// </summary>
                Stats,

                /// <summary>
                /// The name of the map file.
                /// </summary>
                MapFileName,

                /// <summary>
                /// The type of the map file link.
                /// </summary>
                MapFileLinkType,

                /// <summary>
                /// The map file max range.
                /// </summary>
                MapFileMaxRange,

                /// <summary>
                /// The map file promoter range.
                /// </summary>
                MapFilePromoterRange,

                /// <summary>
                /// The map file confidence threshold.
                /// </summary>
                MapFileConfidenceThreshold,
            }

            /// <summary>
            /// Gets the description of the task.
            /// </summary>
            /// <value>The description.</value>
            public override string Description
            {
                get
                {
                    return "Analysis to collate stats on a provided map";
                }
            }

            /// <summary>
            /// Gets the command line option descriptions.
            /// </summary>
            /// <value>The options data.</value>
            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { 
                            Arguments.Stats, 
                            string.Format("CSV of stats to calculate on the given map ({0})", string.Join(", ", Enum.GetNames(typeof(MapStats.Stat)))) 
                        },
                        { Arguments.MapFileName, "Map file name" },
                    };
                }
            }

            /// <summary>
            /// Execute the analysis for which this class is the factor with the given command line arguments.
            /// </summary>
            /// <param name="commandArgs">Command arguments.</param>
            public override void Execute(Args commandArgs)
            {
                var analysis = new MapStats();

                this.ReflectStringArgs(analysis, new Arguments[]
                {
                    Arguments.Stats,
                    Arguments.MapFileName
                });

                this.ReflectOptionalStringArgs(analysis, new Arguments[]
                {
                    Arguments.MapFileLinkType,
                    Arguments.MapFileMaxRange,
                    Arguments.MapFilePromoterRange,
                    Arguments.MapFileConfidenceThreshold,
                });

                analysis.Execute();
            }
        }
    }
}

