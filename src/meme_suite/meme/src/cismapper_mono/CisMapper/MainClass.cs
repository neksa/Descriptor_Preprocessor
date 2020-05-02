//--------------------------------------------------------------------------------
// <copyright file="MainClass.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace CisMapper
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using System.Text.RegularExpressions;
    using Analyses;
    using Genomics;
    using Shared;
    using Tools;

    /// <summary>
    /// Main class. Bundle of all analyses and experiments.
    /// </summary>
    public class MainClass   
    {
        /// <summary>
        /// The experiments.
        /// </summary>
        private static readonly Dictionary<string, Type> Tasks = new List<KeyValuePair<string, Type>>
        {
            // Analyses
            ITaskExecutor.RegisterTask<MapStats.Executor>(),
            ITaskExecutor.RegisterTask<CorrelationMapBuilder.Executor>(),
            ITaskExecutor.RegisterTask<CorrelationMapEligibleGenes.Executor>(),
            ITaskExecutor.RegisterTask<NullMapBuilder.Executor>(),
            ITaskExecutor.RegisterTask<BrowserTrackFromMap.Executor>(),
            ITaskExecutor.RegisterTask<ChIPSeqPeakBias.Executor>(),
            ITaskExecutor.RegisterTask<MapFilter.Executor>(),
            ITaskExecutor.RegisterTask<ConvertMapToGenes.Executor>(),
            ITaskExecutor.RegisterTask<ConvertMapToScoredRegionPairs.Executor>(),
            ITaskExecutor.RegisterTask<ConvertContactsToMap.Executor>(),
            ITaskExecutor.RegisterTask<MergeCorrelationAndNullMaps.Executor>(),
            ITaskExecutor.RegisterTask<PredictLinks.Executor>(),
            ITaskExecutor.RegisterTask<PredictTargets.Executor>(),
            ITaskExecutor.RegisterTask<PredictElements.Executor>(),
        }
            .ToDictionary(x => x.Key, x => x.Value);

        /// <summary>
        /// The entry point of the program, where the program control starts and ends.
        /// </summary>
        /// <param name="args">The command-line arguments.</param>
        public static void Main(string[] args)
        {
            System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("en-US");
            System.Threading.Thread.CurrentThread.CurrentCulture = ci;
            System.Threading.Thread.CurrentThread.CurrentUICulture = ci;

            Args a = new Args(args);

            if ((a.Flags.Contains("Help") && a.StringArgs.Count == 0 && a.IntArgs.Count == 0) ||
                (a.Flags.Count == 0 && a.StringArgs.Count == 0 && a.IntArgs.Count == 0))
            {
                ReportHelp("Tasks", Tasks);
                return;
            }

            string taskName = a.StringArgs["Mode"];

            if (Tasks.ContainsKey(taskName))
            {
                var executor = (ITaskExecutor)Activator.CreateInstance(Tasks[taskName]);
                executor.CommandArgs = a;

                if (a.Flags.Contains("Help"))
                {
                    Console.WriteLine(executor.HelpString);
                }
                else
                {
                    executor.Execute(a);
                }

                return;
            }

            switch (a.StringArgs["Mode"])
            {
                default:
                    Console.WriteLine("Invalid argument: {0}", a.StringArgs["Mode"]);
                    break;
            }
        }

        /// <summary>
        /// Reports the help data on all registered items
        /// </summary>
        /// <param name="type">Experiment or analysis type</param>
        /// <param name="items">Experiments or analysis to report help infor on</param>
        private static void ReportHelp(string type, Dictionary<string, Type> items)
        {
            Console.WriteLine(
                "#######################################\n#" + type + "\n#######################################\n{0}\n",
                string.Join(
                    "\n",
                    items.Select(x => 
                    {
                        var e = (ITaskExecutor)Activator.CreateInstance(x.Value);
                        return new 
                        {
                            e.Name,
                            e.Description
                        };
                    }).OrderBy(x => x.Name).Select(e => string.Format("===============================================\n{0}\n-----------------------------------------------\n\t{1}\n", e.Name, e.Description))));
        }
    }
}
