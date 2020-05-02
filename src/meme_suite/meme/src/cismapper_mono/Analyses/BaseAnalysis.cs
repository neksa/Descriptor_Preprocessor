//--------------------------------------------------------------------------------
// <copyright file="BaseAnalysis.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using System.Text.RegularExpressions;
    using Data;
    using Genomics;
    using Shared;
    using Tools;

    /// <summary>
    /// Base experiment data.
    /// </summary>
    public class BaseAnalysis : IUnknown
    {
        private Dictionary<string, Action<string>> elementArgRegistry = new Dictionary<string, Action<string>>();

        /// <summary>
        /// Gets or sets the output files.
        /// </summary>
        /// <value>The output files.</value>
        protected OutputFiles OutputFiles { get; set; }

        /// <summary>
        /// Initializes a new instance of the <see cref="Experiments.BaseExperiment`2"/> class.
        /// </summary>
        /// <param name="outputRoot">Output root.</param>
        /// <param name="plotRoot">Plot root.</param>
        /// <param name="scriptRoot">Script root.</param>
        public BaseAnalysis(string outputRoot, string plotRoot, string scriptRoot)
        {
            this.OutputFiles = new OutputFiles(outputRoot, plotRoot, scriptRoot);
        }

        public Dictionary<string, Action<string>> ElementArgRegistry
        {
            get
            {
                return this.elementArgRegistry;
            }
        }

        /// <summary>
        /// Table file name.
        /// </summary>
        /// <returns>The file name.</returns>
        /// <param name="prefix">File stem prefix prefix.</param>
        /// <param name="core">File stem core.</param>
        /// <param name="suffix">File stem suffix.</param>
        public string TableFile(
            string prefix, 
            string core = default(string),
            string suffix = default(string))
        {
            return this.OutputFiles.TableFile(this.FileStem(prefix, core, suffix));
        }

        /// <summary>
        /// Script file name.
        /// </summary>
        /// <returns>The file name.</returns>
        /// <param name="prefix">File stem prefix prefix.</param>
        /// <param name="core">File stem core.</param>
        /// <param name="suffix">File stem suffix.</param>
        public string ScriptFile(
            string prefix, 
            string core = default(string),
            string suffix = default(string))
        {
            return this.OutputFiles.ScriptFile(this.FileStem(prefix, core, suffix));
        }

        /// <summary>
        /// Plot file name.
        /// </summary>
        /// <returns>The file name.</returns>
        /// <param name="prefix">File stem prefix prefix.</param>
        /// <param name="core">File stem core.</param>
        /// <param name="suffix">File stem suffix.</param>
        public string PlotFile(
            string prefix, 
            string core = default(string),
            string suffix = default(string))
        {
            return this.OutputFiles.PlotFile(this.FileStem(prefix, core, suffix));
        }

        /// <summary>
        /// Supplementary file name.
        /// </summary>
        /// <returns>The plot file.</returns>
        /// <param name="prefix">Prefix.</param>
        /// <param name="core">Core.</param>
        /// <param name="suffix">Suffix.</param>
        public string SupplementaryPlotFile(
            string prefix, 
            string core = default(string),
            string suffix = default(string))
        {
            return this.OutputFiles.SupplementaryPlotFile(this.FileStem(prefix, core, suffix));
        }

        /// <summary>
        /// Creates a file stem core
        /// </summary>
        /// <returns>The stem core.</returns>
        /// <param name="fields">Fields.</param>
        protected string FileStemCore(string[] fields)
        {
            return string.Join("_", fields);
        }

        /// <summary>
        /// Creates a file stem
        /// </summary>
        /// <returns>The file stem.</returns>
        /// <param name="prefix">File stem prefix prefix.</param>
        /// <param name="core">File stem core.</param>
        /// <param name="suffix">File stem suffix.</param>
        private string FileStem(
            string prefix, 
            string core = default(string),
            string suffix = default(string))
        {
            return string.Join(
                "_",
                new List<string>
            {
                prefix,
                core,
                suffix
            }
                .Where(x => !string.IsNullOrEmpty(x)));
        }
    }
}

