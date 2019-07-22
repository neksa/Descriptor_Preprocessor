//--------------------------------------------------------------------------------
// <copyright file="OutputFiles.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Data
{
    using System;
    using System.Collections.Generic;
    using Shared;

    /// <summary>
    /// Files.
    /// </summary>
    public class OutputFiles
    {
        /// <summary>
        /// The output root.
        /// </summary>
        private string outputRoot;

        /// <summary>
        /// The plot output root.
        /// </summary>
        private string plotOutputRoot;

        /// <summary>
        /// The supplementary plot output root.
        /// </summary>
        private string supplementaryPlotOutputRoot;

        /// <summary>
        /// The script output root.
        /// </summary>
        private string scriptOutputRoot;

        /// <summary>
        /// The output roots.
        /// </summary>
        private Dictionary<OutputFileTypes, string> outputFormats;

        /// <summary>
        /// Initializes a new instance of the <see cref="Data.Files"/> class.
        /// </summary>
        /// <param name="root">Root.</param>
        /// <param name="plotRoot">Plot root.</param>
        /// <param name="scriptRoot">Script root.</param>
        public OutputFiles(string root, string plotRoot, string scriptRoot)
        {
            this.outputRoot = root;
            this.plotOutputRoot = plotRoot;
            this.supplementaryPlotOutputRoot = plotRoot + "/SupPlot/";
            this.scriptOutputRoot = scriptRoot;
        }

        /// <summary>
        /// Output file types.
        /// </summary>
        public enum OutputFileTypes
        {
            /// <summary>
            /// Table file output
            /// </summary>
            Table,

            /// <summary>
            /// Script file.
            /// </summary>
            Script,

            /// <summary>
            /// The script output file.
            /// </summary>
            ScriptOutput,

            /// <summary>
            /// Plot file output
            /// </summary>
            Plot,

            /// <summary>
            /// Supplementary plot output.
            /// </summary>
            SupplementaryPlot
        }

        /// <summary>
        /// Gets the output formats.
        /// </summary>
        /// <value>The output formats.</value>
        protected Dictionary<OutputFileTypes, string> OutputFormats
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.outputFormats,
                    () => new Dictionary<OutputFileTypes, string>
                    {
                        { OutputFileTypes.Table, this.outputRoot + "{0}.tsv" },
                        { OutputFileTypes.Plot, this.plotOutputRoot + "{0}.pdf" },
                        { OutputFileTypes.SupplementaryPlot, this.supplementaryPlotOutputRoot + "{0}.pdf" },
                        { OutputFileTypes.Script, this.scriptOutputRoot + "{0}.R" },
                        { OutputFileTypes.ScriptOutput, this.scriptOutputRoot + "{0}.Rout" },
                });
            }
        }

        /// <summary>
        /// Creates a table file name.
        /// </summary>
        /// <returns>The file.</returns>
        /// <param name="fileStem">File stem.</param>
        public string TableFile(string fileStem)
        {
            return this.FormatOutputFileName(OutputFileTypes.Table, fileStem);
        }

        /// <summary>
        /// Creates a script file name.
        /// </summary>
        /// <returns>The file.</returns>
        /// <param name="fileStem">File stem.</param>
        public string ScriptFile(string fileStem)
        {
            return this.FormatOutputFileName(OutputFileTypes.Script, fileStem);
        }

        /// <summary>
        /// Creates a script output file name.
        /// </summary>
        /// <returns>The output file.</returns>
        /// <param name="fileStem">File stem.</param>
        public string ScriptOutputFile(string fileStem)
        {
            return this.FormatOutputFileName(OutputFileTypes.ScriptOutput, fileStem);
        }

        /// <summary>
        /// Creates a plot file name.
        /// </summary>
        /// <returns>The file.</returns>
        /// <param name="fileStem">File stem.</param>
        public string PlotFile(string fileStem)
        {
            return this.FormatOutputFileName(OutputFileTypes.Plot, fileStem);
        }

        /// <summary>
        /// Creates a plot file name.
        /// </summary>
        /// <returns>The file.</returns>
        /// <param name="fileStem">File stem.</param>
        public string SupplementaryPlotFile(string fileStem)
        {
            return this.FormatOutputFileName(OutputFileTypes.SupplementaryPlot, fileStem);
        }

        /// <summary>
        /// Formats the name of the output file.
        /// </summary>
        /// <returns>The output file name.</returns>
        /// <param name="type">Type.</param>
        /// <param name="fileStem">File stem.</param>
        private string FormatOutputFileName(
            OutputFileTypes type,
            string fileStem)
        {
            return string.Format(OutputFormats[type], fileStem);
        }
    }
}

