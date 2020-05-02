//--------------------------------------------------------------------------------
// <copyright file="R.cs" 
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
    using System.Diagnostics;
    using System.IO;
    using System.Linq;
    using System.Text.RegularExpressions;
    using Genomics;
    using Shared;

    /// <summary>
    /// Wrapper of basic R interfaces
    /// </summary>
    public static class R
    {
        /// <summary>
        /// The executable.
        /// </summary>
        private static string executable;

        /// <summary>
        /// The executable lock.
        /// </summary>
        private static object executableLock = new object();

        /// <summary>
        /// The library names.
        /// </summary>
        private static Dictionary<Libraries, string> libraryNames = new Dictionary<Libraries, string>
        {
            { Libraries.GGPlot, "ggplot2" },
            { Libraries.ColorBrewer, "RColorBrewer" },
            { Libraries.GridExtra, "gridExtra" },
            { Libraries.LARS, "lars" },
            { Libraries.Plsdof, "plsdof" },
            { Libraries.Lattice, "lattice" },
            { Libraries.Boot, "boot" },
            { Libraries.KlaR, "klaR" },
            { Libraries.PRoc, "pROC" },
            { Libraries.Scales, "scales" },
        };

        /// <summary>
        /// R Libraries
        /// </summary>
        public enum Libraries
        {
            /// <summary>
            /// The GG plot library
            /// </summary>
            GGPlot,

            /// <summary>
            /// The color brewer.
            /// </summary>
            ColorBrewer,

            /// <summary>
            /// The grid extra.
            /// </summary>
            GridExtra,

            /// <summary>
            /// The scales.
            /// </summary>
            Scales,

            /// <summary>
            /// The lattice.
            /// </summary>
            Lattice,

            /// <summary>
            /// The LARS least angle regression.
            /// </summary>
            LARS,

            /// <summary>
            /// The bootstrap library.
            /// </summary>
            Boot,

            /// <summary>
            /// The plsdof.
            /// </summary>
            Plsdof,

            /// <summary>
            /// The kla r.
            /// </summary>
            KlaR,

            /// <summary>
            /// The P roc.
            /// </summary>
            PRoc,
        }

        /// <summary>
        /// Gets the executable.
        /// </summary>
        /// <value>The executable.</value>
        public static string Executable
        {
            get
            {
                lock (executableLock)
                {
                    return Helpers.CheckInit(
                        ref executable,
                        () =>
                        {
                            PlatformID p = Environment.OSVersion.Platform;
                            if (p == PlatformID.MacOSX || p == PlatformID.Unix)
                            {
                                //return "/opt/local/bin/R";
                                return "R";
                            }
                            else
                            {
                                return "R.exe";
                            }
                        });
                }
            }
        }

        /// <summary>
        /// Loads the library.
        /// </summary>
        /// <returns>The library.</returns>
        /// <param name="library">Library to load.</param>
        public static string LoadLibrary(Libraries library)
        {
            return string.Format("library({0})\n", libraryNames[library]);
        }

        /// <summary>
        /// Create an output file name for a batch run of a script
        /// </summary>
        /// <returns>The script batch output file.</returns>
        /// <param name="scriptName">Script name.</param>
        public static string MakeScriptOutputFile(string scriptName)
        {
            return scriptName.Replace(".R", ".Rout");
        }

        /// <summary>
        /// Reads a the first column from a table file into the named variable
        /// </summary>
        /// <returns>The column.</returns>
        /// <param name="variable">Variable name to load the column into.</param>
        /// <param name="filename">Filename of the data.</param>
        public static string ReadColumn(string variable, string filename)
        {
            return variable + " <- read.table('" + filename + "')[,1]\n";
        }

        /// <summary>
        /// Reads a table file into the name variable
        /// </summary>
        /// <returns>The column.</returns>
        /// <param name="variable">Variable name to load data into.</param>
        /// <param name="filename">Filename of the data.</param>
        /// <param name="hasHeader">If the table has a header line.</param>
        public static string ReadTable(string variable, string filename, bool hasHeader)
        {
            return variable + " <- read.table('" + filename + "', header=" + (hasHeader ? "T" : "F") + ")\n";
        }

        /// <summary>
        /// Gets the first data column from a csv file
        /// </summary>
        /// <returns>The csv file.</returns>
        /// <param name="fileName">File name.</param>
        /// <param name="parser">Parser of the data.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public static T[] FromCsvFile<T>(string fileName, Func<string, T> parser)
        {
            return Helpers.GetFileDataLines(fileName, true).Select(line => parser(line.Split(',')[1])).ToArray();
        }

        /// <summary>
        /// Gets the first data column from a csv file
        /// </summary>
        /// <returns>The csv file.</returns>
        /// <param name="fileName">File name.</param>
        public static double[] FromCsvFile(string fileName)
        {
            return FromCsvFile(
                fileName, 
                (field) => field == "Inf" ? double.PositiveInfinity : 
                    field == "-Inf" ? double.NegativeInfinity : 
                    double.Parse(field));
        }

        /// <summary>
        /// Writes the table.
        /// </summary>
        /// <returns>The table writing line.</returns>
        /// <param name="tableName">Variable name.</param>
        /// <param name="fileName">File name.</param>
        public static string WriteTable(string variable, string fileName)
        {
            return string.Format("write.table({0}, file='{1}')", variable, fileName);
        }

        /// <summary>
        /// Writes the single value.
        /// </summary>
        /// <returns>The single value.</returns>
        /// <param name="variable">Variable.</param>
        /// <param name="fileName">File name.</param>
        public static string WriteSingleValue(string variable, string fileName)
        {
            return WriteTable(variable, fileName);
        }

        /// <summary>
        /// Reads the single value.
        /// </summary>
        /// <returns>The single value.</returns>
        /// <param name="fileName">File name.</param>
        public static double ReadSingleValue(string fileName)
        {
            return double.Parse(Helpers.GetFileDataLines(fileName, true).First().Split(' ').Last());
        }

        public static double[] ReadDataLine(string fileName)
        {
            return Helpers.GetFileDataLines(fileName, true)
                .First()
                .Split(' ')
                .Where((x, i) => i > 0)
                .Select(x => double.Parse(x))
                .ToArray();
        }

        /// <summary>
        /// Writes the csv.
        /// </summary>
        /// <returns>The csv writing line.</returns>
        /// <param name="variable">Variable name.</param>
        /// <param name="fileName">File name.</param>
        public static string WriteCsv(string variable, string fileName)
        {
            return string.Format("write.csv({0}, file='{1}')", variable, fileName);
        }

        /// <summary>
        /// Access the named column of the variable.
        /// </summary>
        /// <param name="variable">Variable.</param>
        /// <param name="field">Field.</param>
        public static string Column(string variable, string field)
        {
            return string.Format("{0}${1}", variable, field);
        }

        /// <summary>
        /// Convert a column to a factor
        /// </summary>
        /// <returns>The column.</returns>
        /// <param name="variable">Variable.</param>
        /// <param name="field">Field.</param>
        public static string FactorColumn(string variable, string field)
        {
            return Column(variable, field) + " <- factor(" + Column(variable, field) + ");";
        }

        /// <summary>
        /// Binds the columns.
        /// </summary>
        /// <returns>The columns.</returns>
        /// <param name="tableName">Table name.</param>
        /// <param name="vars">Variables.</param>
        /// <param name="varFields">Variable fields.</param>
        public static string BindColumns(string tableName, string[] vars, string[] varFields)
        {
            if (vars.Length != varFields.Length)
            {
                throw new Exception("Error: cannot bind columns. Vars and fields do not match");
            }

            return string.Format(
                "{0} <- data.frame(cbind({1}))\ncolnames({0}) <- c({2})",
                tableName,
                JoinEnumerableCsv(vars.Select((v, ci) => R.Column(v, varFields[ci]))),
                JoinEnumerableCsv(varFields.Select(x => "'" + x + "'")));
        }

        /// <summary>
        /// Writes the and run script returning the lines of the output file
        /// </summary>
        /// <returns>The and run script.</returns>
        /// <param name="scriptLines">Script lines.</param>
        /// <param name="scriptFile">Script file.</param>
        public static string[] WriteAndRunScript(string[] scriptLines, string scriptFile)
        {
            using (TextWriter tr = Helpers.CreateStreamWriter(scriptFile))
            {
                tr.WriteLine(string.Join("\n", scriptLines));
            }

            string outputFile = scriptFile + "out";
            RunScript(scriptFile, outputFile);

            using (TextReader tr = new StreamReader(outputFile))
            {
                return tr.ReadToEnd().Split('\n').Where(line => !string.IsNullOrWhiteSpace(line)).ToArray();
            }
        }

        /// <summary>
        /// Runs the script.
        /// </summary>
        /// <param name="scriptFile">Script file.</param>
        /// <param name="outputFile">Output file.</param>
        public static void RunScript(string scriptFile, string outputFile)
        {
            ExecuteProcess(Executable, string.Format("CMD BATCH --vanilla --slave {0} {1}", scriptFile, outputFile));
        }

        /// <summary>
        /// Joins the array into a csv.
        /// </summary>
        /// <returns>The array csv.</returns>
        /// <param name="data">Data.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public static string JoinEnumerableCsv<T>(IEnumerable<T> data)
        {
            return string.Join(", ", data.Select(x => x.ToString()));
        }

        /// <summary>
        /// Joins the lines.
        /// </summary>
        /// <returns>The lines.</returns>
        /// <param name="lines">Lines.</param>
        public static string JoinLines(IEnumerable<string> lines)
        {
            if (lines == null)
            {
                return string.Empty;
            }

            return string.Join("\n", lines);
        }

        /// <summary>
        /// Breaks the specified start, end and intervals.
        /// </summary>
        /// <param name="start">Start.</param>
        /// <param name="end">End.</param>
        /// <param name="intervals">Intervals.</param>
        public static double[] Breaks(double start, double end, int intervals)
        {
            double range = end - start;
            double interval = range / intervals;

            return Enumerable.Range(0, intervals + 1).Select(x => start + x * interval).ToArray();
        }

        /// <summary>
        /// Executes the process.
        /// </summary>
        /// <param name="command">Command to execute.</param>
        /// <param name="parameters">Parameters to the command.</param>
        private static void ExecuteProcess(string command, string parameters)
        {
            ProcessStartInfo info = new ProcessStartInfo();
            info.UseShellExecute = false;
            info.CreateNoWindow = false;

            info.WorkingDirectory = ".";
            info.FileName = command;
            info.Arguments = parameters;

            System.Diagnostics.Process proc = System.Diagnostics.Process.Start(info);
            proc.WaitForExit();
        }

        /// <summary>
        /// Variables to use in scripts.
        /// </summary>
        public static class Variables
        {
            /// <summary>
            /// The index of the plot.
            /// </summary>
            private static int index = 0;

            /// <summary>
            /// Gets the x variable.
            /// </summary>
            /// <value>The x.</value>
            public static string X 
            {
                get 
                { 
                    return "x";
                }
            }

            /// <summary>
            /// Gets the h variable.
            /// </summary>
            /// <value>The h.</value>
            public static string H
            {
                get
                {
                    return "h";
                }
            }

            /// <summary>
            /// Gets the d variable.
            /// </summary>
            /// <value>The d.</value>
            public static string D
            {
                get
                {
                    return "d";
                }
            }

            /// <summary>
            /// Gets the m variable.
            /// </summary>
            /// <value>The m.</value>
            public static string M
            {
                get
                {
                    return "m";
                }
            }

            /// <summary>
            /// Gets a unique plot variable. 
            /// </summary>
            /// <value>The plot variable.</value>
            public static string Plot
            {
                get
                {
                    return "plot" + index++;
                }
            }

            /// <summary>
            /// Gets a unique legend variable. 
            /// </summary>
            /// <value>The legend variable.</value>
            public static string Legend
            {
                get
                {
                    return "legend" + index++;
                }
            }

            /// <summary>
            /// Gets a unique table variable.
            /// </summary>
            /// <value>The table.</value>
            public static string Table
            {
                get
                {
                    return "table" + index++;
                }
            }
        }
    }
}