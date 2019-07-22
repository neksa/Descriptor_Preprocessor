//--------------------------------------------------------------------------------
// <copyright file="Regression.cs" 
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
    using System.Linq;
    using System.IO;
    using System.Text.RegularExpressions;
    using Shared;
    using Genomics;

    /// <summary>
    /// Regresses C# data in R and manages data transfer into and out of R.
    /// </summary>
    public static class Regression
    {
        private static readonly string ScriptPath = "source/CisMapper/Data/RegressionHelpers.R";

        /// <summary>
        /// TSS column header names
        /// </summary>
        /// <returns>The tss column names.</returns>
        /// <param name="chipDataSets">Chip data sets.</param>
        public static List<string> MakeTssColumnNames(List<string> chipDataSets)
        {
            return chipDataSets.Select(x => string.Format("Tss_{0}", x)).ToList();
        }

        /// <summary>
        /// DRM column header names
        /// </summary>
        /// <returns>The tss column names.</returns>
        /// <param name="chipDataSets">Chip data sets.</param>
        public static List<string> MakeLocusColumnNames(List<string> chipDataSets)
        {
            return chipDataSets.Select(x => string.Format("Drm_{0}", x)).ToList();
        }

        /// <summary>
        /// Regress a set of columns versus expression in a given table data file
        /// </summary>
        /// <returns>The columns.</returns>
        /// <param name="tableFile">Data file name of regression input data</param>
        /// <param name="columns">Column names to regress</param>
        public static double RegressColumns(
            string tableFile,
            IEnumerable<string> columns)
        {
            string scriptName = tableFile.Replace("tsv", "R");

            string outputFile = RegressColumnsCore(
                "regress.columns.r2",
                scriptName.Replace(".R", ".yhat"),
                scriptName,
                tableFile,
                tableFile.Replace(".tsv", ".pdf"),
                columns);

            return ReadR2s(scriptName, new string[] { outputFile } )[0];
        }

        /// <summary>
        /// Regresses the columns logit.
        /// </summary>
        /// <returns>The columns logit.</returns>
        /// <param name="tableFile">Table file.</param>
        /// <param name="columns">Columns.</param>
        public static double RegressColumnsLogit(
            string tableFile,
            IEnumerable<string> columns)
        {
            string scriptName = tableFile.Replace("tsv", "R");

            string outputFile = RegressColumnsCore(
                "regress.columns.logit.r2",
                scriptName.Replace(".R", "." + 976 + ".R2s"),
                scriptName,
                tableFile,
                tableFile.Replace("tsv", "pdf"),
                columns);

            return ReadR2s(scriptName, new string[] { outputFile } )[0];
        }

        /// <summary>
        /// Regresses the columns logit lasso.
        /// </summary>
        /// <returns>The columns logit lasso.</returns>
        /// <param name="scriptName">Script name.</param>
        /// <param name="tableFile">Table file.</param>
        /// <param name="columns">Columns.</param>
        public static double RegressColumnsLogitLasso(
            string tableFile,
            IEnumerable<string> columns)
        {
            string scriptName = tableFile.Replace("tsv", "R");

            string outputFile = RegressColumnsCore(
                "regress.columns.logit.lasso",
                scriptName.Replace(".R", "." + 976 + ".R2s"),
                scriptName,
                tableFile,
                tableFile.Replace("tsv", "pdf"),
                columns);

            return ReadR2s(scriptName, new string[] { outputFile } )[0];
        }

        /// <summary>
        /// Regresses the columns lasso.
        /// </summary>
        /// <returns>The columns lasso.</returns>
        /// <param name="tableFile">Table file.</param>
        /// <param name = "columnTypes"></param>
        /// <param name="columns">Columns.</param>
        public static Tuple<double, double> RegressColumnsLasso(
            string tableFile,
            string columnTypes,
            IEnumerable<string> columns)
        {
            string scriptName = tableFile.Replace("tsv", "R");

            string outputFile = RegressColumnsCore(
                "regress.columns.lasso",
                scriptName.Replace(".R", "." + 976 + ".R2s"),
                scriptName,
                tableFile,
                tableFile.Replace(".tsv", "_" + columnTypes + ".pdf"),
                columns);

            return  new Tuple<double, double>(
                ReadR2s(scriptName, new string[] { outputFile })[0],
                ReadErrors(scriptName, new string[] { outputFile })[0]);
        }

        /// <summary>
        /// Regresses the columns using cross-validated lasso with an internal CV to tune lambda.
        /// </summary>
        /// <returns>The columns lasso.</returns>
        /// <param name="tableFile">Table file.</param>
        /// <param name = "columnTypes"></param>
        /// <param name="columns">Columns.</param>
        public static Tuple<double, double> RegressColumnsLassoCV(
            string tableFile,
            string columnTypes,
            IEnumerable<string> columns)
        {
            string scriptName = tableFile.Replace("tsv", "R");
            string predictionFile = scriptName.Replace(".R", ".yhat");

            string outputFile = RegressColumnsCore(
                "regress.columns.lasso.cv",
                predictionFile,
                scriptName,
                tableFile,
                tableFile.Replace(".tsv", "_" + columnTypes + ".pdf"),
                columns);

            return new Tuple<double, double>(
                ReadR2s(scriptName, new string[] { outputFile })[0],
                ReadErrors(scriptName, new string[] { outputFile })[0]);
        }



        /// <summary>
        /// Predicts the columns.
        /// </summary>
        /// <returns>The columns.</returns>
        /// <param name="tableFile">Table file.</param>
        /// <param name="columns">Columns.</param>
        public static Dictionary<string, double> PredictColumns(
            string tableFile,
            IEnumerable<string> columns)
        {
            string scriptName = tableFile.Replace("tsv", "R");
            string predictionFile = scriptName.Replace(".R", ".yhat");

            RegressColumnsCore(
                "regress.columns.predict",
                predictionFile,
                scriptName,
                tableFile,
                tableFile.Replace("tsv", "pdf"),
                columns);

            return GetPredictions(predictionFile);
        }


        /// <summary>
        /// Regresses the columns core.
        /// </summary>
        /// <returns>The columns core.</returns>
        /// <param name="function">Function.</param>
        /// <param name="predictionFile">Prediction file.</param>
        /// <param name="scriptName">Script name.</param>
        /// <param name="tableFile">Table file.</param>
        /// <param name="columns">Columns.</param>
        public static string RegressColumnsCore(
            string function,
            string predictionFile,
            string scriptName, 
            string tableFile,
            string plotFile,
            IEnumerable<string> columns)
        {
            using (TextWriter tw = Helpers.CreateStreamWriter(scriptName))
            {
                WriteSource(tw);

                tw.WriteLine("columns <- c(\"{0}\");", string.Join("\", \"", columns));
                tw.WriteLine("file <- \"" + tableFile + "\"");

                tw.WriteLine("r2filename <- \"" + predictionFile + "\"");
                tw.WriteLine("plotfile <- \"" + plotFile + "\"");

                tw.WriteLine(function + "(file, columns, r2filename, plotfile)");

                tw.WriteLine("q(save=\"no\")");
            }

            string outputFile = R.MakeScriptOutputFile(scriptName);
            R.RunScript(scriptName, outputFile);

            return outputFile;
        }

        /// <summary>
        /// Gets the predictions.
        /// </summary>
        /// <returns>The predictions.</returns>
        /// <param name="predictionFile">Prediction file.</param>
        public static Dictionary<string, double> GetPredictions(string predictionFile)
        {
            return Helpers.GetFileDataLines(predictionFile, true)
                .Select(line =>
                {
                    var fields = line.Split(',').Select(f => f.Trim()).ToArray();

                    return new 
                    {
                        Name = fields[0].Replace("\"", ""),
                        Value = double.Parse(fields[1])
                    };
                })
                .ToDictionary(x => x.Name, x => x.Value);
        }

        /// <summary>
        /// Read R2 data from all or a selected batch
        /// </summary>
        /// <returns>The r2s.</returns>
        /// <param name="scriptName">Script name.</param>
        /// <param name="batchCount">Batch count.</param>
        public static double[] ReadR2s(
            string scriptName, 
            string[] outputFiles)
        {
            // | Set Index | R-squared | 
            Regex setRs = new Regex("Rsquareds\\|.*\\|.*\\|");

            return ReadRegressionField(scriptName, outputFiles, setRs);
        }



        /// <summary>
        /// Reads the errors.
        /// </summary>
        /// <returns>The errors.</returns>
        /// <param name="scriptName">Script name.</param>
        /// <param name="outputFiles">Output files.</param>
        public static double[] ReadErrors(
            string scriptName, 
            string[] outputFiles)
        {
            // | Set Index | R-squared | 
            Regex error = new Regex("Error\\|.*\\|.*\\|");

            return ReadRegressionField(scriptName, outputFiles, error);
        }

        /// <summary>
        /// Reads the regression field.
        /// </summary>
        /// <returns>The regression field.</returns>
        /// <param name="scriptName">Script name.</param>
        /// <param name="outputFiles">Output files.</param>
        /// <param name="pattern">Pattern.</param>
        private static double[] ReadRegressionField(
            string scriptName, 
            string[] outputFiles, 
            Regex pattern)
        {
            Func<TextReader, IEnumerable<double>> parseFile = tr => 
            {
                return tr.ReadToEnd().Split('\n')
                    .Where(line => pattern.IsMatch(line))
                    .Select(line => double.Parse(line.Trim().Split('|')[2]));
            };

            return outputFiles
                .Select(x => 
                {
                    using (TextReader tr = new StreamReader(x))
                    {
                        return parseFile(tr);
                    }
                })
                .SelectMany(x => x).ToArray();
        }

        /// <summary>
        /// Writes the source invocation for the R source
        /// </summary>
        /// <param name="tw">Tw.</param>
        public  static void WriteSource(TextWriter tw)
        {
            tw.WriteLine(SourceHelpers());
        }

        /// <summary>
        /// Sources the helpers.
        /// </summary>
        /// <returns>The helpers.</returns>
        public static string SourceHelpers()
        {
            return "source(\"" + ScriptPath + "\");";
        }
    }
}

