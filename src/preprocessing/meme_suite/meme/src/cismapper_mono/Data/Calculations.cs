//--------------------------------------------------------------------------------
// <copyright file="Calculations.cs" 
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
    /// Calculations that roundtrip from C# -> R -> C#.
    /// </summary>
    public static class Calculations
    {
        /// <summary>
        /// Convert z scores to log pvalues using R
        /// </summary>
        /// <returns>The to log pvalue.</returns>
        /// <param name="z">The z coordinate.</param>
        public static double[] ZscoreToLogPvalue(IEnumerable<double> z)
        {
            string dataFile = Tables.ToNsvFile(z);
            //string scriptName = "./zscore.R";
            string scriptName = Path.GetTempFileName() + ".R";
            string outputFile = scriptName + "out.tsv";

            R.WriteAndRunScript(new string[]
            {
                R.ReadColumn(R.Variables.X, dataFile),
                string.Format("{0} <- pnorm(-abs({1}), log=T) + log(2)", R.Variables.D, R.Variables.X),
                R.WriteCsv(R.Variables.D, outputFile),
            },
                scriptName);

            return R.FromCsvFile(outputFile);
        }

        /// <summary>
        /// Convert z scores to bonferroni corrected pvalues using R
        /// </summary>
        /// <returns>The to log pvalue.</returns>
        /// <param name="z">The z coordinate.</param>
        public static double[] ZscoreToCorrectedPvalue(double[] z, int n)
        {
            //string dataFile = "./zdatac.tsv";
            string dataFile = Path.GetTempFileName() + ".tsv";
            using (TextWriter tw = Helpers.CreateStreamWriter(dataFile))
            {
                tw.WriteLine(string.Join("\n", z.Select(x => -Math.Abs(x))));
            }

            //string scriptName = "./zscorec.R";
            string scriptName = Path.GetTempFileName() + ".R";
            using (TextWriter tw = Helpers.CreateStreamWriter(scriptName))
            {
                tw.WriteLine("x <- read.table('" + dataFile.Replace(@"\", @"\\") + "')");
                tw.WriteLine("d <- p.adjust(exp(pnorm(x[,1], log=T) + log(2)), method=\"bonferroni\")");
                tw.WriteLine("for (i in 1:dim(x)[1]) { print(d[i]) }");
                //tw.WriteLine(string.Join("\n", z.Select(x => string.Format("bcp = p.adjust(exp(pnorm(" + -Math.Abs(x) + ", log=T) + log(2)), method=\"bonferroni\", n=" + n + "); print(sprintf(\"bcp=%e=\", bcp));"))));
            }

            string outputFile = Path.GetTempFileName();
            //string outputFile = scriptName + "out";
            R.RunScript(scriptName, outputFile);

            Regex lp = new Regex("\\[1\\] ");

            using (TextReader tr = new StreamReader(outputFile))
            {
                return tr.ReadToEnd().Split('\n').Where(line => lp.IsMatch(line))
                    .Select(line =>
                    {
                        return double.Parse(line.Split(' ')[1].Trim());
                    }).ToArray();
            }
        }

        /// <summary>
        /// Performs a sign test on the first column of a matrix relative to the other columns
        /// </summary>
        /// <returns>The ranked matrix.</returns>
        /// <param name="matrixFileName">Matrix file name.</param>
        public static double[] TestRankedMatrix(string matrixFileName)
        {
            string scriptName = matrixFileName + ".R";
            using (TextWriter tw = Helpers.CreateStreamWriter(scriptName))
            {
                Regression.WriteSource(tw);

                tw.WriteLine("test.ranked.matrix('" + matrixFileName + "')");
            }

            string outputFile = scriptName + "out";
            R.RunScript(scriptName, outputFile);

            Regex lp = new Regex("\\[1\\] ");

            using (TextReader tr = new StreamReader(outputFile))
            {
                return tr.ReadToEnd().Split('\n').Where(line => lp.IsMatch(line))
                    .Select(line =>
                    {
                        //Console.WriteLine(line);
                        return double.Parse(line.Split(' ')[1]);
                    }).ToArray();
            }
        }

        public static double TestSetOverlap(int overlap, int sample1, int sample2, int total)
        {
            int[] contingency = { overlap, sample1 - overlap, sample2 - overlap, total - sample1 - sample2 + overlap };

            string scriptName = Path.GetTempFileName() + ".R";
            using (TextWriter tw = Helpers.CreateStreamWriter(scriptName))
            {
                tw.WriteLine("print(fisher.test(matrix(c(" + string.Join(",", contingency) + "), 2, 2))$p.value)");
            }

            string outputFile = scriptName + "out";
            R.RunScript(scriptName, outputFile);

            Regex lp = new Regex("\\[1\\] ");

            using (TextReader tr = new StreamReader(outputFile))
            {
                return tr.ReadToEnd().Split('\n').Where(line => lp.IsMatch(line))
                    .Select(line =>
                    {
                        //Console.WriteLine(line);
                        return double.Parse(line.Split(' ')[1]);
                    }).First();
            }
        }

        public enum TestSide
        {
            TwoSided,
            Less,
            Greater,
        }

        public class TestResult
        {
            public double TestStatistic { get; set; }

            public double Pvalue { get; set; }
        }

        public static Dictionary<TestSide, string> Alternatives = new Dictionary<TestSide, string>
        {
            { TestSide.TwoSided, "'two.sided'" },
            { TestSide.Less, "'less'" },
            { TestSide.Greater, "'greater'" },
        };

        /// <summary>
        /// Kolmogorov-Smirnov p-value
        /// </summary>
        /// <returns>The test.</returns>
        /// <param name="tableFile">Table file.</param>
        /// <param name="dataColumn">Data column.</param>
        /// <param name="factorColumn">Factor column.</param>
        /// <param name="factor1">Factor1.</param>
        /// <param name="factor2">Factor2.</param>
        public static TestResult KSTest(
            TestSide side,
            string tableFile, 
            string dataColumn, 
            string factorColumn, 
            string factor1, 
            string factor2)
        {
            string tempFile = Path.GetTempFileName();
            string pvalueFile = tempFile + ".tsv";
            R.WriteAndRunScript(new string[]
            {
                R.ReadTable(R.Variables.X, tableFile, true),
                "y <- x[which(x$Confidence == 'Linked'),]",
                R.FactorColumn("y", factorColumn),
                string.Format("x1 <- y${0}[which(y${1} == '{2}')]", dataColumn, factorColumn, factor1),
                string.Format("x2 <- y${0}[which(y${1} == '{2}')]", dataColumn, factorColumn, factor2),
                string.Format("result <- ks.test(x1,x2,alternative={0})", Alternatives[side]),
                "m <- matrix(c(result$statistic, result$p.value), 1, 2)",
                R.WriteTable("m", pvalueFile)
            },
                tempFile + ".R");

            var dataLine = R.ReadDataLine(pvalueFile);

            return new TestResult
            {
                TestStatistic = dataLine[0],
                Pvalue = dataLine[1]
            };
        }
    }
}

