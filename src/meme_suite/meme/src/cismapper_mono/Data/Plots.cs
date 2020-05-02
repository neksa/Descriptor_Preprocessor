//--------------------------------------------------------------------------------
// <copyright file="Plots.cs" 
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
    /// Wraps R plotting to take in C# data and produce pdf output.
    /// </summary>
    public static class Plots
    {
        /// <summary>
        /// Gets the legend from a ggplot object returned in a variable called 'legend'
        /// </summary>
        /// <returns>The legend.</returns>
        /// <param name="ggplotVar">Ggplot variable.</param>
        public static string GetLegend(string legendVar, string ggplotVar)
        {
            return string.Join("\n", new string[]
            {
                "g_legend<-function(a.gplot){ ",
                "    tmp <- ggplot_gtable(ggplot_build(a.gplot)) ",
                "    leg <- which(sapply(tmp$grobs, function(x) x$name) == \"guide-box\")",
                "    legend <- tmp$grobs[[leg]] ",
                "    return(legend)} ",
                string.Format("{0} <- g_legend({1}) ", legendVar, ggplotVar),
            });
        }

        /// <summary>
        /// Applies theme options and overrides
        /// </summary>
        /// <returns>The options.</returns>
        /// <param name="manualOverrides">Manual overrides.</param>
        /// <param name="description">Description.</param>
        public static string ThemeOptions(Dictionary<string, string> manualOverrides, PlotDescription description = null)
        {
            var basicOptions = new Dictionary<string, string>
            {
                { "panel.background", "element_blank()" },
                { "panel.grid.major", "element_line(color='grey60')" },
                { "panel.grid.minor", "element_blank()" },
                { "axis.text.x", "element_text(color=\"#000000\", size=14)" },
                { "axis.text.y", "element_text(color=\"#000000\", size=14)" },
                { "axis.title.x", "element_text(face=\"bold\", size=18)" },
                { "axis.title.y", "element_text(face=\"bold\", size=18)" },
                { "legend.background", "element_blank()" },
                { "legend.position", "'bottom'" },
                { "legend.direction", "'horizontal'" },
                { "legend.box", "'vertical'" },
                { "legend.title", "element_blank()" },
            };

            Action<Dictionary<string, string>> OverrideOptions = (overrides) =>
            {
                foreach (var v in overrides)
                {
                    if (basicOptions.ContainsKey(v.Key))
                    {
                        if (string.IsNullOrEmpty(v.Value))
                        {
                            basicOptions.Remove(v.Key);
                        }
                        else
                        {
                            basicOptions[v.Key] = v.Value;
                        }
                    }
                    else if (!string.IsNullOrEmpty(v.Value))
                    {
                        basicOptions.Add(v.Key, v.Value);
                    }
                }
            };

            if (description != null)
            {
                OverrideOptions(description.ThemeOverrides);
            }

            OverrideOptions(manualOverrides);

            return "    theme(" + string.Join(", ", basicOptions.Select(x => x.Key + "=" + x.Value)) + ")";
        }

        /// <summary>
        /// Places the legend at the bottom of the plot.
        /// </summary>
        /// <returns>The legend.</returns>
        public static string BottomLegend(PlotDescription description)
        {
            return ThemeOptions(new Dictionary<string, string>
            {
                { "legend.position", "'bottom'" },
                { "legend.direction", "'vertical'" },
                { "legend.box", "'horizontal'" }, 
            }
                .Concat(LegendTitle(description))
                .ToDictionary(x => x.Key, x => x.Value), 
                description);
        }

        /// <summary>
        /// Places the legend at the bottom of the plot.
        /// </summary>
        /// <returns>The legend.</returns>
        public static string BottomLegendHV(PlotDescription description)
        {
            return ThemeOptions(new Dictionary<string, string>
            {
                { "legend.position", "'bottom'" },
                { "legend.direction", "'horizontal'" },
                { "legend.box", "'vertical'" }, 
            }
                .Concat(LegendTitle(description))
                .ToDictionary(x => x.Key, x => x.Value), 
                description);
        }

        /// <summary>
        /// Places the legend at the bottom of the plot.
        /// </summary>
        /// <returns>The legend.</returns>
        public static string BottomLegendHV()
        {
            return "    theme(legend.position=\"bottom\",  legend.direction='horizontal', legend.box='vertical')";
        }

        /// <summary>
        /// Places the legend at the bottom of the plot.
        /// </summary>
        /// <returns>The legend.</returns>
        public static string BottomLegendHorizontal(PlotDescription description)
        {
            return ThemeOptions(new Dictionary<string, string>
            {
                { "legend.position", "'bottom'" },
                { "legend.direction", "'horizontal'" },
                { "legend.box", "'horizontal'" }, 
            }
                .Concat(LegendTitle(description))
                .ToDictionary(x => x.Key, x => x.Value), 
                description);
        }

        /// <summary>
        /// Places the legend at the bottom of the plot.
        /// </summary>
        /// <returns>The legend.</returns>
        public static string BottomLegendHorizontal()
        {
            return "    theme(legend.position=\"bottom\",  legend.direction='horizontal', legend.box='horizontal')";
        }

        /// <summary>
        /// Hides the legend.
        /// </summary>
        /// <returns>The legend.</returns>
        public static string HideLegend(PlotDescription description = null)
        {
            return ThemeOptions(new Dictionary<string, string>
            {
                { "legend.position", "'none'" },
            },
                description);
        }

        /// <summary>
        /// Starts the plot.
        /// </summary>
        /// <returns>The plot.</returns>
        /// <param name="plotFile">Plot file.</param>
        /// <param name="width">Width.</param>
        /// <param name="height">Height.</param>
        public static string StartPlot(string plotFile, double width, double height)
        {
            return string.Format("pdf(file='{0}', width={1}, height={2})", plotFile, width, height);
        }

        /// <summary>
        /// Starts the plot.
        /// </summary>
        /// <returns>The plot.</returns>
        /// <param name="description">Description.</param>
        public static string StartPlot(PlotDescription description)
        {
            return StartPlot(description.PlotFile, description.Width, description.Height);
        }

        /// <summary>
        /// Starts the legend.
        /// </summary>
        /// <returns>The legend.</returns>
        /// <param name="description">Description.</param>
        public static string StartLegend(PlotDescription description)
        {
            return StartPlot(description.LegendFile, description.LegendWidth, description.LegendHeight);
        }

        /// <summary>
        /// Runs a ggplot with the given lines and legend data.
        /// </summary>
        /// <returns>The plot.</returns>
        /// <param name="description">Description.</param>
        /// <param name="plotLines">Plot lines.</param>
        public static string[] GGPlot(PlotDescription description, string[] plotLines, string dataVar = default(string))
        {
            string plotVar = R.Variables.Plot;

            return new string[]
            {
                StartPlot(description.PlotFile, description.Width, description.Height),
                CreateGGPlot(plotVar, dataVar),
                R.JoinLines(plotLines),
                description.SeparateLegend ? HideLegend(description) : FinishGGPlot(description),
                plotVar,
                EndPlot()
            };
        }

        /// <summary>
        /// Creates a ggplot legend in a separate file if the description calls for it.
        /// </summary>
        /// <returns>The legend.</returns>
        /// <param name="description">Description.</param>
        /// <param name="plotLines">Plot lines.</param>
        public static string[] GGLegend(PlotDescription description, string[] plotLines, string dataVar = default(string))
        {
            if (description.LegendFile != null)
            {
                string plotVar = R.Variables.Plot;
                string legendVar = R.Variables.Legend;

                return new string[]
                {
                    StartPlot(description.LegendFile, description.LegendWidth, description.LegendHeight),
                    CreateGGPlot(plotVar, dataVar),
                    R.JoinLines(plotLines),
                    description.LegendDirection == PlotDescription.Direction.Horizontal ? BottomLegendHV(description) : FinishGGPlot(description),

                    GetLegend(legendVar, plotVar),
                    legendVar,
                    string.Format("grid.draw({0})", legendVar),
                    EndPlot()
                };
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// GGs the plot.
        /// </summary>
        /// <returns>The plot.</returns>
        /// <param name="variable">Variable.</param>
        public static string CreateGGPlot(string plotVariable, string dataVariable = default(string))
        {
            if (string.IsNullOrEmpty(dataVariable))
            {
                dataVariable = string.Empty;
            }

            return string.Format("{0} <- ggplot({1}) + ", plotVariable, dataVariable);
        }

        /// <summary>
        /// Labels the axes.
        /// </summary>
        /// <returns>The axes.</returns>
        /// <param name="description">Description.</param>
        public static string LabelAxes(PlotDescription description)
        {
            string labels = string.Empty;

            Action<string, Func<PlotDescription, AxisDescription>> addLabel = (name, axis) =>
            {
                var a = axis(description);
                if (a != null && !string.IsNullOrEmpty(a.Label))
                {
                    labels += LabelAxis(name, a.Label, !a.UnquotedLabel);
                }
            };

            addLabel("x", d => d.XAxis);
            addLabel("y", d => d.YAxis);

            return labels;
        }

        /// <summary>
        /// Color the specified description.
        /// </summary>
        /// <param name="description">Description.</param>
        public static string Color(PlotDescription description)
        {
            return description.Colors == null ? 
                description.RgbColors == null ? 
                    string.Empty :
                    "scale_color_manual(values=c('" + string.Join(", '", description.RgbColors) + ")) + " : 
                "scale_color_manual(values=c('" + string.Join("', '", description.Colors) + "')) + ";
        }

        /// <summary>
        /// Colors the fill.
        /// </summary>
        /// <returns>The fill.</returns>
        /// <param name="description">Description.</param>
        public static string ColorFill(PlotDescription description)
        {
            return description.FillColors == null ? 
                description.FillRgbColors == null ?
                    string.Empty : 
                    "scale_fill_manual(values=c('" + string.Join(", '", description.FillRgbColors) + ")) + " :
                "scale_fill_manual(values=c('" + string.Join("', '", description.FillColors) + "')) + ";
        }

        /// <summary>
        /// Labels the axis.
        /// </summary>
        /// <returns>The axis.</returns>
        /// <param name="axis">Axis.</param>
        /// <param name="label">Label.</param>
        /// <param name="quoted">If set to <c>true</c> quoted.</param>
        private static string LabelAxis(string axis, string label, bool quoted)
        {
            string delim = string.Empty;
            if (quoted)
            {
                delim = "'";
            }

            return string.Format("    {0}lab({1}{2}{3}) + ", axis, delim, label, delim);
        }

        private static List<KeyValuePair<string, string>> LegendTitle(PlotDescription description)
        {
            return string.IsNullOrEmpty(description.FillLabel) ? 
                new List<KeyValuePair<string, string>>() :
                new List<KeyValuePair<string, string>> { new KeyValuePair<string, string>("legend.title", "'" + description.FillLabel + "'") };
        }

        /// <summary>
        /// Finishes the GG plot.
        /// </summary>
        /// <returns>The GG plot.</returns>
        public static string FinishGGPlot(PlotDescription description)
        {
            return ThemeOptions(new Dictionary<string, string>
            {
            }
                .Concat(LegendTitle(description))
                .ToDictionary(x => x.Key, x => x.Value), 
                description);
        }

        /// <summary>
        /// Finishes the GG plot.
        /// </summary>
        /// <returns>The GG plot.</returns>
        public static string FinishGGPlot()
        {
            return ThemeOptions(new Dictionary<string, string>{ });
        }

        /// <summary>
        /// Ends the plot.
        /// </summary>
        /// <returns>The plot.</returns>
        public static string EndPlot()
        {
            return "dev.off()";
        }

        /// <summary>
        /// Create script lines based on the YNames fields
        /// </summary>
        /// <returns>The to Y names.</returns>
        /// <param name="description">Description.</param>
        /// <param name="codeGenerator">Code generator.</param>
        public static string[] ApplyToYNames(
            PlotDescription description, 
            Func<string, int, string> codeGenerator)
        {
            return description.YNames.Select((name, i) => codeGenerator(name, i)).ToArray();
        }

        /// <summary>
        /// Plots the histone profile.
        /// </summary>
        /// <param name="filename">Filename.</param>
        /// <param name="tissue">Cell line.</param>
        /// <param name="expressionType">Expression type.</param>
        /// <param name="histoneType">Histone type.</param>
        public static void PlotHistoneProfile(string filename, string tissue, string expressionType, string histoneType)
        {
            string createPng = 
                "png(file=\"results/HistoneProfile/profiles/" + tissue + "." + expressionType + "." + histoneType + ".png\", width=600, height=400)\n";

            string createPdf = 
                "pdf(file=\"results/HistoneProfile/profiles/" + tissue + "." + expressionType + "." + histoneType + ".pdf\", width=7.5, height=5)\n";

            string load = 
                "x <- read.table(\"" + filename + "\")\n" +
                "ccc <- c(rgb(0,1,0,0.8), rgb(0,1,1,0.8))\n";

            string plot = 
                "ymax <- 0.75; if (max(x[,2:3]) > ymax) { ymax <- 0.75 }\n" +
                "plot(x[,1], x[,2], col=ccc[1], ylim=c(0, ymax), xlim=c(-10000, 10000), xlab=\"Distance from TSS\", ylab=\"Fraction of Peaks\")\n" +
                "points(x[,1], x[,3], col=ccc[2])\n" +
                "title(\"Promoter Distribution of " + histoneType + " peaks in " + tissue + " " + expressionType + " Map\")\n" +
                "legend(\"topleft\", title=\"TSS Set\", c(\"Variantly Expressed\", \"Map\"), fill=ccc)\n" +
                "dev.off()\n";

            string scriptName = "../temp/plotHistoneProfile" + tissue + expressionType + histoneType + ".R";

            using (TextWriter tw = Helpers.CreateStreamWriter(scriptName))
            {
                tw.WriteLine(load + createPng + plot + createPdf + plot);


                tw.WriteLine("q(save=\"no\")");
            }

            R.RunScript(scriptName, "../temp/junk");
        }

        /// <summary>
        /// Plots the distributions.
        /// </summary>
        /// <param name="description">Description.</param>
        public static void Distributions(PlotDescription description)
        {
            description.Validate();

            const string Mids = "mids";
            const string Frequency = "frequency";

            var createHistograms = ApplyToYNames(
                description, 
                (name, i) => Functions.Histogram(
                    R.Variables.M + i, 
                    R.Variables.X, 
                    name, 
                    description.XAxis.Breaks, true));
                
            var createFrequencies = ApplyToYNames(
                description,
                (name, i) => R.BindColumns(
                    "data" + i, 
                    new string[] { R.Variables.M + i, R.Variables.M + i }, 
                    new string[] { Mids, Frequency }));

            var plotCurves = ApplyToYNames(
                description, (name, i) => string.Format(
                    "    geom_line(aes(x={0}, y={1}, colour={2}), {3}) +", 
                    Mids, 
                    Frequency, 
                    GetColor(description, i), 
                    "data" + i));

            string[] plotLines = GGPlot(description, new string[]
            {
                R.JoinLines(plotCurves),
                LabelAxes(description),
            });

            string[] legendLines = GGLegend(description, plotCurves);

            R.WriteAndRunScript(new string[] 
                {
                    R.LoadLibrary(R.Libraries.GGPlot),
                    R.LoadLibrary(R.Libraries.GridExtra),
                    R.ReadTable(R.Variables.X, description.TableFile, true),
                    R.JoinLines(createHistograms),
                    R.JoinLines(createFrequencies),
                    R.JoinLines(plotLines),
                    R.JoinLines(legendLines),
                },
                description.ScriptFile);
        }

        /// <summary>
        /// Densities the specified description.
        /// </summary>
        /// <param name="description">Description.</param>
        public static void Densities(PlotDescription description)
        {
            description.Validate();

            if (description.PlotVariable == null)
            {
                description.PlotVariable = R.Variables.X;
            }

            string densities = string.Format(
                "    geom_density(aes(x={0}, fill=factor({1}), colour=factor({1})), alpha={2}, {3}) +", 
                description.XName,
                description.Factors[0],
                description.Alpha,
                description.PlotVariable);

            var corePlot = new string[]
            {
                densities,
                LabelAxes(description),
                XCoordinateTransform(description.XAxis),
                Color(description),
                ColorFill(description),
            };

            string[] plotLines = GGPlot(description, corePlot);

            string[] legendLines = GGLegend(description, corePlot);

            R.WriteAndRunScript(new string[] 
            {
                R.LoadLibrary(R.Libraries.GGPlot),
                R.LoadLibrary(R.Libraries.GridExtra),
                R.LoadLibrary(R.Libraries.Scales),
                R.ReadTable(R.Variables.X, description.TableFile, true),
                R.JoinLines(plotLines),
                R.JoinLines(legendLines),
            },
                description.ScriptFile);
        }

        /// <summary>
        /// Creates a scatter plot.
        /// </summary>
        /// <param name="description">Description.</param>
        public static void Scatter(PlotDescription description)
        {
            description.Validate();

            if (description.PlotVariable == null)
            {
                description.PlotVariable = R.Variables.X;
            }

            var corePlot = ApplyToYNames(
                               description, (name, i) => string.Format(
                    "    geom_point(aes(x={0}, y={1}, colour={2}{3}), {4}) +", 
                               description.XName, 
                               name, 
                               description.Factors[0], 
                    description.Size != 0 ? ", size=" + description.Size : string.Empty,
                               description.PlotVariable))
                .Concat(new string[]
                {
                    LabelAxes(description),
                    YCoordinateTransform(description.YAxis),
                    XCoordinateTransform(description.XAxis),
                    Color(description),
                })
                .ToArray();

            string[] plotLines = GGPlot(description, corePlot);

            string[] legendLines = GGLegend(description, corePlot);

            R.WriteAndRunScript(new string[] 
            {
                R.LoadLibrary(R.Libraries.GGPlot),
                R.LoadLibrary(R.Libraries.Scales),
                R.LoadLibrary(R.Libraries.GridExtra),
                R.ReadTable(R.Variables.X, description.TableFile, true),
                R.FactorColumn(R.Variables.X, description.Factors[0]),
                R.JoinLines(plotLines),
                R.JoinLines(legendLines),
            },
                description.ScriptFile);
        }

        /// <summary>
        /// Creates a boxplot of a set of values organized by a factor
        /// </summary>
        /// <param name="description">Description.</param>
        /// <param name="violin">If set to <c>true</c> violin.</param>
        /// <param name="hideOutliers">If set to <c>true</c> hide outliers.</param>
        public static void BoxPlot(PlotDescription description, bool violin, bool hideOutliers = false)
        {
            description.Validate();

            if (description.PlotVariable == null)
            {
                description.PlotVariable = R.Variables.X;
            }

            string function = violin ? "geom_violin" : "geom_boxplot";

            var plotLines = new string[]
            {
                description.Factors.Length == 1 ?
                    string.Format(
                        "    {0}(aes({1}, y={2}), width=0.8, {3}{4}) + ", 
                        function, 
                        description.FormatFactor(0), 
                        description.YNames[0], 
                        R.Variables.X,
                        hideOutliers ? ", outlier.shape=NA" : string.Empty) :
                    string.Format(
                        "    {0}(aes({1}, fill={2}, y={3}), width=0.8, {4}{5}) + ",
                        function, 
                        description.FormatFactor(0), 
                    //description.FactorLevels != null && description.FactorLevels.Count >= 2 && description.FactorLevels[1] != null ?
                        description.FormatFactor(1), 
                        description.YNames[0], 
                        description.PlotVariable,
                        hideOutliers ? ", outlier.shape=NA" : string.Empty),
                YCoordinateTransform(description.YAxis),
                XCoordinateTransform(description.XAxis),
                LabelAxes(description),

                description.FlipCoordinates ? "    coord_flip() + " : string.Empty,
                ColorFill(description),
            };

            string[] corePlot = GGPlot(description, plotLines);
            string[] legendLines = GGLegend(description, plotLines);
            

            R.WriteAndRunScript(new string[] 
            {
                R.LoadLibrary(R.Libraries.GGPlot),
                R.LoadLibrary(R.Libraries.Scales),
                R.LoadLibrary(R.Libraries.GridExtra),
                R.ReadTable(R.Variables.X, description.TableFile, true),
                R.FactorColumn(R.Variables.X, description.Factors[0]),
                description.Factors.Length > 1 ? R.FactorColumn(R.Variables.X, description.Factors[1]) : string.Empty,
                R.JoinLines(corePlot),
                description.LegendFile != null ? R.JoinLines(legendLines) : string.Empty,   
            },
                description.ScriptFile);
        }

        public static void LinePointPlot(PlotDescription description)
        {
            description.Validate();

            int pointSize = description.Size != 0 ? description.Size : 2;
            int lineSize = description.Size != 0 ? description.Size / 2 : 1;
            Func<string, string> plot = (function) => description.Factors.Length == 0 ?
                string.Format(
                    "    {0}(aes(x={1}, y={2}), size={3}, {4}) + ", 
                    function, 
                    description.XName, 
                    description.YNames[0], 
                    function == "geom_point" ? pointSize : lineSize,
                    R.Variables.X) :
                string.Format(
                    "    {0}(aes(x={1}, y={2}, color={3}), size={4}, {5}) + ",
                    function,
                    description.XName,
                    description.YNames[0], 
                    description.FormatFactor(0), 
                    function == "geom_point" ? pointSize : lineSize,
                    R.Variables.X);

            var plotLines = new string[]
            {
                plot("geom_line"),
                plot("geom_point"),
                description.YNames.Length == 3 && description.Factors.Length > 0 ?
                    string.Format(
                    "    {0}(aes(x={1}, ymin={2}, ymax={3}, fill={4}), alpha=0.5, {5}) + ",
                        "geom_ribbon",
                        description.XName,
                        description.YNames[1], 
                        description.YNames[2], 
                        description.FormatFactor(0), 
                        R.Variables.X) :
                string.Empty,

                //description.YAxis.Min != 0 && description.YAxis.Max != 0 ?
                //string.Format("coord_cartesian(ylim = c({0}, {1})) + ", description.YAxis.Min, description.YAxis.Max) :
                //string.Empty,
                YCoordinateTransform(description.YAxis),
                XCoordinateTransform(description.XAxis),
                LabelAxes(description),
                ColorFill(description),

                description.FlipCoordinates ? "    coord_flip() + " : string.Empty,
            };

            string[] corePlot = GGPlot(description, plotLines);
            string[] legendLines = GGLegend(description, plotLines);


            R.WriteAndRunScript(new string[] 
            {
                R.LoadLibrary(R.Libraries.GGPlot),
                R.LoadLibrary(R.Libraries.Scales),
                R.LoadLibrary(R.Libraries.GridExtra),
                R.ReadTable(R.Variables.X, description.TableFile, true),
                description.Factors.Length > 0 ? R.FactorColumn(R.Variables.X, description.Factors[0]) : string.Empty,
                description.Factors.Length > 1 ? R.FactorColumn(R.Variables.X, description.Factors[1]) : string.Empty,
                R.JoinLines(corePlot),
                description.LegendFile != null ? R.JoinLines(legendLines) : string.Empty,   
            },
                description.ScriptFile);
        }


        /// <summary>
        /// Creates a boxplot of a set of values organized by a factor
        /// </summary>
        /// <param name="description">Description.</param>
        public static void FactorBarPlot(PlotDescription description, bool whiskers)
        {
            description.Validate();

            if (description.PlotVariable == null)
            {
                description.PlotVariable = R.Variables.X;
            }

            var plotLines = new string[]
            {
                description.Factors.Length == 1 ?
                string.Format(
                    "    geom_bar(aes(x=factor({0}), y={1}), position='dodge', stat='identity', width=0.8, {2}) + ", 
                    description.Factors[0], 
                    description.YNames[0], 
                    description.PlotVariable) :
                string.Format(
                    "    geom_bar(aes(x={0}, fill={1}, y={2}), position='dodge', stat='identity', width=0.6, {3}) + ",
                    description.FormatFactor(0), 
                    description.FormatFactor(1), 
                    description.YNames[0], 
                    description.PlotVariable),
                whiskers ? string.Format(
                    "    geom_errorbar(aes(x={0}, ymin={1}-{2}, ymax={1}+{2}, fill={3}), position='dodge', stat='identity', width=0.6, {4}) + ",
                    description.FormatFactor(0), 
                    description.YNames[0], 
                    description.YNames[1], 
                    description.FormatFactor(1), 
                    description.PlotVariable) : string.Empty,
                LabelAxes(description),

                YCoordinateTransform(description.YAxis),

                ColorFill(description),
            };

            string[] corePlot = GGPlot(description, plotLines);
            string[] coreLegend = GGLegend(description, plotLines);

            R.WriteAndRunScript(new string[] 
            {
                R.LoadLibrary(R.Libraries.GGPlot),
                R.LoadLibrary(R.Libraries.GridExtra),
                R.LoadLibrary(R.Libraries.Scales),
                R.ReadTable(R.Variables.X, description.TableFile, true),
                R.FactorColumn(R.Variables.X, description.Factors[0]),
                description.Factors.Length > 1 ? R.FactorColumn(R.Variables.X, description.Factors[1]) : string.Empty,
                R.JoinLines(corePlot),
                description.LegendFile != null ? R.JoinLines(coreLegend) : string.Empty,   
            },
                description.ScriptFile);
        }

        /// <summary>
        /// Creates a CDF of the given X value organized by factor
        /// </summary>
        /// <param name="description">Description.</param>
        public static void FactorCDFPlot(PlotDescription description)
        {
            description.Validate();

            string[] corePlot = GGPlot(description, new string[]
            {
                description.Factors.Length == 1 ?
                string.Format(
                    "    stat_ecdf(aes(x={0}, color={1}), {2}) + ", 
                    description.XName, 
                    description.Factors[0], 
                    R.Variables.X) :
                LabelAxes(description),
            });

            R.WriteAndRunScript(new string[] 
            {
                R.LoadLibrary(R.Libraries.GGPlot),
                R.LoadLibrary(R.Libraries.GridExtra),
                R.ReadTable(R.Variables.X, description.TableFile, true),
                R.FactorColumn(R.Variables.X, description.Factors[0]),
                R.JoinLines(corePlot),
            },
                description.ScriptFile);
        }

        /// <summary>
        /// Gets the color.
        /// </summary>
        /// <returns>The color.</returns>
        /// <param name="description">Description.</param>
        /// <param name="index">Index.</param>
        private static string GetColor(PlotDescription description, int index)
        {
            if (description.ColorNames == null)
            {
                return string.Format("'{0}'", index + 1);
            }
            else
            {
                return description.ColorNames[index];
            }
        }

        private static string CoordinateTransform(AxisDescription axis, string axis_name)
        {
            if (axis.Breaks == null && axis.DiscreteBreaks == null && axis.BreakLabels == null &&
                axis.CoordinateTransform == AxisDescription.CoordTransform.None)
            {
                return string.Empty;
            }

            var breakData = string.Empty;
            Func<string[], string[], string> FormatBreaks = (breaks, labels) => labels != null ?
                    string.Format("breaks = c({0}), labels = c({1})",
                        string.Join(", ", breaks),
                        string.Join(", ", labels)) :
                    string.Format("breaks = c({0})",
                        string.Join(", ", breaks));

            if (axis.Breaks != null)
            {
                breakData = FormatBreaks(axis.Breaks.Select(x => x.ToString()).ToArray(), axis.BreakLabels);
            }
            else if (axis.DiscreteBreaks != null)
            {
                breakData = FormatBreaks(axis.DiscreteBreaks, axis.BreakLabels);
            }
            else if (axis.BreakLabels != null)
            {
                breakData = string.Format("breaks = c('{0}')",
                    string.Join("', '", axis.BreakLabels));
            }

            var transformData = axis.CoordinateTransform == AxisDescription.CoordTransform.None ||
                axis.CoordinateTransform == AxisDescription.CoordTransform.Discrete ? 
                    string.Empty :
                    string.Format("trans=log{0}_trans()", axis.CoordinateTransformBase);

            var limitData = (axis.Min != 0 || axis.Max != 0) && axis.Max > axis.Min ?
                string.Format("limits=c({0}, {1})", axis.Min, axis.Max) :
                string.Empty;

            return string.Format(
                "    scale_{0}_{1}({2}) + ", 
                axis_name,
                axis.CoordinateTransform == AxisDescription.CoordTransform.Discrete ? "discrete" : "continuous",
                string.Join(", ", new string[]
                {
                    transformData,
                    breakData,
                    limitData,
                }.Where(x => !string.IsNullOrEmpty(x))));
        }

        public static string YCoordinateTransform(AxisDescription axis)
        {
            return CoordinateTransform(axis, "y");
        }

        public static string XCoordinateTransform(AxisDescription axis)
        {
            return CoordinateTransform(axis, "x");
        }
    }
}

