//--------------------------------------------------------------------------------
// <copyright file="PlotDescription.cs" 
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

    /// <summary>
    /// Manages description of a plot.
    /// </summary>
        public class PlotDescription
    {
        private Dictionary<string, string> themeOverrides;

        /// <summary>
        /// Plot file type.
        /// </summary>
        public enum PlotFileType
        {
            /// <summary>
            /// Pdf output type.
            /// </summary>
            Pdf,

            /// <summary>
            /// Eps output type.
            /// </summary>
            Eps,
        }

        /// <summary>
        /// Gets or sets the table file underlying the plot.
        /// </summary>
        /// <value>The table file.</value>
        public string TableFile { get; set; }

        /// <summary>
        /// Gets or sets the script file.
        /// </summary>
        /// <value>The script file.</value>
        public string ScriptFile { get; set; }

        /// <summary>
        /// Gets or sets the plot output file.
        /// </summary>
        /// <value>The plot file.</value>
        public string PlotFile { get; set; }

        /// <summary>
        /// Gets or sets the plot file width (in inches).
        /// </summary>
        /// <value>The width.</value>
        public double Width { get; set; }

        /// <summary>
        /// Gets or sets the plot file height (in inches).
        /// </summary>
        /// <value>The height.</value>
        public double Height { get; set; }

        /// <summary>
        /// Gets the type of the file.
        /// </summary>
        /// <value>The type of the file.</value>
        public PlotFileType FileType
        { 
            get
            {
                return PlotFileType.Pdf;
            }
        }

        /// <summary>
        /// Gets or sets the size of the point/line, etc
        /// </summary>
        /// <value>The size.</value>
        public int Size { get; set; }

        /// <summary>
        /// Gets or sets the name of the table column to use as the x values.
        /// </summary>
        /// <value>The name of the X.</value>
        public string XName { get; set; }

        /// <summary>
        /// Gets or sets the names of the table columns to use as the y values.
        /// </summary>
        /// <value>The Y names.</value>
        public string[] YNames { get; set; }

        /// <summary>
        /// Data columns to turn into and/or reference as a factors.
        /// </summary>
        /// <value>The factor.</value>
        public string[] Factors { get; set; }

        public string[][] FactorLevels { get; set; }

        public string[][] FactorLabels { get; set; }

        public string FormatFactor(int i)
        {
            if (!this.HasFactorLevels(i) && !this.HasFactorLabels(i))
            {
                return string.Format("factor({0})", this.Factors[i]);
            }
            else if (this.HasFactorLevels(i) && !this.HasFactorLabels(i))
            {
                return string.Format(
                    "factor({0}, ordered=T, levels=c({1}))", 
                    this.Factors[i], 
                    "'" + string.Join("', '", this.FactorLevels[i]) + "'");
            }
            else if (this.HasFactorLevels(i) && this.HasFactorLabels(i))
            {
                return string.Format(
                    "factor({0}, ordered=T, levels=c({1}), labels=c({2}))",
                    this.Factors[i], 
                    "'" + string.Join("', '", this.FactorLevels[i]) + "'",
                    "'" + string.Join("', '", this.FactorLabels[i]) + "'");

            }
            else
            {
                throw new Exception(string.Format("Factor labels provided without factor levels for factor {0}", this.Factors[i]));
            }
        }

        private bool HasFactorLevels(int i)
        {
            return this.FactorLevels != null && i < this.FactorLevels.Length && this.FactorLevels[i] != null;
        }

        private bool HasFactorLabels(int i)
        {
            return this.FactorLabels != null && i < this.FactorLabels.Length && this.FactorLabels[i] != null;
        }

        /// <summary>
        /// Gets or sets the fill colors.
        /// </summary>
        /// <value>The fill colors.</value>
        public string[] FillColors { get; set; }

        /// <summary>
        /// Gets or sets the fill rgb colors.
        /// </summary>
        /// <value>The fill rgb colors.</value>
        public string[] FillRgbColors { get; set; }

        /// <summary>
        /// Gets or sets the colors.
        /// </summary>
        /// <value>The colors.</value>
        public string[] Colors { get; set; }

        /// <summary>
        /// Gets or sets the rgb colors.
        /// </summary>
        /// <value>The rgb colors.</value>
        public string[] RgbColors { get; set; }

        /// <summary>
        /// Gets or sets the names of the table columns to use as the y error values.
        /// Must correspond to YNames.
        /// </summary>
        /// <value>The Y errors.</value>
        public string[] YErrors { get; set; }

        /// <summary>
        /// Gets or sets the color names.
        /// </summary>
        /// <value>The color names.</value>
        public string[] ColorNames { get; set; }

        /// <summary>
        /// Gets or sets the fill names.
        /// </summary>
        /// <value>The fill names.</value>
        public string[] FillNames { get; set; }

        /// <summary>
        /// Gets or sets the label of the color legend. 
        /// </summary>
        /// <value>The color label.</value>
        public string ColorLabel { get; set; }

        /// <summary>
        /// Gets or sets the label of the fill legend.
        /// </summary>
        /// <value>The fill label.</value>
        public string FillLabel { get; set; }

        /// <summary>
        /// Gets or sets the legend plot file (optional).
        /// </summary>
        /// <value>The legend file.</value>
        public string LegendFile { get; set; }

        /// <summary>
        /// Gets or sets the width of the legend file (required if using a legend file).
        /// </summary>
        /// <value>The width of the legend.</value>
        public double LegendWidth { get; set; }

        /// <summary>
        /// Gets or sets the height of the legend (required if using a legend file).
        /// </summary>
        /// <value>The height of the legend.</value>
        public double LegendHeight { get; set; }

        public string PlotVariable { get; set; }

        public enum Direction 
        {
            Horizontal,
            Vertical,
        }

        public Direction LegendDirection { get; set; }

        /// <summary>
        /// Gets a value indicating whether this <see cref="Data.PlotDescription"/> separate legend.
        /// </summary>
        /// <value><c>true</c> if separate legend; otherwise, <c>false</c>.</value>
        public bool SeparateLegend
        {
            get
            {
                return !string.IsNullOrEmpty(this.LegendFile);
            }
        }

        /// <summary>
        /// Gets or sets the X axis.
        /// </summary>
        /// <value>The X axis.</value>
        public AxisDescription XAxis { get; set; }

        /// <summary>
        /// Gets or sets the Y axis.
        /// </summary>
        /// <value>The Y axis.</value>
        public AxisDescription YAxis { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Data.PlotDescription"/> flip coordinates.
        /// </summary>
        /// <value><c>true</c> if flip coordinates; otherwise, <c>false</c>.</value>
        public bool FlipCoordinates { get; set; }

        /// <summary>
        /// Gets or sets the alpha.
        /// </summary>
        /// <value>The alpha.</value>
        public double Alpha { get; set; }

        /// <summary>
        /// Gets or sets the theme overrides.
        /// </summary>
        /// <value>The theme overrides.</value>
        public Dictionary<string, string> ThemeOverrides
        {
            get
            {
                return this.themeOverrides ?? new Dictionary<string, string>();
            }

            set
            {
                this.themeOverrides = value;
            }
        }

        /// <summary>
        /// Validate this instance.
        /// </summary>
        public void Validate()
        {
            this.ValidateField(this.TableFile, "table file");
            this.ValidateField(this.ScriptFile, "script file");
            this.ValidateField(this.PlotFile, "plot file");
            this.ValidateField(this.Width, "plot width");
            this.ValidateField(this.Height, "plot height");

            this.ValidateField(this.XAxis, "x-axis");
            this.ValidateField(this.YAxis, "y-axis");


            if (this.SeparateLegend)
            {
                this.ValidateField(this.LegendWidth, "legend width");
                this.ValidateField(this.LegendWidth, "legend height");
            }
        }

        /// <summary>
        /// Validates the field.
        /// </summary>
        /// <param name="field">Field.</param>
        /// <param name="name">Name.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public void ValidateField<T>(T field, string name)
        {
            if (field == null)
            {
                throw new Exception("Error: missing plot description field: " + name);
            }
        }

        /// <summary>
        /// Validates the field.
        /// </summary>
        /// <param name="field">Field.</param>
        /// <param name="name">Name.</param>
        public void ValidateField(double field, string name)
        {
            if (field == 0.0)
            {
                throw new Exception("Error: invalid zero value in plot description field: " + name);

            }
        }
    }
}