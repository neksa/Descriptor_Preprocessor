//--------------------------------------------------------------------------------
// <copyright file="AxisDescription.cs" 
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
    /// Describes a plot axis.
    /// </summary>
    public class AxisDescription
    {
        /// <summary>
        /// Gets or sets the minimum.
        /// </summary>
        /// <value>The minimum.</value>
        public double Min { get; set; }

        /// <summary>
        /// Gets or sets the max.
        /// </summary>
        /// <value>The max.</value>
        public double Max { get; set; }

        /// <summary>
        /// Gets or sets the label.
        /// </summary>
        /// <value>The label.</value>
        public string Label { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether this <see cref="Data.AxisDescription"/> unquoted label.
        /// </summary>
        /// <value><c>true</c> if unquoted label; otherwise, <c>false</c>.</value>
        public bool UnquotedLabel { get; set; }

        /// <summary>
        /// Gets or sets the breaks.
        /// </summary>
        /// <value>The breaks.</value>
        public double[] Breaks { get; set; }

        /// <summary>
        /// Gets or sets the discrete breaks.
        /// </summary>
        /// <value>The discrete breaks.</value>
        public string[] DiscreteBreaks { get; set; }

        /// <summary>
        /// Gets or sets the break labels.
        /// </summary>
        /// <value>The break labels.</value>
        public string[] BreakLabels { get; set; }

        public enum CoordTransform
        {
            None,
            Log2,
            Log10,
            Discrete,
        }

        private Dictionary<CoordTransform, double> coordTransformBase = new Dictionary<CoordTransform, double>
        {
            { CoordTransform.None, double.NaN },
            { CoordTransform.Discrete, double.NaN },
            { CoordTransform.Log2, 2 },
            { CoordTransform.Log10, 10 },
        };

        public CoordTransform CoordinateTransform { get; set; }

        public double CoordinateTransformBase
        {
            get
            {
                return this.coordTransformBase[this.CoordinateTransform];
            }
        }
    }
}