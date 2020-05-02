//--------------------------------------------------------------------------------
// <copyright file="FActorDescription.cs" 
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

    /// <summary>
    /// Factor description.
    /// </summary>
    public class FactorDescription
    {
        /// <summary>
        /// Gets or sets the levels.
        /// </summary>
        /// <value>The levels.</value>
        public string[] Levels { get; set; }

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

        /// <summary>
        /// Gets or sets the labels.
        /// </summary>
        /// <value>The labels.</value>
        public string[] Labels { get; set; }

        /// <summary>
        /// Gets or sets the colors.
        /// </summary>
        /// <value>The colors.</value>
        public string[] Colors { get; set; }

        /// <summary>
        /// Gets the color map.
        /// </summary>
        /// <value>The color map.</value>
        public string[] ColorMap 
        {
            get
            {
                return this.Labels.Select((x, i) => string.Format("{0}'='{1}", x, this.Colors[i])).ToArray();
            }
        }

        /// <summary>
        /// Gets the color rgb map.
        /// </summary>
        /// <value>The color rgb map.</value>
        public string[] ColorRgbMap
        {
            get
            {
                return this.Labels.Select((x, i) => string.Format("{0}'={1}", x, this.Colors[i])).ToArray();
            }
        }
    }
}

