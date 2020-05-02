//--------------------------------------------------------------------------------
// <copyright file="Functions.cs" 
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
    /// Function calls in R.
    /// </summary>
    public static class Functions
    {
        /// <summary>
        /// Histogram the specified outputVariable, inputVariable and breaks.
        /// </summary>
        /// <param name="outputVariable">Output variable.</param>
        /// <param name="inputVariable">Input variable.</param>
        /// <param name="breaks">Breaks.</param>
        public static string Histogram(
            string outputVariable,
            string inputVariable,
            string inputColumn,
            double[] breaks,
            bool setFrequency)
        {
            string histCommand = string.Format(
                                     "{0} <- hist({1}${2}, breaks=c({3})); ", 
                                     outputVariable,
                                     inputVariable,
                                     inputColumn,
                                     R.JoinEnumerableCsv(breaks));

            if (setFrequency)
            {
                histCommand += string.Format("\n{0}$frequency = {0}$counts/sum({0}$counts);\n", outputVariable);
            }

            return histCommand;
        }
    }
}

