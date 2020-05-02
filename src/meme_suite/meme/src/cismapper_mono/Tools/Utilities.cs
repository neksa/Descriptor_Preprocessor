//--------------------------------------------------------------------------------
// File: Utilities.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Tools
{
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;

    /// <summary>
    /// Static data processing methods
    /// </summary>
    public static class Utilities
    {
		/// <summary>
		/// Paste together cell type and expression type used for generating a map or regression score
		/// </summary>
		/// <returns>The source.</returns>
		/// <param name="tissue">Cell line.</param>
		/// <param name="expressionType">Expression type.</param>
		public static string ScoreSource(string tissue, string expressionType)
		{
			return tissue + "." + expressionType;
		}

        /// <summary>
        /// TFs the score source.
        /// </summary>
        /// <returns>The score source.</returns>
        /// <param name="tfName">Tf name.</param>
        /// <param name="tissue">Cell line.</param>
        /// <param name="expressionType">Expression type.</param>
        public static string TFScoreSource(string tfName, string tissue, string expressionType)
        {
            return string.Join(".", new List<string> { tfName, tissue, expressionType });
        }

		/// <summary>
		/// Gets the name of the drm score file.
		/// </summary>
		/// <returns>The drm score file name.</returns>
		/// <param name="sourceFileName">Source file name.</param>
        /// <param name="scoreSource">the score source type.</param>
        public static string GetLocicoreFileName(string sourceFileName, string scoreSource)
		{
            return GetScoreFileName("Locicores." + scoreSource, sourceFileName.Split('/').Last().Replace(".", "_"));
		}

		/// <summary>
		/// Gets the name of the tss score file.
		/// </summary>
		/// <returns>The tss score file name.</returns>
		/// <param name="sourceFileName">Source file name.</param>
        /// <param name="scoreSource">the score source type.</param>
        public static string GetTssScoreFileName(string sourceFileName, string scoreSource)
		{
            return GetScoreFileName("TssScores." + scoreSource, sourceFileName.Split('/').Last().Replace(".", "_"));
		}

		/// <summary>
		/// Gets the name of the score file.
		/// </summary>
		/// <returns>The score file name.</returns>
		/// <param name="scoreType">Score type.</param>
		/// <param name="sourceFileName">Source file name.</param>
		public static string GetScoreFileName(string scoreType, string sourceFileName)
		{
			return string.Format("../temp/GenomicFeatures/{0}.{1}.tsv", scoreType, sourceFileName);
		}

		/// <summary>
		/// Gets the ChIP data file names.
		/// </summary>
		/// <returns>The ch IP data file names.</returns>
		/// <param name="chipFileName">Chip file name.</param>
		public static IEnumerable<string> GetChIPDataFileNames(string chipFileName)
		{  
			using (TextReader tr = new StreamReader(chipFileName))
			{
				return tr.ReadToEnd()
						.Split('\n')
						.Where(line => !string.IsNullOrWhiteSpace(line))
						.Select(line => line.Trim());
			}
		}

		/// <summary>
		/// Gets the ChIP data set names.
		/// </summary>
		/// <returns>The ch IP data set names.</returns>
		/// <param name="chipFileName">Chip file name.</param>
		public static IEnumerable<string> GetChIPDataSetNames(string chipFileName)
		{  
            return GetChIPDataFileNames(chipFileName).Select(x => Path.GetFileNameWithoutExtension(x).Split('/').Last().Replace(".", "_"));
		}
    }
}

