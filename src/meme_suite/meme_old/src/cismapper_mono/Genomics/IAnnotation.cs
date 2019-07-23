//--------------------------------------------------------------------------------
// File: IAnnotation.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System.Collections.Generic;

	/// <summary>
	/// Annotation set data
	/// </summary>
	abstract public class IAnnotation
    {
        /// <summary>
        /// Gets or sets the valid gene types.
        /// </summary>
        /// <value>The valid gene types.</value>
	abstract public string[] ValidTranscriptTypes { get; set; }

	/// <summary>
	/// Gets the valid transcripts.
	/// </summary>
	/// <value>The valid transcripts.</value>
	abstract public HashSet<string> ValidTranscripts { get; }

        /// <summary>
        /// Gets the valid genes.
        /// </summary>
        /// <value>The valid genes.</value>
        abstract public HashSet<string> ValidGenes { get; }

        /// <summary>
        /// Gets the ordered gene locations.
        /// </summary>
        /// <returns>The ordered gene locations.</returns>
        /// <param name="types">Types.</param>
	abstract public Dictionary<string, List<Genomics.Location>> OrderedTranscriptLocations { get; }

        /// <summary>
        /// Gets the indexed transcript locations.
        /// </summary>
        /// <value>The indexed transcript locations.</value>
	abstract public Dictionary<string, Dictionary<int, List<Genomics.Location>>> IndexedTranscriptLocations { get; }

        /// <summary>
        /// Gets maximum feature size
        /// </summary>
        /// <value>The size of the max feature.</value>
	abstract public int MaxFeatureSize { get; }

        /// <summary>
        /// Gets the bin size used for indexing
        /// </summary>
        /// <value>The size of the index.</value>
	abstract public int IndexSize { get; }
    }
}
