using System;

namespace Genomics
{
    public class FastqEntry<AlphabetType> where AlphabetType : Alphabet
    {
        public FastqEntry()
        {
        }

        /// <summary>
        /// Gets the identifier.
        /// </summary>
        /// <value>The identifier.</value>
        public string Id { get; set; }

        /// <summary>
        /// Gets the sequence.
        /// </summary>
        /// <value>The sequence.</value>
        public string Sequence { get; set; }

        /// <summary>
        /// Gets the quality.
        /// </summary>
        /// <value>The quality.</value>
        public string Quality { get; set; }
    }
}

