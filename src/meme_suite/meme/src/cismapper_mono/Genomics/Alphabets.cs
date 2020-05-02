
namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Shared;

    /// <summary>
    /// Alphabet.
    /// </summary>
    public abstract class Alphabet
    {
        protected Alphabet(char[] letters, uint[] values)
        {
            this.Letters = letters;
            this.Values = values;
        }

        /// <summary>
        /// Gets the letters.
        /// </summary>
        /// <value>The letters.</value>
        public char[] Letters { get; private set; }

        /// <summary>
        /// Gets the values.
        /// </summary>
        /// <value>The values.</value>
        public uint[] Values { get; private set; }

        /// <summary>
        /// The letter values.
        /// </summary>
        protected static Dictionary<char, uint> letterValues;

        /// <summary>
        /// Gets the letter value.
        /// </summary>
        /// <value>The letter value.</value>
        public Dictionary<char, uint> LetterValue
        {
            get
            {
                return Helpers.CheckInit(
                    ref letterValues,
                    () => this.Letters
                    .Select((x, i) => new { Key = x, Value = this.Values[i] })
                    .ToDictionary(x => x.Key, x => x.Value));
            }
        }

        /// <summary>
        /// Gets the letter set.
        /// </summary>
        /// <value>The letter set.</value>
        public HashSet<char> LetterSet
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.letterSet,
                    () => new HashSet<char>(this.Letters));
            }
        }

        /// <summary>
        /// Validates the sequence.
        /// </summary>
        /// <returns><c>true</c>, if sequence was validated, <c>false</c> otherwise.</returns>
        /// <param name="sequence">Sequence.</param>
        public bool ValidateSequence(string sequence)
        {
            return sequence.All(this.LetterSet.Contains);
        }

        /// <summary>
        /// The letter set.
        /// </summary>
        protected HashSet<char> letterSet;
    }

    /// <summary>
    /// Nucelotide alphabet.
    /// </summary>
    public abstract class NucleotideAlphabet : Alphabet
    {
        /// <summary>
        /// The reverse complement letter.
        /// </summary>
        private Dictionary<char, char> reverseComplementLetter;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.NucleotideAlphabet"/> class.
        /// </summary>
        /// <param name="letters">Letters.</param>
        /// <param name="reverseComplementLetters">Reverse complement letters.</param>
        /// <param name="values">Values.</param>
        protected NucleotideAlphabet(char[] letters, char[] reverseComplementLetters, uint[] values)
            : base(letters, values)
        {
            this.ReverseComplementLetters = reverseComplementLetters;
        }

        /// <summary>
        /// Reverses the complement.
        /// </summary>
        /// <returns>The complement.</returns>
        /// <param name="sequence">Sequence.</param>
        public string ReverseComplement(string sequence)
        {
            return string.Join("", sequence.Reverse().Select(x => this.ReverseComplementLetter[x]));
        }

        /// <summary>
        /// Gets the reverse complement letters.
        /// </summary>
        /// <value>The reverse complement letters.</value>
        public char[] ReverseComplementLetters { get; private set; }

        /// <summary>
        /// Gets the reverse complement letter.
        /// </summary>
        /// <value>The reverse complement letter.</value>
        public Dictionary<char, char> ReverseComplementLetter
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.reverseComplementLetter,
                    () => this.Letters
                        .Select((x, i) => new { Key = x, Value = this.ReverseComplementLetters[i] })
                        .ToDictionary(x => x.Key, x => x.Value));
            }
        }
    }

    /// <summary>
    /// Alphabets.
    /// </summary>
    public static class Alphabets
    {
        /// <summary>
        /// Get the specified type.
        /// </summary>
        /// <param name="type">Type.</param>
        public static Alphabet Get(Type type)
        {
            //if (Enum.Parse(Alphabets.Type, type
            return alphabets[type];
        }

        /// <summary>
        /// Gets the quality alphabet.
        /// </summary>
        /// <value>The quality alphabet.</value>
        public static Alphabet QualityAlphabet
        {
            get
            {
                return Get(typeof(Genomics.QualityAlphabet));
            }
        }

        /// <summary>
        /// The nucleotide alphabets.
        /// </summary>
        private static readonly Dictionary<Type, Alphabet> alphabets = new Dictionary<Type, Alphabet>
        {
            { typeof(DnaAlphabet), new DnaAlphabet() },
            { typeof(DnaAmbiguousAlphabet), new DnaAmbiguousAlphabet() },
            { typeof(RnaAlphabet), new RnaAlphabet() },
            { typeof(QualityAlphabet), new QualityAlphabet() },
        };
    }

    public static class NucleotideAlphabets
    {
        /// <summary>
        /// Get the specified type.
        /// </summary>
        /// <param name="type">Type.</param>
        public static NucleotideAlphabet Get(Type type)
        {
            //if (Enum.Parse(Alphabets.Type, type
            return nucleotideAlphabets[type];
        }


        /// <summary>
        /// The alphabets.
        /// </summary>
        private static readonly Dictionary<Type, NucleotideAlphabet> nucleotideAlphabets = new Dictionary<Type, NucleotideAlphabet>
        {
            { typeof(DnaAlphabet), new DnaAlphabet() },
            { typeof(DnaAmbiguousAlphabet), new DnaAmbiguousAlphabet() },
            { typeof(RnaAlphabet), new RnaAlphabet() },
        };

    }

    /// <summary>
    /// Rna alphabet.
    /// </summary>
    public class RnaAlphabet : NucleotideAlphabet
    {
        /// <summary>
        /// The letters.
        /// </summary>
        private static readonly char[] letters = { 'A', 'C', 'G', 'U' };

        /// <summary>
        /// The reverse complement letters.
        /// </summary>
        private static readonly char[] reverseComplementLetters = { 'U', 'G', 'C', 'A' };

        /// <summary>
        /// The values.
        /// </summary>
        private static readonly uint[] values = { 1, 2, 4, 8, };

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.RnaAlphabet"/> class.
        /// </summary>
        public RnaAlphabet()
            : base(letters, reverseComplementLetters, values)
        {
        }
    }

    /// <summary>
    /// Dna alphabet.
    /// </summary>
    public class DnaAlphabet : NucleotideAlphabet
    {
        /// <summary>
        /// The letters.
        /// </summary>
        private static readonly char[] letters = { 'A', 'C', 'G', 'T' };

        /// <summary>
        /// Sets the reverse complement letters.
        /// </summary>
        /// <value>The reverse complement letters.</value>
        private static readonly char[] reverseComplementLetters = { 'T', 'G', 'C', 'A' };

        /// <summary>
        /// The values.
        /// </summary>
        private static readonly uint[] values = { 1, 2, 4, 8, };

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.DnaAlphabet"/> class.
        /// </summary>
        public DnaAlphabet()
            : base(letters, reverseComplementLetters, values)
        {
        }
    }

    /// <summary>
    /// Dna alphabet.
    /// </summary>
    public class DnaAmbiguousAlphabet : NucleotideAlphabet
    {
        /// <summary>
        /// The letters.
        /// </summary>
        private static readonly char[] letters = { 'A', 'C', 'G', 'T', 'N' };

        /// <summary>
        /// Sets the reverse complement letters.
        /// </summary>
        /// <value>The reverse complement letters.</value>
        private static readonly char[] reverseComplementLetters = { 'T', 'G', 'C', 'A', 'A' };

        /// <summary>
        /// The values.
        /// </summary>
        private static readonly uint[] values = { 1, 2, 4, 8, 15 };

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.DnaAmbiguousAlphabet"/> class.
        /// </summary>
        public DnaAmbiguousAlphabet()
            : base(letters, reverseComplementLetters, values)
        {
        }
    }

    /// <summary>
    /// Quality alphabet.
    /// </summary>
    public class QualityAlphabet : Alphabet
    {
        /// <summary>
        /// The letters.
        /// </summary>
        private static readonly char[] letters = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~".ToArray();

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.QualityAlphabet"/> class.
        /// </summary>
        public QualityAlphabet()
            : base(letters, letters.Select((x, i) => (uint)i).ToArray())
        {
        }
    }
}

