using System;

namespace Genomics
{
    using System.IO;
    using System.Linq;
    using Shared;
    using System.IO.Compression;

    public class FastqParser<AlphabetType> where AlphabetType : Alphabet
    {
        public FastqParser()
        {
        }

        public static StreamReader GetReaderForFile(string filename)
        {
            if (Path.GetExtension(filename) == ".gz")
            {
                return new StreamReader(new GZipStream(
                    new FileStream(filename, FileMode.Open), 
                    CompressionMode.Decompress));
            }
            else
            {
                return new StreamReader(filename);
            }
        }

        private static Alphabet alphabet;

        /// <summary>
        /// Gets the alphabet.
        /// </summary>
        /// <value>The alphabet.</value>
        public static Alphabet Alphabet
        {
            get
            {
                return Helpers.CheckInit(
                    ref alphabet,
                    () => Alphabets.Get(typeof(AlphabetType)));
            }
        }

        public static FastqEntry<AlphabetType> FastqEntryFromReader(TextReader tr)
        {

            if (tr.Peek() == -1)
            {
                throw new MissingHeader();
            }

            var header = tr.ReadLine();

            if (header.Length < 2 || header.First() != '@')
            {
                throw new HeaderFormat();
            }

            if (tr.Peek() == -1)
            {
                throw new MissingSequence();
            }

            var sequence = tr.ReadLine();

            if (!Alphabet.ValidateSequence(sequence))
            {
                throw new SequenceFormat();
            }

            if (tr.Peek() == -1)
            {
                throw new MissingDelimiter();
            }

            var delimiter = tr.ReadLine();

            if (delimiter != "+")
            {
                throw new InvalidDelimiter();
            }

            if (tr.Peek() == -1)
            {
                throw new MissingQuality();
            }

            var quality = tr.ReadLine();

            if (!Alphabets.QualityAlphabet.ValidateSequence(quality))
            {
                throw new QualityFormat();
            }

            if (quality.Length != sequence.Length)
            {
                throw new LengthMismatch();
            }

            var entry = new FastqEntry<AlphabetType>
            {
                Id = header,
                Sequence = sequence,
                Quality = quality,
            };

            return entry;
        }


        public class FastqParserException : Exception
        {
        }

        public class MissingHeader : FastqParserException
        {
        }

        public class MissingSequence : FastqParserException
        {
        }

        public class MissingDelimiter : FastqParserException
        {
        }

        public class InvalidDelimiter : FastqParserException
        {
        }

        public class MissingQuality : FastqParserException
        {
        }

        public class LengthMismatch : FastqParserException
        {
        }

        public class HeaderFormat : FastqParserException
        {
        }

        public class SequenceFormat : FastqParserException
        {
        }

        public class QualityFormat : FastqParserException
        {
        }
    }
}

