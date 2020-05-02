//--------------------------------------------------------------------------------
// File: TRFScorer.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Shared
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.IO;
    using System.Xml.Linq;

    /// <summary>
    /// Static functions for general use
    /// </summary>
    public static class Helpers
    {
        /// <summary>
        /// Creates and instance of a stream writer and ensures the existence of the write path
        /// </summary>
        /// <param name="file"></param>
        /// <returns></returns>
        public static StreamWriter CreateStreamWriter(string file, bool append = false)
        {
            EnsureFilePathExists(file);

            return new StreamWriter(file, append);
        }

        /// <summary>
        /// Appends the stream writer.
        /// </summary>
        /// <returns>The stream writer.</returns>
        /// <param name="file">File.</param>
        public static StreamWriter AppendStreamWriter(string file)
        {
            EnsureFilePathExists(file);

            return File.AppendText(file);
        }

        /// <summary>
        /// Converts the enum.
        /// </summary>
        /// <returns>The enum.</returns>
        /// <param name="eFrom">E from.</param>
        /// <typeparam name="TFrom">The 1st type parameter.</typeparam>
        /// <typeparam name="TTo">The 2nd type parameter.</typeparam>
        public static TTo ConvertEnum<TFrom, TTo>(TFrom eFrom)
        {
            return (TTo)Enum.Parse(typeof(TTo), eFrom.ToString());
        }

        /// <summary>
        /// Enums the is subset.
        /// </summary>
        /// <returns><c>true</c>, if is subset was enumed, <c>false</c> otherwise.</returns>
        /// <typeparam name="TSuperset">The 1st type parameter.</typeparam>
        /// <typeparam name="TSubset">The 2nd type parameter.</typeparam>
        public static bool EnumIsSubset<TSuperset, TSubset>() 
            where TSuperset : struct, IConvertible 
            where TSubset : struct, IConvertible
        {
            EnumTypeException<TSuperset>.VerifyEnum();
            EnumTypeException<TSubset>.VerifyEnum();

            foreach (var value in Enum.GetValues(typeof(TSubset)).Cast<TSubset>())
            {
                TSuperset e;
                if (!Enum.TryParse(value.ToString(), out e))
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Enum type exception.
        /// </summary>
        public class EnumTypeException<TEnum> : Exception
        {
            public static void VerifyEnum()
            {
                if (!typeof(TEnum).IsEnum)
                {
                    throw new EnumTypeException<TEnum>();
                }
            }

            /// <summary>
            /// Initializes a new instance of the <see cref="Shared.Helpers+EnumTypeException`1"/> class.
            /// </summary>
            public EnumTypeException()
                : base("Non enum type " + typeof(TEnum).Name + " used as enum")
            {}
        }

        /// <summary>
        /// Ensures the file path exists.
        /// </summary>
        /// <param name="file">File.</param>
        public static void EnsureFilePathExists(string file)
        {
            var dirname = Path.GetDirectoryName(file);
            if (!string.IsNullOrEmpty(dirname) && !Directory.Exists(dirname))
            {
                Directory.CreateDirectory(dirname);
            }
        }

        /// <summary>
        /// Gets the file data lines.
        /// </summary>
        /// <returns>The file data lines.</returns>
        /// <param name="file">File.</param>
        /// <param name="hasHeader">If set to <c>true</c> has header.</param>
        public static IEnumerable<string> GetFileDataLines(string file, bool hasHeader)
        {
            using (TextReader tr = new StreamReader(file))
            {
                if (hasHeader)
                {
                    tr.ReadLine();
                }

                return tr.ReadToEnd().Split('\n').Where(line => !string.IsNullOrWhiteSpace(line));
            }
        }

        /// <summary>
        /// Creates a series of 1-based numbers in sequence
        /// </summary>
        /// <param name="n">N.</param>
        public static List<int> Sequence(int n)
        {
            List<int> l = new List<int>(n);
            for (int i = 0; i < n; i++)
            {
                l.Add(i + 1);
            }

            return l;
        }

        /// <summary>
        /// Checks an initializes a field backing a property
        /// </summary>
        /// <returns>The init.</returns>
        /// <param name="item">Item.</param>
        /// <param name="initializer">Initializer.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
        public static T CheckInit<T>(ref T item, Func<T> initializer) where T : class
		{
			if (item == null)
			{
				item = initializer();
			}

			return item;
		}
    }

	/// <summary>
	/// Dictionary that reports missing requested values on miss
	/// </summary>
	public class CheckLookupDictionary<TKey, TValue>
	{
        /// <summary>
        /// Initializes a new instance of the <see cref="Shared.CheckLookupDictionary`2"/> class.
        /// </summary>
        /// <param name="d">D.</param>
		public CheckLookupDictionary(Dictionary<TKey, TValue> d)
		{
			this.d = d;
		}
		
        /// <summary>
        /// Gets the <see cref="Shared.CheckLookupDictionary`2"/> with the specified i.
        /// </summary>
        /// <param name="i">The index.</param>
		public virtual TValue this [TKey i]
		{
			get
			{
				if (d.ContainsKey(i))
				{
					return d[i];
				}
				else
				{
					throw new Exception(string.Format("Missing argument {0}", i));
				}
			}
		}

        /// <summary>
        /// Tries to get an optional parameter
        /// </summary>
        /// <returns><c>true</c>, if optional was tryed, <c>false</c> otherwise.</returns>
        /// <param name="i">The index.</param>
        /// <param name="v">V.</param>
		public virtual bool TryOptional(TKey i, out TValue v)
		{
			if (d.ContainsKey(i))
			{
				v = d[i];
				return true;
			}
			v = default(TValue);
			return false;
		}
		
        /// <summary>
        /// Gets an optional parameter.
        /// </summary>
        /// <returns>The optional.</returns>
        /// <param name="i">The index.</param>
		public virtual TValue GetOptional(TKey i)
		{
			if (d.ContainsKey(i))
			{
				return d[i];
			}
			return default(TValue);
		}

        /// <summary>
        /// Checks if the key is in the dictionary
        /// </summary>
        /// <returns><c>true</c>, if key was containsed, <c>false</c> otherwise.</returns>
        /// <param name="i">The index.</param>
        public virtual bool ContainsKey(TKey i)
        {
            return d.ContainsKey(i);
        }

        public virtual void Add(TKey k, TValue v)
        {
            d.Add(k, v);
        }

        public virtual void Remove(TKey k)
        {
            d.Remove(k);
        }

        /// <summary>
        /// Gets the count.
        /// </summary>
        /// <value>The count.</value>
        public int Count
        {
            get
            {
                return d.Count;
            }
        }
		
        /// <summary>
        /// Internal dictionary
        /// </summary>
		public Dictionary<TKey, TValue> d;
	}

    /// <summary>
    /// Enum dictionary.
    /// </summary>
    public class EnumDictionary<TValue> : CheckLookupDictionary<string, TValue>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Shared.EnumDictionary`1"/> class.
        /// </summary>
        /// <param name="d">D.</param>
        public EnumDictionary(Dictionary<string, TValue> d)
            : base(d)
        {
        }

        /// <summary>
        /// Gets the <see cref="Shared.EnumDictionary`1"/> with the specified i.
        /// </summary>
        /// <param name="i">The index.</param>
        public TValue this [System.Enum i]
        {
            get
            {
                return base[i.ToString()];
            }
        }

        /// <summary>
        /// Gets an optional parameter.
        /// </summary>
        /// <returns>The optional.</returns>
        /// <param name="i">The index.</param>
        public TValue GetOptional(System.Enum i)
        {
            if (d.ContainsKey(i.ToString()))
            {
                return d[i.ToString()];
            }
            return default(TValue);
        }

        /// <summary>
        /// Checks if the key is in the dictionary
        /// </summary>
        /// <returns><c>true</c>, if key was containsed, <c>false</c> otherwise.</returns>
        /// <param name="i">The index.</param>
        public virtual bool ContainsKey(System.Enum i)
        {
            return d.ContainsKey(i.ToString());
        }

        public virtual void Add(System.Enum i, TValue v)
        {
            d.Add(i.ToString(), v);
        }


        public virtual void Remove(System.Enum i)
        {
            d.Remove(i.ToString());
        }
    }

	/// <summary>
	/// Dictionary for getting values from XML data
	/// </summary>
	public class XmlDictionary<TValue>
	{
		/// <summary>
		/// Core xml document
		/// </summary>
		private XDocument document;

		/// <summary>
		/// String parser/converter
		/// </summary>
		private Func<string, TValue> typeConverter;

		/// <summary>
		/// The default value for missing values
		/// </summary>
		private TValue defaultValue;

		/// <summary>
        /// Initializes a new instance of the <see cref="Shared.XmlDictionary`1"/> class.
		/// </summary>
		/// <param name="document">Document.</param>
		/// <param name="typeConverter">Type converter.</param>
		/// <param name="defaultValue">Default value.</param>
		public XmlDictionary(XDocument document, Func<string, TValue> typeConverter, TValue defaultValue)
		{
			this.document = document;
			this.typeConverter = typeConverter;
			this.defaultValue = defaultValue;
		}

		/// <summary>
		/// Gets the <see cref="Base.XmlDictionary`1"/> with the specified i. Throws with informative exception on failed lookup
		/// </summary>
		/// <param name="i">The index.</param>
		public virtual TValue this [string i]
		{
			get
			{
				TValue v;
				if (TryOptional(i, out v))
				{
					return v;
				}
				else
				{
					throw new Exception(string.Format("Missing xml tag: {0}", i));
				}
			}
		}

        /// <summary>
        /// Gets the <see cref="Shared.XmlDictionary`1"/> with the specified e.
        /// </summary>
        /// <param name="e">E.</param>
        public virtual TValue this [Enum e]
        {
            get
            {
                return this[e.ToString()];
            }
        }

		/// <summary>
		/// Gets value at i if present, default if not; reports on missed lookup
		/// </summary>
		/// <returns><c>true</c>, if optional was tryed, <c>false</c> otherwise.</returns>
		/// <param name="i">The index.</param>
		/// <param name="v">V.</param>
		public virtual bool TryOptional(string i, out TValue v)
		{
			var elem = Lookup(i);

			if (elem == null)
			{
				v = this.defaultValue;
				return false;
			}
			else
			{
				v = this.typeConverter(elem.Value);
				return true;
			}
		}

		/// <summary>
		/// Containses the key.
		/// </summary>
		/// <returns><c>true</c>, if key was containsed, <c>false</c> otherwise.</returns>
		/// <param name="i">The index.</param>
		public virtual bool ContainsKey(string i)
		{
			return this.document.Descendants(i).Any();
		}

		/// <summary>
		/// Gets value at i if present, default if not
		/// </summary>
		/// <returns>The optional.</returns>
		/// <param name="i">The index.</param>
		public virtual TValue GetOptional(string i)
		{
			TValue v;
			this.TryOptional(i, out v);
			return v;
		}

		/// <summary>
		/// Lookup the specified path.
		/// </summary>
		/// <param name="path">Path.</param>
		private XElement Lookup(string path)
		{
			XElement root = document.Root;
			IEnumerable<XElement> descendants = null;
			foreach (string tag in path.Split('/'))
			{
				descendants = root.Descendants(tag);
				if (!descendants.Any())
				{
					return null;
				}

				root = descendants.First();
			}

			return root;
		}
	}
}

