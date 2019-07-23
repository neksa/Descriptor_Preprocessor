//--------------------------------------------------------------------------------
// File: ConfigBase.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Shared
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Xml.Linq;

    /// <summary>
    /// Basic xml configuration file object representation.
    /// </summary>
	public class ConfigBase
	{
		// * Private Members ************************************************

		// * Constructors ************************************************

		/// <summary>
		/// Initialize new document of config data
		/// </summary>
		/// <param name="xmlFile">Xml file.</param>
		public ConfigBase(string xmlFile)
		{
			Document = XDocument.Load(xmlFile);

			StringValues = new XmlDictionary<string>(Document, s => s, string.Empty); 
			CsvValues = new XmlDictionary<List<string>>(Document, s => s.Split(',').ToList(), new List<string> {} ); 
            CsvIntValues = new XmlDictionary<List<int>>(Document, s => s.Split(',').Select(v => int.Parse(v)).ToList(), new List<int> {} );
			IntValues = new XmlDictionary<int>(Document, s => int.Parse(s), 0); 
			DoubleValues = new XmlDictionary<double>(Document, s => double.Parse(s), 0.0); 
			NameValuePairs = new XmlDictionary<Dictionary<string, string>>(
				Document,
				s => s.Split(';').Select(nvp => 
		            {
						var fields = nvp.Split(',');
						return new
						{
							Name = fields[0],
							Value = fields[1]
						};
					}).ToDictionary(x => x.Name, x => x.Value),
				new Dictionary<string, string>{});

			NameValueLists = new XmlDictionary<ILookup<string, string>>(
				Document,
				s => s.Split(';').Select(nvp => 
         			{
						var fields = nvp.Split(',');
						return new
						{
							Name = fields[0],
							Value = fields[1]
						};
					}).ToLookup(x => x.Name, x => x.Value),
				new Dictionary<string, string>{}.ToLookup(x => x.Key, x => x.Value));
		}

		// * Public Properties ************************************************

		/// <summary>
		/// String value document lookup
		/// </summary>
		public XmlDictionary<string> StringValues { get; set; }
		
		/// <summary>
		/// Csv string value document lookup
		/// </summary>
		public XmlDictionary<List<string>> CsvValues;

        /// <summary>
        /// CSV string value of an int list
        /// </summary>
        public XmlDictionary<List<int>> CsvIntValues;
		
		/// <summary>
		/// Int value document lookup
		/// </summary>
		public XmlDictionary<int> IntValues;
		
		/// <summary>
		/// Double value document lookup
		/// </summary>
		public XmlDictionary<double> DoubleValues;

		/// <summary>
		/// Name pair values in XML
		/// </summary>
		public XmlDictionary<Dictionary<string, string>> NameValuePairs;

		/// <summary>
		/// Names lists of values
		/// </summary>
		public XmlDictionary<ILookup<string, string>> NameValueLists;

        // * Protected Properties ************************************************

		/// <summary>
		/// Core document of config data
		/// </summary>
		/// <param name="xmlFile">Xml file.</param>
		protected XDocument Document { get; set; }

		// * Public Methods ************************************************

        /// <summary>
        /// Validate the specified argumentEnum.
        /// </summary>
        /// <param name="argumentEnum">Argument enum.</param>
        public void Validate(Type argumentEnum)
        {
            HashSet<string> names = new HashSet<string>(Enum.GetNames(argumentEnum));

            var missingConfiguration = names.Where(x => !this.NameValuePairs.ContainsKey(x)).ToList();

            if (missingConfiguration.Count > 0)
            {
                throw new Exception("Missing configuration members found: " + string.Join(", ", missingConfiguration));
            }
        }

		// * Protected Methods ************************************************
	}
}

