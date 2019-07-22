//--------------------------------------------------------------------------------
// File: BedExpressionFile.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using Genomics;
    using Shared;

    /// <summary>
    /// Wrapper for a bed file containing expression data
    /// </summary>
	public class BedExpressionFile : BedFile 
    {
		/// <summary>
		/// Tests if a transcript name should be included in the expression set
		/// </summary>
		private Func<string, bool> isValidTranscriptName;

		/// <summary>
		/// Writes out a Locus bed file from a list of locations
		/// </summary>
		/// <param name="locations">Locations.</param>
		/// <param name="filename">Filename.</param>
        public static void ToFileBedExpression(List<Genomics.Location> locations, string filename)
		{
			using (TextWriter tw = new StreamWriter(filename))
			{
				//tw.WriteLine(string.Join("\t", locations.Select(x => string.Join("\t", new string[] { x.Chromosome, "MapBuilder", "transcript", x.Start.ToString(), x.End.ToString(), "0", x.Strand, x.Name, x.Score } ))));
			}
		}

		/// <summary>
		/// Initializes a new instance of the <see cref="Genomics.BedExpressionFile"/> class.
		/// </summary>
		/// <param name="filename">Filename.</param>
		/// <param name="layout">Layout.</param>
		/// <param name="annotation">Annotation.</param>
		public BedExpressionFile(string filename, Layout layout, IAnnotation annotation)
			: base(filename, layout)
		{
			this.Annotation = annotation;
		}

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.BedExpressionFile"/> class.
        /// </summary>
        /// <param name="filename">Filename.</param>
        /// <param name="layout">Layout.</param>
		public BedExpressionFile(string filename, Layout layout)
			: base(filename, layout)
        {
        }

		/// <summary>
		/// Queries the interface.
		/// </summary>
		/// <returns>The interface.</returns>
		/// <param name="t">T.</param>
		override protected object QueryInterface(Type t)
		{
			if (t == typeof(IExpressionData))
			{
				return new ExpressionDataProxy(this) as IExpressionData;
			}

			return base.QueryInterface(t);
		}

		/// <summary>
		/// Gets a value indicating whether this instance is valid transcript name.
		/// </summary>
		/// <value><c>true</c> if this instance is valid transcript name; otherwise, <c>false</c>.</value>
		protected Func<string, bool> IsValidTranscriptName
		{
			get
			{
				return Helpers.CheckInit(
					ref this.isValidTranscriptName,
					() =>
					{
						if (this.Annotation == null)
						{
							return transcriptName => true;
						}
						else
						{
							return transcriptName => this.Annotation.ValidTranscripts.Contains(transcriptName);
						}
					});
			}
		}

		/// <summary>
		/// Gets or sets the annotation.
		/// </summary>
		/// <value>The annotation.</value>
		protected IAnnotation Annotation { get; set; }

		/// <summary>
		/// Parses the fields.
		/// </summary>
		/// <param name="fields">Fields.</param>
		/// <param name="layout">Layout.</param>
		/// <param name="data">Data.</param>
		/// <param name="entryCount">Entry count.</param>
        override protected void ParseFields(string[] fields, Layout layout, List<Tuple<Genomics.Location, string>> data, ref int entryCount)
		{
			if (this.IsValidTranscriptName(fields[layout.Name]))
			{
				base.ParseFields(fields, layout, data, ref entryCount);
			}
		}

        /// <summary>
        /// Gets the transcripts.
        /// </summary>
        /// <value>The transcripts.</value>
	    public Dictionary<string, Genomics.Location> Transcripts
        {
            get
            {
				return this.Locations;
            }
        }    

        private class ExpressionDataProxy : IExpressionData
		{
			private readonly BedExpressionFile data;

			public ExpressionDataProxy(BedExpressionFile data)
			{
				this.data = data;
			}

            override public Dictionary<string, Location> Transcripts { get { return this.data.Transcripts; } }
		}
    }


}

