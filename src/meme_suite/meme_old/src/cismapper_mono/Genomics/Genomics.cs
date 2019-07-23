//--------------------------------------------------------------------------------
// File: Genomics.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// A genomic locus
    /// </summary>
    public class Location 
    {
        /// <summary>
        /// Gets or sets location the name.
        /// </summary>
        /// <value>The name.</value>
        public string Name { get; set; }

        /// <summary>
        /// Chromosome of the location
        /// </summary>
        /// <value>The chromosome.</value>
        public string Chromosome { get; set; }

        /// <summary>
        /// Location start
        /// </summary>
        /// <value>The start.</value>
        public int Start { get; set; }

        /// <summary>
        /// Location end
        /// </summary>
        /// <value>The end.</value>
        public int End { get; set; }

        /// <summary>
        /// Location center
        /// </summary>
        /// <value>The middle.</value>
        public int Mid { get { return (this.Start + this.End) / 2; } }

        /// <summary>
        /// Gets the length.
        /// </summary>
        /// <value>The length.</value>
        public int Length { get { return (this.DirectionalEnd - this.DirectionalStart); } }

        /// <summary>
        /// Location's strand
        /// </summary>
        /// <value>The strand.</value>
        public string Strand { get; set; } 

        /// <summary>
        /// Gets or sets the score associated with the location
        /// </summary>
        /// <value>The score.</value>
        public double Score { get; set; } 

        /// <summary>
        /// Raw data backing the location
        /// </summary>
        /// <value>The data.</value>
        public string[] Data { get; set; }

        /// <summary>
        /// Allow attaching an array of integers to this location
        /// </summary>
        /// <value>The int data.</value>
        public int[] IntData{ get; set; }

        /// <summary>
        /// Arbitrary named fields associated with the location
        /// </summary>
        /// <value>The additional fields.</value>
        public Dictionary<string, string> AdditionalFields { get; set; }

        /// <summary>
        /// Alternate name for this location
        /// </summary>
        /// <value>The name of the alternate.</value>
        public string AlternateName { get; set; }

        /// <summary>
        /// Gets the directional start.
        /// </summary>
        /// <value>The directional start.</value>
        public int DirectionalStart 
        {
            get
            {
                if (Strand == "-")
                {
                    return End - 1;
                }
                return Start;
            }
        }

        /// <summary>
        /// Gets the directional end.
        /// </summary>
        /// <value>The directional end.</value>
        public int DirectionalEnd 
        {
            get
            {
                if (Strand == "-")
                {
                    return Start;
                }
                return End - 1;
            }
        }

        /// <summary>
        /// Determines overlap with another location
        /// </summary>
        /// <param name="l">L.</param>
        public bool Overlaps(Location l)
        {
            //Console.WriteLine("{0} - {1}, {2} - {3}", this.Start, this.End, l.Start, l.End );
            return 
                this.Chromosome == l.Chromosome &&
                (Between(this.Start, l.Start,    l.End) ||
                 Between(this.End,   l.Start,    l.End) ||
                 Between(l.Start,    this.Start, this.End) ||
                 Between(l.End,      this.Start, this.End));
        }

        public static bool Overlaps(Location l1, Location l2)
        {
            return l1.Overlaps(l2);
        }

        /// <summary>
        /// Determines overlap with the directional start another location
        /// </summary>
        /// <param name="l">L.</param>
        public bool OverlapsDirectionalStart(Location l)
        {
            //Console.WriteLine("{0} - {1}, {2} - {3}", this.Start, this.End, l.Start, l.End );
            return 
                this.Chromosome == l.Chromosome &&
                (Between(this.Start, l.DirectionalStart,    l.DirectionalStart) ||
                    Between(this.End,   l.DirectionalStart,    l.DirectionalStart) ||
                    Between(l.DirectionalStart,    this.Start, this.End));
        }

        public static bool OverlapsDirectionalStart(Location l1, Location l2)
        {
            return l1.OverlapsDirectionalStart(l2);
        }

        /// <summary>
        /// Is the specified x between [start and end).
        /// </summary>
        /// <param name="x">The x coordinate.</param>
        /// <param name="start">Start.</param>
        /// <param name="end">End.</param>
        private static bool Between(int x, int start, int end)
        {
            //if (x >= start && x < end)
            //{
            //    Console.WriteLine("BETWEEN: {0}, {1} - {2}", x, start, end);
            //}

            return x >= start && x < end;
        }
    }
}

