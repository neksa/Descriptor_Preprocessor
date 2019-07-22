//--------------------------------------------------------------------------------
// <copyright file="MapLink.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using System.Linq;
    using Shared;
    using Tools;
    using System.Text.RegularExpressions;

    /// <summary>
    /// A single link in a regulatory mmap
    /// </summary>
    public class MapLink
    {
        private Location locusLocation;
        private Location transcriptLocation;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.MapLink"/> class.
        /// </summary>
        public MapLink()
        {
        }

        /// <summary>
        /// Gets or sets the confidence score.
        /// </summary>
        /// <value>The confidence score.</value>
        public double ConfidenceScore { get; set; }

        /// <summary>
        /// Gets or sets the correlation.
        /// </summary>
        /// <value>The correlation.</value>
        public double Correlation { get; set; }

        /// <summary>
        /// Gets or sets the length of the link.
        /// </summary>
        /// <value>The length of the link.</value>
        public int LinkLength { get; set; }

        /// <summary>
        /// Gets or sets the name of the transcript.
        /// </summary>
        /// <value>The name of the transcript.</value>
        public Tss TranscriptName { get; set; }

        /// <summary>
        /// Gets or sets the chromosome.
        /// </summary>
        /// <value>The chromosome.</value>
        public string Chromosome { get; set; }

        /// <summary>
        /// Gets or sets the tss position.
        /// </summary>
        /// <value>The tss position.</value>
        public int TssPosition { get; set; }

        /// <summary>
        /// Gets the tss location.
        /// </summary>
        /// <value>The tss location.</value>
        public Location TssLocation
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.transcriptLocation,
                    () => new Location
                    {
                        Chromosome = this.Chromosome,
                        Name = this.TranscriptName,
                        AlternateName = this.TssName,
                        Start = this.TssPosition,
                        End = this.TssPosition,
                        Strand = this.Strand,
                    });
            }
        }

        /// <summary>
        /// Gets or sets the size of the locus.
        /// </summary>
        /// <value>The size of the locus.</value>
        public int LocusSize 
        {
            get
            {
                return this.LocusEnd - this.LocusStart;
            }
        }

        /// <summary>
        /// Gets the locus start.
        /// </summary>
        /// <value>The locus start.</value>
        public int LocusStart
        {
            get
            {
                return this.LocusName.Start;
            }
        }

        /// <summary>
        /// Gets or sets the locus end.
        /// </summary>
        /// <value>The locus end.</value>
        public int LocusEnd 
        {
            get
            {
                return this.LocusName.End;
            }
        }

        /// <summary>
        /// Gets or sets the histone namw.
        /// </summary>
        /// <value>The histone namw.</value>
        public string HistoneName { get; set; }

        /// <summary>
        /// Gets or sets the name of the gene.
        /// </summary>
        /// <value>The name of the gene.</value>
        public string GeneName { get; set; }

        /// <summary>
        /// Gets or sets the name of the tss.
        /// </summary>
        /// <value>The name of the tss.</value>
        public string TssName { get; set; }

        /// <summary>
        /// Gets or sets the name of the locus.
        /// </summary>
        /// <value>The name of the locus.</value>
        public Locus LocusName { get; set; }

        /// <summary>
        /// Gets or sets the strand of the target gene/transcript/tss.
        /// </summary>
        /// <value>The strand.</value>
        public string Strand { get; set; }

        /// <summary>
        /// Gets the locus position.
        /// </summary>
        /// <value>The locus position.</value>
        public int LocusPosition
        {
            get
            {
                if (this.Strand == "+")
                {
                    return this.LinkLength;
                }
                else if (this.Strand == "-")
                {
                    return -this.LinkLength;
                }
                else
                {
                    throw new Exception("Missing strand in link data");
                }
            }
        }

        public Location LocusLocation
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.locusLocation,
                    () => new Location
                {
                    Chromosome = this.LocusName.Chr,
                    Name = this.LocusName,
                    Start = this.LocusStart,
                    End = this.LocusEnd,
                    Strand = this.Strand,
                });
            }
        }

        /// <summary>
        /// Gets the length of the abs link.
        /// </summary>
        /// <value>The length of the abs link.</value>
        public int AbsLinkLength
        {
            get
            {
                return Math.Abs(this.LinkLength);
            }
        }

        /// <summary>
        /// Serves as a hash function for a <see cref="Genomics.MapLink+LinkData"/> object.
        /// </summary>
        /// <returns>A hash code for this instance that is suitable for use in hashing algorithms and data structures such as
        /// a hash table.</returns>
        public override int GetHashCode()
        {
            return (this.TssName + "_" + this.LocusName).GetHashCode();
        }

        /// <summary>
        /// Determines whether the specified <see cref="System.Object"/> is equal to the current <see cref="Genomics.MapLink+LinkData"/>.
        /// </summary>
        /// <param name="obj">The <see cref="System.Object"/> to compare with the current <see cref="Genomics.MapLink+LinkData"/>.</param>
        /// <returns><c>true</c> if the specified <see cref="System.Object"/> is equal to the current
        /// <see cref="Genomics.MapLink+LinkData"/>; otherwise, <c>false</c>.</returns>
        public override bool Equals(object obj)
        {
            if (!(obj is MapLink))
            {
                return false;
            }

            var other = (MapLink)obj;
            return this.TssName == other.TssName && this.LocusName == other.LocusName;
        }

        public static bool operator ==(MapLink l1, MapLink l2)
        {
            return l1.Equals(l2);
        }

        public static bool operator !=(MapLink l1, MapLink l2)
        {
            return !l1.Equals(l2);
        }

        /// <summary>
        /// Test if two links have the same transcript and locus.
        /// </summary>
        public class SameTranscriptAndLocus : IEqualityComparer<MapLink>
        {
            /// <summary>
            /// Equals the specified a and b.
            /// </summary>
            /// <returns>true if the link has the same trancript and locus endpoints</returns>
            /// <param name="a">Link a.</param>
            /// <param name="b">Link b.</param>
            public bool Equals(MapLink a, MapLink b)
            {
                return a.Equals(b);
            }

            /// <Docs>The object for which the hash code is to be returned.</Docs>
            /// <para>Returns a hash code for the specified object.</para>
            /// <returns>A hash code for the specified object.</returns>
            /// <summary>
            /// Gets the hash code.
            /// </summary>
            /// <param name="a">The alpha component.</param>
            public int GetHashCode(MapLink a)
            {
                return a.GetHashCode();
            }
        }

        public class Locus : EqualityComparer<Locus>
        {
            private readonly string locusName;
            private string chr;
            private int start = -1;
            private int end = -1;

            private Regex nameFormat = new Regex("chr..?:[0-9]+-[0-9]+");

            public Locus(string s)
            {
                this.locusName = s;
            }

            public static implicit operator string(Locus s)
            {
                return s.locusName;
            }

            public static implicit operator Locus(string s)
            {
                return new Locus(s);
            }

            public override bool Equals(object obj)
            {
                return obj is Locus && ((Locus)obj).locusName == this.locusName;
            }

            public override int GetHashCode()
            {
                return this.locusName.GetHashCode();
            }

            public override string ToString()
            {
                return this.locusName;
            }

            public override bool Equals(Locus x, Locus y)
            {
                return x.Equals(y);
            }

            public override int GetHashCode(Locus obj)
            {
                return obj.GetHashCode();
            }

            public static bool operator ==(Locus c1, Locus c2)
            {
                return c1.locusName.Equals(c2.locusName);
            }

            public static bool operator !=(Locus c1, Locus c2)
            {
                return !c1.locusName.Equals(c2.locusName);
            }

            public string Chr
            {
                get
                {
                    return Helpers.CheckInit(
                        ref this.chr,
                        () => this.locusName.Split(':')[0]
                    );
                }
            }

            public int Start
            {
                get
                {
                    if (this.start == -1)
                    {
                        this.start = int.Parse(this.locusName.Split(':')[1].Split('-')[0]);
                    }

                    return this.start;
                }
            }

            public int End
            {
                get
                {
                    if (this.end == -1)
                    {
                        this.end = int.Parse(this.locusName.Split(':')[1].Split('-')[1]);
                    }

                    return this.end;
                }
            }
        }

        public class Tss : EqualityComparer<Tss>
        {
            private readonly string tssName;

            public Tss(string s)
            {
                this.tssName = s;
            }

            public static implicit operator string(Tss s)
            {
                return s.tssName;
            }

            public static implicit operator Tss(string s)
            {
                return new Tss(s);
            }

            public override bool Equals(object obj)
            {
                return obj is Tss && ((Tss)obj).tssName == this.tssName;
            }

            public override int GetHashCode()
            {
                return this.tssName.GetHashCode();
            }

            public override string ToString()
            {
                return this.tssName;
            }

            public override bool Equals(Tss x, Tss y)
            {
                return x.Equals(y);
            }

            public override int GetHashCode(Tss obj)
            {
                return obj.GetHashCode();
            }

            public static bool operator ==(Tss t1, Tss t2)
            {
                return t1.tssName.Equals(t2.tssName);
            }

            public static bool operator !=(Tss t1, Tss t2)
            {
                return !t1.tssName.Equals(t2.tssName);
            }
        }
    }
}
