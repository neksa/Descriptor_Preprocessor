//--------------------------------------------------------------------------------
// File: BaseMapData.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.IO;
    using Genomics;
    using Shared;

    /// <summary>
    /// Base experiment data that includes a TSS-Locus map on top of the elements of BaseData.
    /// </summary>
	public class BaseMapData : BaseData
    {
        /// <summary>
        /// The name of the map file.
        /// </summary>
		protected string MapFileName;

        /// <summary>
        /// The Locus set.
        /// </summary>
		private HashSet<string> locusSet;

        /// <summary>
        /// The max range.
        /// </summary>
        private int maxRange = -1;

        /// <summary>
        /// The nearest neighbor map.
        /// </summary>
        //private Dictionary<string, Dictionary<string, MapLink>> nearestNeighborMap;

	/// <summary>
	/// The link data from the map
	/// </summary>
	private TssRegulatoryMap links;

        /// <summary>
        /// Initializes a new instance of the <see cref="Genomics.BaseMapData"/> class.
        /// </summary>
        /// <param name="xmlFile">Xml file.</param>
        /// <param name="mapFileName">Map file name.</param>
		public BaseMapData(string outputRoot, string plotRoot, string scriptRoot, string xmlFile, string mapFileName)
			: base(outputRoot, plotRoot, scriptRoot, xmlFile)
        {
			this.MapFileName = mapFileName;
        }

        /// <summary>
        /// Loads a map.
        /// </summary>
        /// <returns>The map.</returns>
        /// <param name="mapFileName">Map file name.</param>
        /// <param name="configData">Config data.</param>
        /// <param name="maxRange">Max range.</param>
        public static TssRegulatoryMap LoadMap(
            string mapFileName,
            List<int> promoterRange,
            int maxRange,
            double confidenceThreshold)
        {
            var mapProperties = new MapLinkFilter
            {
                PromoterUpstreamRange = promoterRange[0],
                PromoterDownstreamRange = promoterRange.Count == 1 ? promoterRange[0] : promoterRange[1],
                MaximumLinkLength = maxRange,
                ConfidenceThreshold = confidenceThreshold,
            };

            return TssRegulatoryMap.LoadMap(mapFileName, mapProperties);
        }

		/// <summary>
		/// Gets the map.
		/// </summary>
		/// <value>The map.</value>
        public virtual TssRegulatoryMap Map
		{
			get
			{
                return Helpers.CheckInit(
                    ref this.links,
                    () => LoadMap(this.MapFileName, this.ConfigData.PromoterRange, this.MaxRange, 1));
			}
		}

        /// <summary>
        /// Clears the map.
        /// </summary>
        public virtual void ClearMap()
        {
            this.links = null;
            GC.Collect();
        }

        /// <summary>
        /// Gets the max range.
        /// </summary>
        /// <value>The max range.</value>
        public int MaxRange 
        {
            get
            {
                if (this.maxRange == -1)
                {
                    this.maxRange = 1000000;

                    int tempRange = 0;
                    if (ConfigData.IntValues.TryOptional("MaxRange", out tempRange))
                    {
                        maxRange = tempRange;
                    }
                }

                return this.maxRange;
            }

            set
            {
                this.maxRange = value;
            }
        }

		

		/// <summary>
		/// Gets the Locus set.
		/// </summary>
		/// <value>The Locus set.</value>
		public HashSet<string> LocusSet
		{
			get
			{
				return Helpers.CheckInit(
					ref this.locusSet,
                    () => new HashSet<string>(this.Map.Select(x => x.Value.Keys.Select(y => (string)y)).SelectMany(x => x)));
			}
		}
    }
}

