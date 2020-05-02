// <copyright file="MapFileElement.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;
    using Genomics;
    using Shared;

    public class MapFileElement : AnalysisElement<TssRegulatoryMap>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.MapFileElement"/> class.
        /// </summary>
        /// <param name="propertyPrefix">Property prefix.</param>
        public MapFileElement(string propertyPrefix)
            : base(propertyPrefix)
        {
        }

        /// <summary>
        /// Gets or sets the name of the map file.
        /// </summary>
        /// <value>The name of the map file.</value>
        public string MapFileName { get; set; }

        /// <summary>
        /// Gets or sets the type of the map file link.
        /// </summary>
        /// <value>The type of the map file link.</value>
        public string MapFileLinkType { get; set; }

        /// <summary>
        /// Gets or sets the map file max range.
        /// </summary>
        /// <value>The map file max range.</value>
        public string MapFileMaxRange { get; set; }

        /// <summary>
        /// Gets or sets the map file promoter range.
        /// </summary>
        /// <value>The map file promoter range.</value>
        public string MapFilePromoterRange { get; set; }

        /// <summary>
        /// Gets or sets the map file confidence threshold.
        /// </summary>
        /// <value>The map file confidence threshold.</value>
        public string MapFileConfidenceThreshold { get; set; }

        /// <summary>
        /// Gets or sets the element.
        /// </summary>
        /// <value>The element.</value>
        public override TssRegulatoryMap Element
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.element,
                    () =>
                    {
                        var filter = new MapLinkFilter();

                        if (!string.IsNullOrEmpty(this.MapFileLinkType))
                        {
                            filter.LinkTypeFilter = (MapLinkFilter.LinkType)Enum.Parse(typeof(MapLinkFilter.LinkType), this.MapFileLinkType);
                        }

                        if (!string.IsNullOrEmpty(this.MapFileMaxRange))
                        {
                            filter.MaximumLinkLength = int.Parse(this.MapFileMaxRange);
                        }

                        if (!string.IsNullOrEmpty(this.MapFilePromoterRange))
                        {
                            var range = this.MapFilePromoterRange.Split(',');

                            filter.PromoterUpstreamRange = int.Parse(range[0]);
                            if (range.Length > 1)
                            {
                                filter.PromoterDownstreamRange = int.Parse(range[1]);
                            }
                            else
                            {
                                filter.PromoterDownstreamRange = filter.PromoterUpstreamRange;
                            }
                        }
                        
                        if (!string.IsNullOrEmpty(this.MapFileConfidenceThreshold))
                        {
                            filter.ConfidenceThreshold = double.Parse(this.MapFileConfidenceThreshold);
                        }

                        return TssRegulatoryMap.LoadMap(this.MapFileName, filter);
                    });
            }

        }

        /// <summary>
        /// Register the specified analysis.
        /// </summary>
        /// <param name="analysis">Analysis.</param>
        public override void Register(BaseAnalysis analysis)
        {
            this.RegisterProperty(analysis, "MapFileName", x => this.MapFileName = x);
            this.RegisterProperty(analysis, "MapFileLinkType", x => this.MapFileLinkType = x);
            this.RegisterProperty(analysis, "MapFileMaxRange", x => this.MapFileMaxRange = x);
            this.RegisterProperty(analysis, "MapFilePromoterRange", x => this.MapFilePromoterRange = x);
            this.RegisterProperty(analysis, "MapFileConfidenceThreshold", x => this.MapFileConfidenceThreshold = x);
        }
    }
}

