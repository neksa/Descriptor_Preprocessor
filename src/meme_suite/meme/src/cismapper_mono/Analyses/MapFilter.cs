
namespace Analyses
{
    using System;
    using System.Linq;
    using System.Collections.Generic;
    using Genomics;
    using Shared;

    public class MapFilter 
    {
        public MapFilter(string mapFileName)
        {
            this.MapFileName = mapFileName;
        }

        public string MapFileName { get; set; }
        public string OutputFile { get; set; }
        public string HistoneName { get; set; }
        public IExpressionData ExpressionData { get; set; }

        public TssRegulatoryMap Map { get; set; }

        public void BestWorst(double fraction)
        {
            this.Map = TssRegulatoryMap.LoadMap(this.MapFileName, new MapLinkFilter { });

            var bestLinkMap = this.Map.GetBestNeighborMap(1);

            var orderedLinks = bestLinkMap.Links.OrderBy(x => x.ConfidenceScore).ToList();

            int linkCount = (int)(orderedLinks.Count * fraction);

            var topLinks = orderedLinks
                .Take(linkCount)
                .OrderBy(x => x.AbsLinkLength)
                .ToList();

            var bottomLinks = orderedLinks
                .OrderBy(x => -x.ConfidenceScore)
                .Take(orderedLinks.Count - linkCount)
                .ToList();

            List<MapLink> bestLinks = new List<MapLink>();
            List<MapLink> worstLinks = new List<MapLink>();

            foreach (var link in topLinks)
            {
                int linkIndex = int.MaxValue;
                for (int i = 0; i < bottomLinks.Count; i++)
                {
                    var bottomLink = bottomLinks[i];
                    if (Math.Sign(link.LinkLength) == Math.Sign(bottomLink.LinkLength) &&
                        link.AbsLinkLength < bottomLink.AbsLinkLength * 1.01 &&
                        link.AbsLinkLength > bottomLink.AbsLinkLength * 0.99)
                    {
                        linkIndex = i;
                        break;
                    }
                }

                if (linkIndex != int.MaxValue)
                {
                    bestLinks.Add(link);
                    worstLinks.Add(bottomLinks[linkIndex]);
                    bottomLinks.RemoveAt(linkIndex);
                }
            }

            var bottomMapFile = this.OutputFile.Replace(".bed", ".bottom.bed");

            NullMapBuilder.WriteMap(new TssRegulatoryMap(bestLinks), this.ExpressionData.Genes, this.HistoneName, this.OutputFile);
            NullMapBuilder.WriteMap(new TssRegulatoryMap(worstLinks), this.ExpressionData.Genes, this.HistoneName, bottomMapFile);
        }

        public void Link(MapLinkFilter.LinkType type )
        {
            this.Map = TssRegulatoryMap.LoadMap(this.MapFileName, new MapLinkFilter { LinkTypeFilter = type });
            NullMapBuilder.WriteMap(this.Map, this.ExpressionData.Genes, this.HistoneName, this.OutputFile);
        }

        public class Executor : IAnalysisExecutor<MapFilter, MapFilter.Executor.Arguments>
        {
            public enum Arguments
            {
                MapFileName,
                FilterType,
                FilterParams,
                OutputFile,
            }

            public enum FilterType
            {
                Link,
                Confidence,
                BestWorst,
            }

            protected override Dictionary<Arguments, string> OptionsData
            {
                get
                {
                    return new Dictionary<Arguments, string>
                    {
                        { Arguments.MapFileName, "Map file to filter" },
                        { Arguments.FilterType, "Operation to perform" },
                        { Arguments.FilterParams, "FilterType-specific params" },
                        { Arguments.OutputFile, "Output map file" },
                    };
                }
            }

            override public void Execute(Args commandArgs)
            {
                var filter = new MapFilter(commandArgs.StringEnumArgs[Arguments.MapFileName]);

                filter.OutputFile = commandArgs.StringEnumArgs[Arguments.OutputFile];

                switch ((FilterType)Enum.Parse(typeof(FilterType), commandArgs.StringEnumArgs[Arguments.FilterType]))
                {
                    case FilterType.BestWorst:
                        {
                            string[] p = commandArgs.StringEnumArgs[Arguments.FilterParams].Split(',');

                            filter.ExpressionData = IExpressionData.LoadExpressionData(
                                p[1], 
                                "LongPap",
                                "gtf",
                                null);

                            filter.HistoneName = p[2];

                            filter.BestWorst(double.Parse(p[0]));
                            break;
                        }

                    case FilterType.Link:
                        {
                            string[] p = commandArgs.StringEnumArgs[Arguments.FilterParams].Split(',');

                            filter.ExpressionData = IExpressionData.LoadExpressionData(
                                p[1], 
                                "LongPap",
                                "gtf",
                                null);

                            filter.HistoneName = "None";
                            filter.Link((MapLinkFilter.LinkType)Enum.Parse(typeof(MapLinkFilter.LinkType), p[0]));
                            break;
                        }

                }

            }

            public override string Description
            {
                get
                {
                    return "Map filter";
                }
            }
        }
    }
}

