cismapperlibexecdir = ${libexecdir}/${PACKAGE_NAME}-${PACKAGE_VERSION}

if CISMAPPER_AVAILABLE
cismapperlibexec_SCRIPTS = Analyses.dll Data.dll Genomics.dll Shared.dll Tools.dll
endif

Shared_lib_references = \
  -r:System \
	-r:System.Xml \
	-r:System.Xml.Linq
Shared_lib_source = \
	Shared/Properties/AssemblyInfo.cs \
	Shared/Args.cs \
	Shared/Base.cs \
	Shared/ConfigBase.cs \
	Shared/ConfigBase.cs \
	Shared/ITaskExecutor.cs \
	Shared/IUnknown.cs \
	Shared/RegressionConfig.cs

if CISMAPPER_AVAILABLE
Shared.dll: $(Shared_lib_source)
	$(MCS) $(Shared_lib_references) -target:library -out:$@ $^ 
endif

Tools_lib_references = \
  -r:System 
Tools_lib_source = \
	Tools/Properties/AssemblyInfo.cs \
	Tools/Stats.cs \
	Tools/Utilities.cs
Tools_extra = Tools/RegressionHelpers.R

if CISMAPPER_AVAILABLE
Tools.dll: $(Tools_lib_source)
	$(MCS) $(Tools_lib_references) -target:library -out:$@ $^ 
endif

Data_lib_references = \
	-r:Shared.dll \
  -r:System \
	-r:Tools.dll
Data_lib_source = \
	Data/Properties/AssemblyInfo.cs \
	Data/AxisDescription.cs \
	Data/Calculations.cs \
	Data/FactorDescription.cs \
	Data/Functions.cs \
	Data/OutputFiles.cs \
	Data/PlotDescription.cs \
	Data/Plots.cs \
	Data/R.cs \
	Data/Regression.cs \
	Data/Stats.cs \
	Data/Tables.cs

if CISMAPPER_AVAILABLE
Data.dll: Shared.dll Tools.dll $(Data_lib_source)
	$(MCS) $(Data_lib_references) -target:library -out:$@ $(Data_lib_source)
endif

Genomics_lib_references = \
	-r:Data.dll \
	-r:Shared.dll \
  -r:System \
	-r:Tools.dll
Genomics_lib_source = \
	Genomics/Properties/AssemblyInfo.cs \
	Genomics/Alphabets.cs \
	Genomics/BedExpressionFile.cs \
	Genomics/BedFile.cs \
	Genomics/CuffDiffFile.cs \
  Genomics/FastqEntry.cs \
  Genomics/FastqParser.cs \
  Genomics/Genomics.cs \
  Genomics/GtfExpressionFile.cs \
  Genomics/GtfGencodeFile.cs \
  Genomics/GtfLinkingFile.cs \
  Genomics/IAnnotation.cs \
  Genomics/IExpressionData.cs \
  Genomics/JunctionFile.cs \
  Genomics/KmerCollection.cs \
  Genomics/LocusRegulatoryMap.cs \
  Genomics/MapLink.cs \
  Genomics/MapLinkFilter.cs \
  Genomics/RegulatoryMap.cs \
  Genomics/Trie/Trie.cs \
  Genomics/Trie/TrieNode.cs \
  Genomics/TRFScorer.cs \
  Genomics/TssRegulatoryMap.cs

if CISMAPPER_AVAILABLE
Genomics.dll: Data.dll Shared.dll Tools.dll $(Genomics_lib_source)
	$(MCS) $(Genomics_lib_references) -target:library -out:$@ $(Genomics_lib_source)
endif

Analyses_lib_references = \
	-r:Data.dll \
	-r:Genomics.dll \
	-r:Shared.dll \
  -r:System 
Analyses_lib_source = \
	Analyses/Properties/AssemblyInfo.cs \
	Analyses/BaseAnalysis.cs \
	Analyses/BaseData.cs \
	Analyses/BaseMapData.cs \
	Analyses/BrowserTrackFromMap.cs \
	Analyses/ChIPSeqPeakBias.cs \
	Analyses/ConvertContactsToMap.cs \
	Analyses/ConvertMapToGenes.cs \
	Analyses/ConvertMapToScoredRegionPairs.cs \
	Analyses/CorrelationMapBuilder.cs \
	Analyses/CorrelationMapEligibleGenes.cs \
	Analyses/Elements/AnalysisElement.cs \
	Analyses/Elements/LocusFileElement.cs \
	Analyses/Elements/MapFileElement.cs \
	Analyses/IAnalysisExecutor.cs \
	Analyses/MapBuilder.cs \
	Analyses/MapFilter.cs \
	Analyses/MapStats.cs \
	Analyses/MergeCorrelationAndNullMaps.cs \
	Analyses/NullMapBuilder.cs \
	Analyses/PredictElements.cs \
	Analyses/PredictLinks.cs \
	Analyses/PredictTargets.cs \
	Analyses/PredictionBase.cs

if CISMAPPER_AVAILABLE
Analyses.dll: Data.dll Genomics.dll Shared.dll Tools.dll $(Analyses_lib_source)
	$(MCS) $(Analyses_lib_references) -target:library -out:$@ $(Analyses_lib_source)
endif

cismapper_source = \
	CisMapper/AssemblyInfo.cs \
	CisMapper/CisMapper.cs \
	CisMapper/MainClass.cs
cismapper_references = \
  -r:Analyses.dll \
	-r:Data.dll \
	-r:Genomics.dll \
	-r:Microsoft.CSharp \
	-r:Shared.dll \
	-r:System \
	-r:System.Core \
	-r:System.Xml \
	-r:System.Xml.Linq \
	-r:Tools.dll

if CISMAPPER_AVAILABLE
bin_SCRIPTS = cismapper.exe$(EXEEXT)
cismapper.exe$(EXEEXT): Analyses.dll Data.dll Genomics.dll Shared.dll Tools.dll $(cismapper_source)
	$(MCS) $(cismapper_references) -out:$@ $(cismapper_source)
endif

EXTRA_DIST = \
	$(Analyses_lib_source) \
	$(cismapper_source) \
	$(Data_lib_source) \
	$(Genomics_lib_source) \
	$(Shared_lib_source) \
	$(Tools_lib_source) \
	Tools/RegressionHelpers.R

CLEANFILES = \
	Analyses.dll \
	Data.dll \
	Genomics.dll \
	Shared.dll \
	Tools.dll

