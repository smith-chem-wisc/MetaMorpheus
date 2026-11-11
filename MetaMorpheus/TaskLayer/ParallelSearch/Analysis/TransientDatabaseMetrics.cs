#nullable enable
using CsvHelper.Configuration.Attributes;
using EngineLayer;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using TaskLayer.ParallelSearch.Util.Converters;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>
/// Aggregated result from all analyzers
/// Stores results as a dynamic dictionary that can be serialized to CSV
/// </summary>
public class TransientDatabaseMetrics : IEquatable<TransientDatabaseMetrics>
{
    /// <summary>
    /// Aggregated result from all analyzers
    /// Stores results as a dynamic dictionary that can be serialized to CSV
    /// </summary>
    public TransientDatabaseMetrics(string dbname) 
    {
        this.DatabaseName = dbname;
    }

    // Empty constructor for csv helper. 
    public TransientDatabaseMetrics() { }

    public string DatabaseName { get; set; }
    
    /// <summary>
    /// Serialized representation of all analysis results for CSV storage
    /// This uses a custom converter to flatten the dictionary
    /// </summary>
    [Ignore]
    public Dictionary<string, object> Results { get; set; } = new();
    
    /// <summary>
    /// List of errors encountered during analysis
    /// </summary>
    [Ignore]
    public List<string> Errors { get; } = [];

    #region Core Metrics 
    
    // Basic counts
    public int TotalProteins { get; set; }
    public int TransientProteinCount { get; set; }
    public int TransientPeptideCount { get; set; }
    
    // PSM metrics
    public int TargetPsmsAtQValueThreshold { get; set; }
    public int TargetPsmsFromTransientDb { get; set; }
    public int TargetPsmsFromTransientDbAtQValueThreshold { get; set; }
    
    // Peptide metrics
    public int TargetPeptidesAtQValueThreshold { get; set; }
    public int TargetPeptidesFromTransientDb { get; set; }
    public int TargetPeptidesFromTransientDbAtQValueThreshold { get; set; }
    
    // Protein group metrics (0 if parsimony not run)
    public int TargetProteinGroupsAtQValueThreshold { get; set; }
    public int TargetProteinGroupsFromTransientDb { get; set; }
    public int TargetProteinGroupsFromTransientDbAtQValueThreshold { get; set; }

    [Optional] public double MedianPsmsPerProteinGroup { get; set; }
    [Optional] public double MedianPeptidesPerProteinGroup { get; set; }
    [Optional] public double MedianUniquePeptidesPerProteinGroup { get; set; }

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] AllPsmsPerProteinGroup { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] AllPeptidesPerProteinGroup { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] AllUniquePeptidesPerProteinGroup { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] AllSequenceCoverageFractions { get; set; } = Array.Empty<double>();

    [Optional]
    public double MedianSequenceCoverageFraction { get; set; } = double.NaN;

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] AllFragmentSequenceCoverageFractions { get; set; } = Array.Empty<double>();

    #endregion

    #region Organism Specificity Metrics 

    public int PsmTargets { get; set; }
    public int PsmDecoys { get; set; }
    public int PsmBacterialTargets { get; set; }
    public int PsmBacterialDecoys { get; set; }
    public int PsmBacterialAmbiguous { get; set; }
    public int PsmBacterialUnambiguousTargets { get; set; }
    public int PsmBacterialUnambiguousDecoys { get; set; }
    
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmBacterialUnambiguousTargetScores { get; set; } = Array.Empty<double>();
    
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmBacterialUnambiguousDecoyScores { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmPrecursorMassErrors { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmBacterialTargetDeltaScores { get; set; } = Array.Empty<double>();

    public int PeptideTargets { get; set; }
    public int PeptideDecoys { get; set; }
    public int PeptideBacterialTargets { get; set; }
    public int PeptideBacterialDecoys { get; set; }
    public int PeptideBacterialAmbiguous { get; set; }
    public int PeptideBacterialUnambiguousTargets { get; set; }
    public int PeptideBacterialUnambiguousDecoys { get; set; }

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptideBacterialUnambiguousTargetScores { get; set; } = Array.Empty<double>();
    
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptideBacterialUnambiguousDecoyScores { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptidePrecursorMassErrors { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptideBacterialTargetDeltaScores { get; set; } = Array.Empty<double>();


    public int ProteinGroupTargets { get; set; }
    public int ProteinGroupDecoys { get; set; }
    public int ProteinGroupBacterialTargets { get; set; }
    public int ProteinGroupBacterialDecoys { get; set; }
    public int ProteinGroupBacterialUnambiguousTargets { get; set; }
    public int ProteinGroupBacterialUnambiguousDecoys { get; set; }

    #endregion

    #region Fragment Ions

    public double Psm_Bidirectional_MedianTargets { get; set; }
    public double Psm_ComplementaryCount_MedianTargets { get; set; }
    public double Psm_SequenceCoverageFraction_MedianTargets { get; set; }
    public double Psm_Bidirectional_MedianDecoys { get; set; }
    public double Psm_ComplementaryCount_MedianDecoys { get; set; }
    public double Psm_SequenceCoverageFraction_MedianDecoys { get; set; }

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_BidirectionalTargets { get; set; } = Array.Empty<double>();
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_ComplementaryCountTargets { get; set; } = Array.Empty<double>();
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_SequenceCoverageFractionTargets { get; set; } = Array.Empty<double>();
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_BidirectionalDecoys { get; set; } = Array.Empty<double>();
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_ComplementaryCountDecoys { get; set; } = Array.Empty<double>();
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_SequenceCoverageFractionDecoys { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    [Optional] public double[] Psm_FragmentPPMErrors { get; set; } = Array.Empty<double>();

    public double Peptide_Bidirectional_MedianTargets { get; set; }
    public double Peptide_ComplementaryCount_MedianTargets { get; set; }
    public double Peptide_SequenceCoverageFraction_MedianTargets { get; set; }
    public double Peptide_Bidirectional_MedianDecoys { get; set; }
    public double Peptide_ComplementaryCount_MedianDecoys { get; set; }
    public double Peptide_SequenceCoverageFraction_MedianDecoys { get; set; }

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_BidirectionalTargets { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_ComplementaryCountTargets { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_SequenceCoverageFractionTargets { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_BidirectionalDecoys { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_ComplementaryCountDecoys { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_SequenceCoverageFractionDecoys { get; set; } = Array.Empty<double>();

    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    [Optional] public double[] Peptide_FragmentPPMErrors { get; set; } = Array.Empty<double>();

    #endregion

    #region Retention Time

    [Optional] public double Psm_MeanAbsoluteRtError { get; set; }

    [Optional] [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Psm_AllRtErrors { get; set; } = Array.Empty<double>();

    [Optional] public double Peptide_MeanAbsoluteRtError { get; set; }

    [Optional] [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] Peptide_AllRtErrors { get; set; } = Array.Empty<double>();

    #endregion

    #region Deconvolution

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmPrecursorDeconScores { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmPrecursorEnvelopePeakCounts { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PsmPrecursorFractionalIntensities { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptidePrecursorDeconScores { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptidePrecursorEnvelopePeakCounts { get; set; } = Array.Empty<double>();

    [Optional]
    [TypeConverter(typeof(SemiColonDelimitedToDoubleArrayTypeConverter))]
    public double[] PeptidePrecursorFractionalIntensities { get; set; } = Array.Empty<double>();

    #endregion

    /// <summary>
    /// Populates the Results dictionary from the typed properties
    /// Called after CSV deserialization
    /// </summary>
    public void PopulateResultsFromProperties()
    {
        Results.Clear();
        
        // Core metrics
        Results[BasicMetricCollector.TotalProteins] = TotalProteins;
        Results[BasicMetricCollector.TransientProteinCount] = TransientProteinCount;
        Results[BasicMetricCollector.TransientPeptideCount] = TransientPeptideCount;
        Results[BasicMetricCollector.TargetPsmsAtQValueThreshold] = TargetPsmsAtQValueThreshold;
        Results[BasicMetricCollector.TargetPsmsFromTransientDb] = TargetPsmsFromTransientDb;
        Results[BasicMetricCollector.TargetPsmsFromTransientDbAtQValueThreshold] = TargetPsmsFromTransientDbAtQValueThreshold;
        Results[BasicMetricCollector.TargetPeptidesAtQValueThreshold] = TargetPeptidesAtQValueThreshold;
        Results[BasicMetricCollector.TargetPeptidesFromTransientDb] = TargetPeptidesFromTransientDb;
        Results[BasicMetricCollector.TargetPeptidesFromTransientDbAtQValueThreshold] = TargetPeptidesFromTransientDbAtQValueThreshold;
        Results[ProteinGroupCollector.TargetProteinGroupsAtQValueThreshold] = TargetProteinGroupsAtQValueThreshold;
        Results[ProteinGroupCollector.TargetProteinGroupsFromTransientDb] = TargetProteinGroupsFromTransientDb;
        Results[ProteinGroupCollector.TargetProteinGroupsFromTransientDbAtQValueThreshold] = TargetProteinGroupsFromTransientDbAtQValueThreshold;
        
        // Organism specificity
        Results[PsmPeptideSearchCollector.PsmTargets] = PsmTargets;
        Results[PsmPeptideSearchCollector.PsmDecoys] = PsmDecoys;
        Results[PsmPeptideSearchCollector.PsmBacterialTargets] = PsmBacterialTargets;
        Results[PsmPeptideSearchCollector.PsmBacterialDecoys] = PsmBacterialDecoys;
        Results[PsmPeptideSearchCollector.PsmBacterialAmbiguous] = PsmBacterialAmbiguous;
        Results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargets] = PsmBacterialUnambiguousTargets;
        Results[PsmPeptideSearchCollector.PsmBacterialUnambiguousDecoys] = PsmBacterialUnambiguousDecoys;
        Results[PsmPeptideSearchCollector.PsmBacterialUnambiguousTargetScores] = PsmBacterialUnambiguousTargetScores;
        Results[PsmPeptideSearchCollector.PsmBacterialUnambiguousDecoyScores] = PsmBacterialUnambiguousDecoyScores;
        Results[PsmPeptideSearchCollector.PsmBacteriaTargetlDeltaScores] = PsmBacterialTargetDeltaScores;

        Results[PsmPeptideSearchCollector.PeptideTargets] = PeptideTargets;
        Results[PsmPeptideSearchCollector.PeptideDecoys] = PeptideDecoys;
        Results[PsmPeptideSearchCollector.PeptideBacterialTargets] = PeptideBacterialTargets;
        Results[PsmPeptideSearchCollector.PeptideBacterialDecoys] = PeptideBacterialDecoys;
        Results[PsmPeptideSearchCollector.PeptideBacterialAmbiguous] = PeptideBacterialAmbiguous;
        Results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargets] = PeptideBacterialUnambiguousTargets;
        Results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousDecoys] = PeptideBacterialUnambiguousDecoys;
        Results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargetScores] = PeptideBacterialUnambiguousTargetScores;
        Results[PsmPeptideSearchCollector.PeptideBacterialUnambiguousDecoyScores] = PeptideBacterialUnambiguousDecoyScores;
        Results[PsmPeptideSearchCollector.PeptideBacterialTargetDeltaScores] = PeptideBacterialTargetDeltaScores;
        
        Results[ProteinGroupCollector.ProteinGroupTargets] = ProteinGroupTargets;
        Results[ProteinGroupCollector.ProteinGroupDecoys] = ProteinGroupDecoys;
        Results[ProteinGroupCollector.ProteinGroupBacterialTargets] = ProteinGroupBacterialTargets;
        Results[ProteinGroupCollector.ProteinGroupBacterialDecoys] = ProteinGroupBacterialDecoys;
        Results[ProteinGroupCollector.ProteinGroupBacterialUnambiguousTargets] = ProteinGroupBacterialUnambiguousTargets;
        Results[ProteinGroupCollector.ProteinGroupBacterialUnambiguousDecoys] = ProteinGroupBacterialUnambiguousDecoys;

        Results[ProteinGroupCollector.MedianPeptidesPerProteinGroup] = MedianPeptidesPerProteinGroup;
        Results[ProteinGroupCollector.MedianUniquePeptidesPerProteinGroup] = MedianUniquePeptidesPerProteinGroup;
        Results[ProteinGroupCollector.MedianPsmsPerProteinGroup] = MedianPsmsPerProteinGroup;
        Results[ProteinGroupCollector.AllPeptidesPerProteinGroup] = AllPeptidesPerProteinGroup;
        Results[ProteinGroupCollector.AllUniquePeptidesPerProteinGroup] = AllUniquePeptidesPerProteinGroup;
        Results[ProteinGroupCollector.AllPsmsPerProteinGroup] = AllPsmsPerProteinGroup;
        Results[ProteinGroupCollector.AllSequenceCoverageFractions] = AllSequenceCoverageFractions;
        Results["MedianSequenceCoverageFraction"] = MedianSequenceCoverageFraction;
        Results[ProteinGroupCollector.AllFragmentSequenceCoverageFractions] = AllFragmentSequenceCoverageFractions;

        // Decon
        Results[PsmPeptideSearchCollector.PsmPrecusorDeconScores] = PsmPrecursorDeconScores;
        Results[PsmPeptideSearchCollector.PsmPrecursorMassErrors] = PsmPrecursorMassErrors;
        Results[PsmPeptideSearchCollector.PsmPrecursorEnvelopePeakCount] = PsmPrecursorEnvelopePeakCounts;
        Results[PsmPeptideSearchCollector.PsmPrecusorFractionalIntensity] = PsmPrecursorFractionalIntensities;
        Results[PsmPeptideSearchCollector.PeptidePrecusorDeconScores] = PeptidePrecursorDeconScores;
        Results[PsmPeptideSearchCollector.PeptidePrecursorMassErrors] = PeptidePrecursorMassErrors;
        Results[PsmPeptideSearchCollector.PeptidePrecursorEnvelopePeakCount] = PeptidePrecursorEnvelopePeakCounts;
        Results[PsmPeptideSearchCollector.PeptidePrecusorFractionalIntensity] = PeptidePrecursorFractionalIntensities;

        // Fragment Ion metrics - PSM
        Results[FragmentIonCollector.PSM_LongestIonSeriesBidirectionalTargets] = Psm_Bidirectional_MedianTargets;
        Results[FragmentIonCollector.PSM_ComplementaryIonCountTargets] = Psm_ComplementaryCount_MedianTargets;
        Results[FragmentIonCollector.PSM_SequenceCoverageFractionTargets] = Psm_SequenceCoverageFraction_MedianTargets;
        Results[FragmentIonCollector.PSM_LongestIonSeriesBidirectionalDecoys] = Psm_Bidirectional_MedianDecoys;
        Results[FragmentIonCollector.PSM_ComplementaryIonCountDecoys] = Psm_ComplementaryCount_MedianDecoys;
        Results[FragmentIonCollector.PSM_SequenceCoverageFractionDecoys] = Psm_SequenceCoverageFraction_MedianDecoys;
        Results[FragmentIonCollector.PSM_LongestIonSeriesBidirectional_AllTargets] = Psm_BidirectionalTargets;
        Results[FragmentIonCollector.PSM_ComplementaryIonCount_AllTargets] = Psm_ComplementaryCountTargets;
        Results[FragmentIonCollector.PSM_SequenceCoverageFraction_AllTargets] = Psm_SequenceCoverageFractionTargets;
        Results[FragmentIonCollector.PSM_LongestIonSeriesBidirectional_AllDecoys] = Psm_BidirectionalDecoys;
        Results[FragmentIonCollector.PSM_ComplementaryIonCount_AllDecoys] = Psm_ComplementaryCountDecoys;
        Results[FragmentIonCollector.PSM_SequenceCoverageFraction_AllDecoys] = Psm_SequenceCoverageFractionDecoys;
        Results[FragmentIonCollector.PSM_FragmentPPMErrors] = Psm_FragmentPPMErrors;

        // Fragment Ion metrics - Peptide
        Results[FragmentIonCollector.Peptide_LongestIonSeriesBidirectionalTargets] = Peptide_Bidirectional_MedianTargets;
        Results[FragmentIonCollector.Peptide_ComplementaryIonCountTargets] = Peptide_ComplementaryCount_MedianTargets;
        Results[FragmentIonCollector.Peptide_SequenceCoverageFractionTargets] = Peptide_SequenceCoverageFraction_MedianTargets;
        Results[FragmentIonCollector.Peptide_LongestIonSeriesBidirectionalDecoys] = Peptide_Bidirectional_MedianDecoys;
        Results[FragmentIonCollector.Peptide_ComplementaryIonCountDecoys] = Peptide_ComplementaryCount_MedianDecoys;
        Results[FragmentIonCollector.Peptide_SequenceCoverageFractionDecoys] = Peptide_SequenceCoverageFraction_MedianDecoys;
        Results[FragmentIonCollector.Peptide_LongestIonSeriesBidirectional_AllTargets] = Peptide_BidirectionalTargets;
        Results[FragmentIonCollector.Peptide_ComplementaryIonCount_AllTargets] = Peptide_ComplementaryCountTargets;
        Results[FragmentIonCollector.Peptide_SequenceCoverageFraction_AllTargets] = Peptide_SequenceCoverageFractionTargets;
        Results[FragmentIonCollector.Peptide_LongestIonSeriesBidirectional_AllDecoys] = Peptide_BidirectionalDecoys;
        Results[FragmentIonCollector.Peptide_ComplementaryIonCount_AllDecoys] = Peptide_ComplementaryCountDecoys;
        Results[FragmentIonCollector.Peptide_SequenceCoverageFraction_AllDecoys] = Peptide_SequenceCoverageFractionDecoys;
        Results[FragmentIonCollector.Peptide_FragmentPPMErrors] = Peptide_FragmentPPMErrors;

        // Retention Time metrics
        Results[RetentionTimeCollector.PsmMeanAbsoluteRtError] = Psm_MeanAbsoluteRtError;
        Results[RetentionTimeCollector.PsmAllRtErrors] = Psm_AllRtErrors;
        Results[RetentionTimeCollector.PeptideMeanAbsoluteRtError] = Peptide_MeanAbsoluteRtError;
        Results[RetentionTimeCollector.PeptideAllRtErrors] = Peptide_AllRtErrors;
    }

    /// <summary>
    /// Populates typed properties from the Results dictionary
    /// Called before CSV serialization
    /// </summary>
    public void PopulatePropertiesFromResults()
    {
        // Core metrics
        TotalProteins = GetValue<int>(BasicMetricCollector.TotalProteins);
        TransientProteinCount = GetValue<int>(BasicMetricCollector.TransientProteinCount);
        TransientPeptideCount = GetValue<int>(BasicMetricCollector.TransientPeptideCount);
        TargetPsmsAtQValueThreshold = GetValue<int>(BasicMetricCollector.TargetPsmsAtQValueThreshold);
        TargetPsmsFromTransientDb = GetValue<int>(BasicMetricCollector.TargetPsmsFromTransientDb);
        TargetPsmsFromTransientDbAtQValueThreshold = GetValue<int>(BasicMetricCollector.TargetPsmsFromTransientDbAtQValueThreshold);
        TargetPeptidesAtQValueThreshold = GetValue<int>(BasicMetricCollector.TargetPeptidesAtQValueThreshold);
        TargetPeptidesFromTransientDb = GetValue<int>(BasicMetricCollector.TargetPeptidesFromTransientDb);
        TargetPeptidesFromTransientDbAtQValueThreshold = GetValue<int>(BasicMetricCollector.TargetPeptidesFromTransientDbAtQValueThreshold);
        TargetProteinGroupsAtQValueThreshold = GetValue<int>(ProteinGroupCollector.TargetProteinGroupsAtQValueThreshold);
        TargetProteinGroupsFromTransientDb = GetValue<int>(ProteinGroupCollector.TargetProteinGroupsFromTransientDb);
        TargetProteinGroupsFromTransientDbAtQValueThreshold = GetValue<int>(ProteinGroupCollector.TargetProteinGroupsFromTransientDbAtQValueThreshold);
        
        // Organism specificity
        PsmTargets = GetValue<int>(PsmPeptideSearchCollector.PsmTargets);
        PsmDecoys = GetValue<int>(PsmPeptideSearchCollector.PsmDecoys);
        PsmBacterialTargets = GetValue<int>(PsmPeptideSearchCollector.PsmBacterialTargets);
        PsmBacterialDecoys = GetValue<int>(PsmPeptideSearchCollector.PsmBacterialDecoys);
        PsmBacterialAmbiguous = GetValue<int>(PsmPeptideSearchCollector.PsmBacterialAmbiguous);
        PsmBacterialUnambiguousTargets = GetValue<int>(PsmPeptideSearchCollector.PsmBacterialUnambiguousTargets);
        PsmBacterialUnambiguousDecoys = GetValue<int>(PsmPeptideSearchCollector.PsmBacterialUnambiguousDecoys);
        PsmBacterialUnambiguousTargetScores = GetValue<double[]>(PsmPeptideSearchCollector.PsmBacterialUnambiguousTargetScores) ?? Array.Empty<double>();
        PsmBacterialUnambiguousDecoyScores = GetValue<double[]>(PsmPeptideSearchCollector.PsmBacterialUnambiguousDecoyScores) ?? Array.Empty<double>();
        PsmBacterialTargetDeltaScores = GetValue<double[]>(PsmPeptideSearchCollector.PsmBacteriaTargetlDeltaScores) ?? Array.Empty<double>();
        
        PeptideTargets = GetValue<int>(PsmPeptideSearchCollector.PeptideTargets);
        PeptideDecoys = GetValue<int>(PsmPeptideSearchCollector.PeptideDecoys);
        PeptideBacterialTargets = GetValue<int>(PsmPeptideSearchCollector.PeptideBacterialTargets);
        PeptideBacterialDecoys = GetValue<int>(PsmPeptideSearchCollector.PeptideBacterialDecoys);
        PeptideBacterialAmbiguous = GetValue<int>(PsmPeptideSearchCollector.PeptideBacterialAmbiguous);
        PeptideBacterialUnambiguousTargets = GetValue<int>(PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargets);
        PeptideBacterialUnambiguousDecoys = GetValue<int>(PsmPeptideSearchCollector.PeptideBacterialUnambiguousDecoys);
        PeptideBacterialUnambiguousTargetScores = GetValue<double[]>(PsmPeptideSearchCollector.PeptideBacterialUnambiguousTargetScores) ?? Array.Empty<double>();
        PeptideBacterialUnambiguousDecoyScores = GetValue<double[]>(PsmPeptideSearchCollector.PeptideBacterialUnambiguousDecoyScores) ?? Array.Empty<double>();
        PeptideBacterialTargetDeltaScores = GetValue<double[]>(PsmPeptideSearchCollector.PeptideBacterialTargetDeltaScores) ?? Array.Empty<double>();

        // Decon
        PsmPrecursorDeconScores = GetValue<double[]>(PsmPeptideSearchCollector.PsmPrecusorDeconScores) ?? Array.Empty<double>();
        PsmPrecursorMassErrors = GetValue<double[]>(PsmPeptideSearchCollector.PsmPrecursorMassErrors) ?? Array.Empty<double>();
        PsmPrecursorEnvelopePeakCounts = GetValue<double[]>(PsmPeptideSearchCollector.PsmPrecursorEnvelopePeakCount) ?? Array.Empty<double>();
        PsmPrecursorFractionalIntensities = GetValue<double[]>(PsmPeptideSearchCollector.PsmPrecusorFractionalIntensity) ?? Array.Empty<double>();
        PeptidePrecursorDeconScores = GetValue<double[]>(PsmPeptideSearchCollector.PeptidePrecusorDeconScores) ?? Array.Empty<double>();
        PeptidePrecursorMassErrors = GetValue<double[]>(PsmPeptideSearchCollector.PeptidePrecursorMassErrors) ?? Array.Empty<double>();
        PeptidePrecursorEnvelopePeakCounts = GetValue<double[]>(PsmPeptideSearchCollector.PeptidePrecursorEnvelopePeakCount) ?? Array.Empty<double>();
        PeptidePrecursorFractionalIntensities = GetValue<double[]>(PsmPeptideSearchCollector.PeptidePrecusorFractionalIntensity) ?? Array.Empty<double>();

        ProteinGroupTargets = GetValue<int>(ProteinGroupCollector.ProteinGroupTargets);
        ProteinGroupDecoys = GetValue<int>(ProteinGroupCollector.ProteinGroupDecoys);
        ProteinGroupBacterialTargets = GetValue<int>(ProteinGroupCollector.ProteinGroupBacterialTargets);
        ProteinGroupBacterialDecoys = GetValue<int>(ProteinGroupCollector.ProteinGroupBacterialDecoys);
        ProteinGroupBacterialUnambiguousTargets = GetValue<int>(ProteinGroupCollector.ProteinGroupBacterialUnambiguousTargets);
        ProteinGroupBacterialUnambiguousDecoys = GetValue<int>(ProteinGroupCollector.ProteinGroupBacterialUnambiguousDecoys);
        MedianPeptidesPerProteinGroup = GetValue<double>(ProteinGroupCollector.MedianPeptidesPerProteinGroup);
        MedianUniquePeptidesPerProteinGroup = GetValue<double>(ProteinGroupCollector.MedianUniquePeptidesPerProteinGroup);
        MedianPsmsPerProteinGroup = GetValue<double>(ProteinGroupCollector.MedianPsmsPerProteinGroup);
        AllPeptidesPerProteinGroup = GetValue<double[]>(ProteinGroupCollector.AllPeptidesPerProteinGroup) ?? Array.Empty<double>();
        AllUniquePeptidesPerProteinGroup = GetValue<double[]>(ProteinGroupCollector.AllUniquePeptidesPerProteinGroup) ?? Array.Empty<double>();
        AllPsmsPerProteinGroup = GetValue<double[]>(ProteinGroupCollector.AllPsmsPerProteinGroup) ?? Array.Empty<double>();
        AllSequenceCoverageFractions = GetValue<double[]>(ProteinGroupCollector.AllSequenceCoverageFractions) ?? Array.Empty<double>();
        MedianSequenceCoverageFraction = GetValue<double>("MedianSequenceCoverageFraction", double.NaN);
        AllFragmentSequenceCoverageFractions = GetValue<double[]>(ProteinGroupCollector.AllFragmentSequenceCoverageFractions) ?? Array.Empty<double>();

        // Fragment Ion metrics - PSM
        Psm_Bidirectional_MedianTargets = GetValue<double>(FragmentIonCollector.PSM_LongestIonSeriesBidirectionalTargets);
        Psm_ComplementaryCount_MedianTargets = GetValue<double>(FragmentIonCollector.PSM_ComplementaryIonCountTargets);
        Psm_SequenceCoverageFraction_MedianTargets = GetValue<double>(FragmentIonCollector.PSM_SequenceCoverageFractionTargets);
        Psm_Bidirectional_MedianDecoys = GetValue<double>(FragmentIonCollector.PSM_LongestIonSeriesBidirectionalDecoys);
        Psm_ComplementaryCount_MedianDecoys = GetValue<double>(FragmentIonCollector.PSM_ComplementaryIonCountDecoys);
        Psm_SequenceCoverageFraction_MedianDecoys = GetValue<double>(FragmentIonCollector.PSM_SequenceCoverageFractionDecoys);
        Psm_BidirectionalTargets = GetValue<double[]>(FragmentIonCollector.PSM_LongestIonSeriesBidirectional_AllTargets) ?? Array.Empty<double>();
        Psm_ComplementaryCountTargets = GetValue<double[]>(FragmentIonCollector.PSM_ComplementaryIonCount_AllTargets) ?? Array.Empty<double>();
        Psm_SequenceCoverageFractionTargets = GetValue<double[]>(FragmentIonCollector.PSM_SequenceCoverageFraction_AllTargets) ?? Array.Empty<double>();
        Psm_BidirectionalDecoys = GetValue<double[]>(FragmentIonCollector.PSM_LongestIonSeriesBidirectional_AllDecoys) ?? Array.Empty<double>();
        Psm_ComplementaryCountDecoys = GetValue<double[]>(FragmentIonCollector.PSM_ComplementaryIonCount_AllDecoys) ?? Array.Empty<double>();
        Psm_SequenceCoverageFractionDecoys = GetValue<double[]>(FragmentIonCollector.PSM_SequenceCoverageFraction_AllDecoys) ?? Array.Empty<double>();
        Psm_FragmentPPMErrors = GetValue<double[]>(FragmentIonCollector.PSM_FragmentPPMErrors) ?? Array.Empty<double>();

        // Fragment Ion metrics - Peptide
        Peptide_Bidirectional_MedianTargets = GetValue<double>(FragmentIonCollector.Peptide_LongestIonSeriesBidirectionalTargets);
        Peptide_ComplementaryCount_MedianTargets = GetValue<double>(FragmentIonCollector.Peptide_ComplementaryIonCountTargets);
        Peptide_SequenceCoverageFraction_MedianTargets = GetValue<double>(FragmentIonCollector.Peptide_SequenceCoverageFractionTargets);
        Peptide_Bidirectional_MedianDecoys = GetValue<double>(FragmentIonCollector.Peptide_LongestIonSeriesBidirectionalDecoys);
        Peptide_ComplementaryCount_MedianDecoys = GetValue<double>(FragmentIonCollector.Peptide_ComplementaryIonCountDecoys);
        Peptide_SequenceCoverageFraction_MedianDecoys = GetValue<double>(FragmentIonCollector.Peptide_SequenceCoverageFractionDecoys);
        Peptide_BidirectionalTargets = GetValue<double[]>(FragmentIonCollector.Peptide_LongestIonSeriesBidirectional_AllTargets) ?? Array.Empty<double>();
        Peptide_ComplementaryCountTargets = GetValue<double[]>(FragmentIonCollector.Peptide_ComplementaryIonCount_AllTargets) ?? Array.Empty<double>();
        Peptide_SequenceCoverageFractionTargets = GetValue<double[]>(FragmentIonCollector.Peptide_SequenceCoverageFraction_AllTargets) ?? Array.Empty<double>();
        Peptide_BidirectionalDecoys = GetValue<double[]>(FragmentIonCollector.Peptide_LongestIonSeriesBidirectional_AllDecoys) ?? Array.Empty<double>();
        Peptide_ComplementaryCountDecoys = GetValue<double[]>(FragmentIonCollector.Peptide_ComplementaryIonCount_AllDecoys) ?? Array.Empty<double>();
        Peptide_SequenceCoverageFractionDecoys = GetValue<double[]>(FragmentIonCollector.Peptide_SequenceCoverageFraction_AllDecoys) ?? Array.Empty<double>();
        Peptide_FragmentPPMErrors = GetValue<double[]>(FragmentIonCollector.Peptide_FragmentPPMErrors) ?? Array.Empty<double>();

        // Retention Time metrics
        Psm_MeanAbsoluteRtError = GetValue<double>(RetentionTimeCollector.PsmMeanAbsoluteRtError);
        Psm_AllRtErrors = GetValue<double[]>(RetentionTimeCollector.PsmAllRtErrors) ?? Array.Empty<double>();
        Peptide_MeanAbsoluteRtError = GetValue<double>(RetentionTimeCollector.PeptideMeanAbsoluteRtError);
        Peptide_AllRtErrors = GetValue<double[]>(RetentionTimeCollector.PeptideAllRtErrors) ?? Array.Empty<double>();
    }

    /// <summary>
    /// Gets a typed value from the results dictionary with a default fallback
    /// </summary>
    public T? GetValue<T>(string key, T? defaultValue = default)
    {
        if (Results.TryGetValue(key, out var value))
        {
            try
            {
                if (value is T typedValue)
                    return typedValue;
                
                return (T)Convert.ChangeType(value, typeof(T));
            }
            catch
            {
                return defaultValue;
            }
        }
        return defaultValue;
    }

    /// <summary>
    /// Writes the database results to a text file
    /// </summary>
    public async Task WriteToTextFileAsync(string filePath, double qValueThreshold, bool doParsimony)
    {
        await using StreamWriter file = new StreamWriter(filePath);
        await file.WriteLineAsync($"Database: {DatabaseName}");
        await file.WriteLineAsync($"Total proteins in combined database: {TotalProteins}");
        await file.WriteLineAsync($"Total proteins from transient database: {TransientProteinCount}");
        await file.WriteLineAsync($"Total peptides from transient database: {TransientPeptideCount}");
        await file.WriteLineAsync();
        
        await file.WriteLineAsync($"Target PSMs at {qValueThreshold * 100}% FDR: {TargetPsmsAtQValueThreshold}");
        await file.WriteLineAsync($"Target PSMs from transient database: {TargetPsmsFromTransientDb}");
        await file.WriteLineAsync($"Target PSMs from transient database at {qValueThreshold * 100}% FDR: {TargetPsmsFromTransientDbAtQValueThreshold}");
        await file.WriteLineAsync($"PSM Bacterial Targets: {PsmBacterialTargets}");
        await file.WriteLineAsync($"PSM Bacterial Unambiguous Targets: {PsmBacterialUnambiguousTargets}");
        await file.WriteLineAsync();
        
        await file.WriteLineAsync($"Target peptides at {qValueThreshold * 100}% FDR: {TargetPeptidesAtQValueThreshold}");
        await file.WriteLineAsync($"Target peptides from transient database: {TargetPeptidesFromTransientDb}");
        await file.WriteLineAsync($"Target peptides from transient database at {qValueThreshold * 100}% FDR: {TargetPeptidesFromTransientDbAtQValueThreshold}");
        await file.WriteLineAsync($"Peptide Bacterial Targets: {PeptideBacterialTargets}");
        await file.WriteLineAsync($"Peptide Bacterial Unambiguous Targets: {PeptideBacterialUnambiguousTargets}");

        if (doParsimony)
        {
            await file.WriteLineAsync();
            await file.WriteLineAsync($"Target protein groups at {qValueThreshold * 100}% FDR: {TargetProteinGroupsAtQValueThreshold}");
            await file.WriteLineAsync($"Target protein groups with transient database proteins: {TargetProteinGroupsFromTransientDb}");
            await file.WriteLineAsync($"Target protein groups with transient database proteins at {qValueThreshold * 100}% FDR: {TargetProteinGroupsFromTransientDbAtQValueThreshold}");
            await file.WriteLineAsync($"Protein Group Bacterial Targets: {ProteinGroupBacterialTargets}");
            await file.WriteLineAsync($"Protein Group Bacterial Unambiguous Targets: {ProteinGroupBacterialUnambiguousTargets}");
        }
    }

    public bool Equals(TransientDatabaseMetrics? other)
    {
        if (other is null) return false;
        if (other.DatabaseName != DatabaseName) return false;
        if (other.Results.Count != Results.Count) return false;
        return true;
    }

    public override int GetHashCode()
    {
        return HashCode.Combine(DatabaseName, Results.Count);
    }

    public TransientDatabaseMetrics Add(TransientDatabaseMetrics other)
    {
        if (other == null) return this;
        if (this.DatabaseName != other.DatabaseName) throw new MetaMorpheusException("Cannot add metrics from different databases", new ArgumentException("Database names do not match"));

        var metrics = new TransientDatabaseMetrics(this.DatabaseName);
        foreach (var key in Results.Keys)
        {
            if (other.Results.ContainsKey(key))
            {
                var value1 = Results[key];
                var value2 = other.Results[key];
                if (value1 is int int1 && value2 is int int2)
                {
                    metrics.Results[key] = int1 + int2;
                }
                else if (value1 is double double1 && value2 is double double2)
                {
                    metrics.Results[key] = double1 + double2;
                }
                else if (value1 is double[] arr1 && value2 is double[] arr2)
                {
                    metrics.Results[key] = arr1.Concat(arr2).ToArray();
                }
                else
                {
                    throw new MetaMorpheusException("Must add new if case to support new type for addition", new ArgumentException($"Type {value1.GetType()} is not supported for addition"));
                }
            }
            else
            {
                metrics.Results[key] = Results[key];
            }
        }


        return metrics;
    }
}
