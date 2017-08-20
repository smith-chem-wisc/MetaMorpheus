using MzLibUtil;
using System.Collections.Generic;
using System;

namespace EngineLayer
{
    //Common Parameters
    public class CommonParameters
    {
        public CommonParameters()
        {
            MaxMissedCleavages = 2;
            MinPeptideLength = 5;
            MaxPeptideLength = null;
            //Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new AbsoluteTolerance(0.01);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            TotalPartitions = 1;
            LocalizeAll = true;
            
            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            //ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            ConserveMemory = true;
            MaxDegreeOfParallelism = 1;

            // Deconvolution stuff
            DoPrecursorDeconvolution = true;
            UseProvidedPrecursorInfo = true;
            DeconvolutionIntensityRatio = 4;
            DeconvolutionMaxAssumedChargeState = 10;
            DeconvolutionMassTolerance = new PpmTolerance(5);
            ScoreCutoff = 5;

            //GPTMD ListOfModsLocalize
            //ListOfModsLocalize = new List<Tuple<string, string>>();
            
        }


        //MetaMorpheusTask Parameters
        public int? MaxDegreeOfParallelism { get; set; }
        public bool LocalizeAll { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }
        //public MyTask TaskType { get; set; }

        public bool DoPrecursorDeconvolution { get; set; }
        public bool UseProvidedPrecursorInfo { get; set; }
        public double DeconvolutionIntensityRatio { get; set; }
        public int DeconvolutionMaxAssumedChargeState { get; set; }
        public Tolerance DeconvolutionMassTolerance { get; set; }

        //Common Parameters

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public int TotalPartitions { get; set; }

        public Protease Protease { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }

        public Tolerance ProductMassTolerance { get; set; }
        
        public bool ConserveMemory { get; set; }

        public double ScoreCutoff { get; set; }

        public static List<string> AllModLists { get; private set; }

    }

    //Search Parameters
    public class SearchParameters
    {
        public SearchParameters()
        {
            DisposeOfFileWhenDone = true;
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = true;
            DoQuantification = false;
            QuantifyPpmTol = 5;
            SearchTarget = true;
            SearchDecoy = true;
            DoHistogramAnalysis = false;
            DoLocalizationAnalysis = true;
            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;
            //MassDiffAcceptors = GlobalTaskLevelSettings.SearchModesKnown.Take(1).ToList();
        }
        public bool DisposeOfFileWhenDone { get; set; }
        public bool AddCompIons { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchTarget { get; set; }
        public bool SearchDecoy { get; set; }
        public List<MassDiffAcceptor> MassDiffAcceptors { get; set; }
        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }
        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }
        public SearchType SearchType { get; set; }
    }

    //Calibration Parameters
    public class CalibrationParameters
    {
        public CalibrationParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(10);
            NonLinearCalibration = true;
            WriteIntermediateFiles = false;
        }
        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }
    }

    //GPTMD Parameters
    public class GPTMDParameters
    {
        public GPTMDParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(2);
        //    ListOfModsGptmd = GlobalTaskLevelSettings.AllModsKnown.Where(b =>
        //        b.modificationType.Equals("Glycan") ||
        //        b.modificationType.Equals("Mod") ||
        //        b.modificationType.Equals("PeptideTermMod") ||
        //        b.modificationType.Equals("Metal") ||
        //        b.modificationType.Equals("ProteinTermMod")).Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();
        }
        public Tolerance PrecursorMassTolerance { get; set; }
        public List<Tuple<string, string>> ListOfModsGptmd { get; set; }
    }

    //XLSearch Parameters
    public class XLSearchParameters
    {
        public XLSearchParameters()
        {
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = false;
            DoQuantification = false;
            QuantifyPpmTol = 5;
            SearchDecoy = true;
            DoHistogramAnalysis = false;

            DoLocalizationAnalysis = true;
            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;
            CrosslinkerType = CrosslinkerType.DSS;
            CrosslinkSearchTopNum = 50;
            CrosslinkSearchWithAllBeta = false;
            UdXLkerName = null;
            UdXLkerCleavable = false;
            UdXLkerShortMass = null;
            UdXLkerLongMass = null;
            UdXLkerTotalMass = null;
            UdXLkerResidue = 'K';
            XLprecusorMsTl = new PpmTolerance(10);
            XLBetaPrecusorMsTl = new PpmTolerance(10);
            DisposeOfFileWhenDone = true;
        }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchDecoy { get; set; }

        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }

        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }

        public CrosslinkerType CrosslinkerType { get; set; }
        public int CrosslinkSearchTopNum { get; set; }
        public bool CrosslinkSearchWithAllBeta { get; set; }
        public string UdXLkerName { get; set; }
        public bool UdXLkerCleavable { get; set; }
        public double? UdXLkerTotalMass { get; set; }
        public double? UdXLkerShortMass { get; set; }
        public double? UdXLkerLongMass { get; set; }
        public char UdXLkerResidue { get; set; }
        public Tolerance XLprecusorMsTl { get; set; }
        public Tolerance XLBetaPrecusorMsTl { get; set; }
        public bool DisposeOfFileWhenDone { get; set; }
    }

}
