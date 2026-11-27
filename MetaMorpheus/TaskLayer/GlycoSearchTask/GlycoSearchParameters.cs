using EngineLayer.GlycoSearch;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class GlycoSearchParameters
    {
        public GlycoSearchParameters()
        {
            OGlycanDatabasefile = "OGlycan.gdb";
            NGlycanDatabasefile = "NGlycan.gdb";
            SelectedGlycans = new List<(string, string)>();
            GlycoSearchType = GlycoSearchType.OGlycanSearch;
            OxoniumIonFilt = true;
            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 50;
            MaximumOGlycanAllowed = 4;
            MaximumGlycanBoxMass = EngineLayer.GlycanBox.DefaultMaximumGlycanBoxMass;
            MaximumSpectraFilesInParallel = 0;
            DoParsimony = true;
            NoOneHitWonders = false;
            ModPeptidesAreDifferent = false;

            //quantification options
            DoQuantification = false;
            DoMbrAnalysis = true;
            QuantifyPpmTol = 5;
            Normalize = false;

            //output options
            WriteIndividualFiles = false;
            WriteDecoys = true;
            WriteContaminants = true;
            WriteSpectrumLibrary = false;
            DisposeOfFileWhenDone = true;
            WritePrunedDataBase = false;

            ModsToWriteSelection = SearchParameters.DefaultModsToWriteSelection();
        }
        public string OGlycanDatabasefile { get; set; }
        public string NGlycanDatabasefile { get; set; }

        /// <summary>
        /// The individual glycans the user checked, as (database file name, glycan IdWithMotif) pairs.
        /// EMPTY MEANS THE WHOLE DATABASE -- which is what every task written before this existed says, so
        /// old TOMLs and users who never open the tree keep today's behaviour exactly.
        /// </summary>
        /// <remarks>
        /// Stored by composition string, never by Glycan.GlyId: GlyId is a positional index into the loaded
        /// array, so a saved index would quietly point at a different glycan if the .gdb were edited.
        ///
        /// List&lt;(string, string)&gt; is deliberate -- MetaMorpheusTask.tomlConfig already registers a
        /// converter for exactly this type, so it round-trips with no new serialization code.
        ///
        /// The key is composition-based, so two structural isomers in one user-built structure database
        /// (e.g. (N(H(N))) and (N(H)(N))) share it: ticking either searches both. No shipped database holds
        /// isomers; separating them needs a structure-aware key, which is follow-up work.
        /// </remarks>
        public List<(string, string)> SelectedGlycans { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }
        public bool OxoniumIonFilt { get; set; }

        /// <summary>
        /// With <see cref="OxoniumIonFilt"/> on, a custom-monosaccharide diagnostic ion (MonosaccharidesCustom.tsv) counts as
        /// observed only above this fraction of the HexNAc 138.055 ion. Within [0, 1]; default 0.02 (2%). The built-in NeuAc
        /// and HexNAc(1)Hex(1) rules always use 0.02. TOMLs written before this existed omit it and get the default.
        /// </summary>
        public double CustomOxoniumThreshold
        {
            get => _customOxoniumThreshold;
            set
            {
                if (!GlycoPeptides.IsValidCustomOxoniumThreshold(value))
                {
                    throw new EngineLayer.MetaMorpheusException(
                        $"CustomOxoniumThreshold must be between 0 and 1 (a fraction of the 138.055 ion, e.g. 0.02 for 2%), but was {value}.");
                }
                _customOxoniumThreshold = value;
            }
        }
        private double _customOxoniumThreshold = GlycoPeptides.DefaultCustomOxoniumThreshold;

        /// <summary>
        /// The precursor mass offsets a glycopeptide may explain, as in the Classic and Modern searches. A large glycopeptide's
        /// monoisotopic peak is small, so its precursor is often picked one or more 13C peaks heavy (or light); an isotope-offset
        /// type searches each offset as the observed precursor mass minus the offset, and each match keeps the offset's notch.
        /// Only <see cref="SupportedMassDiffAcceptorTypes"/> are allowed. Default <see cref="MassDiffAcceptorType.Exact"/>, the
        /// behaviour before this setting existed; TOMLs written before it omit it and get the default.
        /// </summary>
        public MassDiffAcceptorType MassDiffAcceptorType
        {
            get => _massDiffAcceptorType;
            set
            {
                if (!SupportedMassDiffAcceptorTypes.Contains(value))
                {
                    throw new EngineLayer.MetaMorpheusException(
                        $"Glyco search supports the precursor mass difference acceptors {string.Join(", ", SupportedMassDiffAcceptorTypes)}, but was {value}.");
                }
                _massDiffAcceptorType = value;
            }
        }
        private MassDiffAcceptorType _massDiffAcceptorType = MassDiffAcceptorType.Exact;

        /// <summary>
        /// The acceptors that only allow isotope (missed monoisotopic) offsets of the observed precursor mass. Open, interval and
        /// custom acceptors would let any mass difference through as a "glycan" and are not supported.
        /// </summary>
        public static readonly IReadOnlyList<MassDiffAcceptorType> SupportedMassDiffAcceptorTypes = new[]
        {
            MassDiffAcceptorType.Exact,
            MassDiffAcceptorType.OneMM,
            MassDiffAcceptorType.TwoMM,
            MassDiffAcceptorType.ThreeMM,
            MassDiffAcceptorType.PlusOrMinusThreeMM,
        };

        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }

        /// <summary>
        /// Glycan boxes (the summed glycans placed on one peptide) heavier than this, in Da, are never built or searched.
        /// Applies to O-, N- and N+O searches; in an N-glycan search the box is the single N-glycan.
        /// </summary>
        public double MaximumGlycanBoxMass { get; set; }

        /// <summary>
        /// Most spectra files searched at the same time. 0 decides from the thread budget (MaxThreadsToUsePerFile, divided across the
        /// files) and free memory; 1 searches files one after another.
        /// </summary>
        public int MaximumSpectraFilesInParallel { get; set; }

        public bool DoParsimony { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool ModPeptidesAreDifferent { get; set; }
        
        //quantification options
        public bool DoQuantification { get; set; }
        public bool DoMbrAnalysis { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool Normalize { get; set; }

        //output options
        public bool WriteIndividualFiles { get; set; }
        public bool WriteDecoys { get; set; }
        public bool WriteContaminants { get; set; }
        public bool WriteSpectrumLibrary { get; set; }
        public bool WritePrunedDataBase { get; set; }
        public bool DisposeOfFileWhenDone { get; set; }

        public Dictionary<string, int> ModsToWriteSelection { get; set; }
    }
}