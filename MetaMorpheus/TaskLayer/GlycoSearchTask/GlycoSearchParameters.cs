using EngineLayer.GlycoSearch;
using System.Collections.Generic;
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
        /// </remarks>
        public List<(string, string)> SelectedGlycans { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }
        public bool OxoniumIonFilt { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }

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