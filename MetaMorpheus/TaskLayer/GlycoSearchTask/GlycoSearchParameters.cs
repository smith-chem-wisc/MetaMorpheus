using EngineLayer.GlycoSearch;
using System.Collections.Generic;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class GlycoSearchParameters : SearchParameters
    {
        public GlycoSearchParameters()
        {
            OGlycanDatabasefile = "OGlycan.gdb";
            NGlycanDatabasefile = "NGlycan.gdb";
            GlycoSearchType = GlycoSearchType.OGlycanSearch;
            OxoniumIonFilt = true;
            DecoyType = DecoyType.Reverse;
            GlycoSearchTopNum = 50;
            MaximumOGlycanAllowed = 4;
            MaximumGlycanBoxMass = EngineLayer.GlycanBox.DefaultMaximumGlycanBoxMass;
            DoParsimony = true;
            NoOneHitWonders = false;
            ModPeptidesAreDifferent = false;

            //quantification options
            DoLabelFreeQuantification = false;
            MatchBetweenRuns = true;
            QuantifyPpmTol = 5;
            Normalize = false;

            //output options
            WriteIndividualFiles = false;
            WriteDecoys = true;
            WriteContaminants = true;
            WriteSpectralLibrary = false;
            DisposeOfFileWhenDone = true;
            WritePrunedDatabase = false;

            ModsToWriteSelection = SearchParameters.DefaultModsToWriteSelection();
        }
        public string OGlycanDatabasefile { get; set; }
        public string NGlycanDatabasefile { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }
        public bool OxoniumIonFilt { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }

        /// <summary>
        /// Glycan boxes (the summed glycans placed on one peptide) heavier than this, in Da, are never built or searched.
        /// Applies to O-, N- and N+O searches; in an N-glycan search the box is the single N-glycan.
        /// </summary>
        public double MaximumGlycanBoxMass { get; set; }
    }
}
