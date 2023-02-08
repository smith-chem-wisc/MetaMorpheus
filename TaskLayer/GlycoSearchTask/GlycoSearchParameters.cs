using EngineLayer;
using System.Linq;
using UsefulProteomicsDatabases;
using EngineLayer.GlycoSearch;

namespace TaskLayer
{
    public class GlycoSearchParameters
    {
        public GlycoSearchParameters()
        {

        }
        public GlycoSearchParameters(string oGlycanDatabaseFile = "OGlycan.gdb", string nGlycanDatabaseFile = "NGlycan.gdb", GlycoSearchType glycoSearchType = GlycoSearchType.OGlycanSearch, 
            bool oxoniumIonFit = true, DecoyType decoyType = DecoyType.Reverse, int glycoSearchTopNum = 50, int maxOglycansAllowed = 4, 
            bool doParsimony = false, bool noOneHitWonders = false, bool modPeptidesAreDifferent = false, bool writeIndividualFiles = false, bool writeDecoys = false, bool writeContaminants = false)
        {
            OGlycanDatabasefile = oGlycanDatabaseFile;
            NGlycanDatabasefile = nGlycanDatabaseFile;
            GlycoSearchType = glycoSearchType;
            OxoniumIonFilt = oxoniumIonFit;
            DecoyType = decoyType;
            GlycoSearchTopNum = glycoSearchTopNum;
            MaximumOGlycanAllowed = maxOglycansAllowed;

            DoParsimony = doParsimony;
            NoOneHitWonders = noOneHitWonders;
            ModPeptidesAreDifferent = modPeptidesAreDifferent;

            WriteIndividualFiles = writeIndividualFiles;
            WriteDecoys = writeDecoys;
            WriteContaminants = writeContaminants;
        }
        public string OGlycanDatabasefile { get; set; }
        public string NGlycanDatabasefile { get; set; }
        public GlycoSearchType GlycoSearchType { get; set; }
        public bool OxoniumIonFilt { get; set; }
        public DecoyType DecoyType { get; set; }
        public int GlycoSearchTopNum { get; set; }
        public int MaximumOGlycanAllowed { get; set; }

        public bool DoParsimony { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool ModPeptidesAreDifferent { get; set; }

        public bool WriteIndividualFiles { get; set; }

        public bool WriteDecoys { get; set; }
        public bool WriteContaminants { get; set; }
    }
}