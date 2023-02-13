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