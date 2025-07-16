using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer.GlycoSearch;
using EngineLayer.ModernSearch;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.ModSearch
{
    public class ModSearchEngine : ModernSearchEngine
    {
        protected readonly List<GlycoSpectralMatch>[] GlobalMsms; // Global mod spectral matches

        // For general modification settings
        private readonly int TopN; // DDA top Peak number
        private readonly int MaxModNumber;
        private readonly Tolerance PrecursorSearchMode;
        private readonly MassDiffAcceptor ProductSearchMode;
        private readonly List<int>[] SecondFragmentIndex;

        // For glyco settings
        private readonly bool OxoniumIonFilter;
        private readonly GlycoSearchType GlycoSearchType;
        private readonly string OGlycanDatabaseFile;
        private readonly string NGlycanDatabaseFile;

        public ModSearchEngine(List<GlycoSpectralMatch>[] globalCsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex,
            List<int>[] fragmentIndex, List<int>[] secondFragmentIndex, int currentPartition, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int modSearchTopNum, int maxModNum, bool oxoniumIonFilter, List<string> nestedIds)
            : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, commonParameters, fileSpecificParameters, new OpenSearchMode(), 0, nestedIds)
        {
            this.GlobalMsms = globalCsms;
            this.SecondFragmentIndex = secondFragmentIndex;
            this.TopN = modSearchTopNum;
            this.MaxModNumber = maxModNum;
            this.OxoniumIonFilter = oxoniumIonFilter;
            this.GlycoSearchType = glycoSearchType;
            this.OGlycanDatabaseFile = oglycanDatabase;
            this.NGlycanDatabaseFile = nglycanDatabase;

            this.SecondFragmentIndex = secondFragmentIndex;
            this.ProductSearchMode = new SinglePpmAroundZeroSearchMode(20); //For ScanOxoniumIonFilter only
            this.PrecursorSearchMode = commonParameters.PrecursorMassTolerance;

            //Load glycan databases and build the modBox 

        }



    }
}
