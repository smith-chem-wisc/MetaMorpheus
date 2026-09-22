using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.GlycoSearch
{
    /// <summary>
    /// The glycans and glycan boxes a glyco search tries against every peptide. They depend only on the glycan databases and the
    /// search settings, never on a spectra file, so a task builds them once and hands the same instance to every file's engine.
    /// </summary>
    public sealed class GlycanSearchSpace
    {
        public GlycoSearchType GlycoSearchType { get; }

        /// <summary> The glycan boxes, sorted by mass, for O-glycan and N+O searches; null for an N-glycan search. </summary>
        public GlycanBox[] GlycanBoxes { get; }

        /// <summary> The N-glycans, sorted by mass, for an N-glycan search; null otherwise. </summary>
        public Glycan[] NGlycans { get; }

        /// <summary> GlycanBoxes[i].Mass. </summary>
        internal double[] GlycanBoxMasses { get; }

        /// <summary> NGlycans[i].Mass in Da. </summary>
        internal double[] NGlycanMasses { get; }

        private GlycanSearchSpace(GlycoSearchType glycoSearchType, GlycanBox[] glycanBoxes, Glycan[] nGlycans)
        {
            GlycoSearchType = glycoSearchType;
            GlycanBoxes = glycanBoxes;
            NGlycans = nGlycans;
            GlycanBoxMasses = glycanBoxes?.Select(p => p.Mass).ToArray();
            NGlycanMasses = nGlycans?.Select(p => (double)p.Mass / 1E5).ToArray();
        }

        /// <summary>
        /// Loads the glycan databases and builds the glycan boxes. Also sets the static glycan state the rest of the glyco code reads
        /// (GlycanBox.GlobalOGlycans, GlobalNGlycans, OGlycanBoxes, NOGlycanBoxes and GlycoSpectralMatch.GlycanBoxes), so call it
        /// before, not while, any search that reads them.
        /// </summary>
        public static GlycanSearchSpace Build(string oglycanDatabase, string nglycanDatabase, GlycoSearchType glycoSearchType, int maxOGlycanNum,
            double maxGlycanBoxMass = GlycanBox.DefaultMaximumGlycanBoxMass)
        {
            if (glycoSearchType == GlycoSearchType.OGlycanSearch) //if we do the O-glycan search, we need to load the O-glycan database and generate the glycoBox.
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.OGlycanBoxes = GlycanBox.BuildOGlycanBoxes(maxOGlycanNum, false, maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray(); //generate glycan box for O-glycan search
                GlycoSpectralMatch.GlycanBoxes = GlycanBox.OGlycanBoxes;
                return new GlycanSearchSpace(glycoSearchType, GlycanBox.OGlycanBoxes, null);
            }
            if (glycoSearchType == GlycoSearchType.NGlycanSearch) //because the there is only one glycan in N-glycanpeptide, so we don't need to build the n-glycanBox here.
            {
                // The single N-glycan is the whole box here, so the box mass cap applies to each glycan on its own.
                var nGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == nglycanDatabase).First(), true, false)
                    .Where(p => (double)p.Mass / 1E5 <= maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray();
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
                return new GlycanSearchSpace(glycoSearchType, null, nGlycans);
            }
            if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch) //search both N-glycan and O-glycan is still not tested and build completely yet.
            {
                GlycanBox.GlobalOGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.OGlycanDatabasePaths.Where(p => System.IO.Path.GetFileName(p) == oglycanDatabase).First(), true, true).ToArray();
                GlycanBox.GlobalNGlycans = new Dictionary<int, Glycan>();
                // For N-glycan, we use negative index to distinguish with O-glycan.
                var nGlycans = GlycanDatabase.LoadGlycan(GlobalVariables.NGlycanDatabasePaths.First(p => System.IO.Path.GetFileName(p) == nglycanDatabase),
                        true, false).OrderBy(p => p.Mass);
                int indexForNGlycan = -1;
                foreach (var nGlycan in nGlycans)
                {
                    GlycanBox.GlobalNGlycans.Add(indexForNGlycan, nGlycan);
                    indexForNGlycan--;
                }

                GlycanBox.NOGlycanBoxes = GlycanBox.BuildNOGlycanBoxes(maxOGlycanNum, false, maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray();
                GlycoSpectralMatch.GlycanBoxes = GlycanBox.NOGlycanBoxes;
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
                return new GlycanSearchSpace(glycoSearchType, GlycanBox.NOGlycanBoxes, null);
            }
            return new GlycanSearchSpace(glycoSearchType, null, null);
        }
    }
}
