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
                GlycanBox.GlobalOGlycans = GlycoSearchEngine.LoadGlycanDatabase(GlobalVariables.OGlycanDatabasePaths, oglycanDatabase, "O-glycan", true);
                GlycanBox.OGlycanBoxes = GlycoSearchEngine.CheckedGlycanBoxes( //generate glycan box for O-glycan search
                    GlycanBox.BuildOGlycanBoxes(maxOGlycanNum, false, maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray(),
                    $"the O-glycan database '{oglycanDatabase}'", maxGlycanBoxMass);
                GlycoSpectralMatch.GlycanBoxes = GlycanBox.OGlycanBoxes;
                return new GlycanSearchSpace(glycoSearchType, GlycanBox.OGlycanBoxes, null);
            }
            if (glycoSearchType == GlycoSearchType.NGlycanSearch) //because the there is only one glycan in N-glycanpeptide, so we don't need to build the n-glycanBox here.
            {
                // The single N-glycan is the whole box here, so the box mass cap applies to each glycan on its own.
                var nGlycans = GlycoSearchEngine.LoadGlycanDatabase(GlobalVariables.NGlycanDatabasePaths, nglycanDatabase, "N-glycan", false)
                    .Where(p => (double)p.Mass / 1E5 <= maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray();
                // LoadGlycanDatabase refused an empty file, but the cap can still empty it here, and the search
                // would then skip every scan and report nothing. The O and N+O paths refuse this in
                // CheckedGlycanBoxes; this is the same refusal for the path that builds no boxes.
                if (nGlycans.Length == 0)
                {
                    throw new MetaMorpheusException(
                        $"No glycan in the N-glycan database '{nglycanDatabase}' is within the maximum glycan box mass of {maxGlycanBoxMass} Da, " +
                        "so there is nothing to search for. Raise that maximum, or choose a database of lighter glycans.");
                }
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
                return new GlycanSearchSpace(glycoSearchType, null, nGlycans);
            }
            if (glycoSearchType == GlycoSearchType.N_O_GlycanSearch) //search both N-glycan and O-glycan is still not tested and build completely yet.
            {
                GlycanBox.GlobalOGlycans = GlycoSearchEngine.LoadGlycanDatabase(GlobalVariables.OGlycanDatabasePaths, oglycanDatabase, "O-glycan", true);
                GlycanBox.GlobalNGlycans = new Dictionary<int, Glycan>();
                // For N-glycan, we use negative index to distinguish with O-glycan.
                var nGlycans = GlycoSearchEngine.LoadGlycanDatabase(GlobalVariables.NGlycanDatabasePaths, nglycanDatabase, "N-glycan", false).OrderBy(p => p.Mass);
                int indexForNGlycan = -1;
                foreach (var nGlycan in nGlycans)
                {
                    GlycanBox.GlobalNGlycans.Add(indexForNGlycan, nGlycan);
                    indexForNGlycan--;
                }

                GlycanBox.NOGlycanBoxes = GlycoSearchEngine.CheckedGlycanBoxes(
                    GlycanBox.BuildNOGlycanBoxes(maxOGlycanNum, false, maxGlycanBoxMass).OrderBy(p => p.Mass).ToArray(),
                    $"the O-glycan database '{oglycanDatabase}' and the N-glycan database '{nglycanDatabase}'", maxGlycanBoxMass);
                GlycoSpectralMatch.GlycanBoxes = GlycanBox.NOGlycanBoxes;
                //TO THINK: Glycan Decoy database.
                //DecoyGlycans = Glycan.BuildTargetDecoyGlycans(NGlycans);
                return new GlycanSearchSpace(glycoSearchType, GlycanBox.NOGlycanBoxes, null);
            }
            return new GlycanSearchSpace(glycoSearchType, null, null);
        }
    }
}
