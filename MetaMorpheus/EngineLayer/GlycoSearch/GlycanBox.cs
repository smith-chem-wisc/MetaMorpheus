using System.Collections.Generic;
using System.Linq;
using System;
using System.Threading;
using EngineLayer.GlycoSearch;
using MzLibUtil;

namespace EngineLayer
{

    /// <summary>
    /// A defined combination of glycans to modify on one peptide. Ex. if we have 3 glycans on one peptide (g1,g2,g3), the GlycanBoxMass is the sum of the three glycans.(glycanBox: [g1,g2,g3])
    /// </summary>
    public class GlycanBox:ModBox
    {
        /// <summary>
        /// The global list of O-glycans loaded from the glycan database file.
        /// </summary>
        public static Glycan[] GlobalOGlycans { get; set; }
        public static Dictionary<int, Glycan> GlobalNGlycans { get; set; }

        public int[] OGlycanIds { get; set; }
        public int NGlycanId { get; set; }

        /// <summary>
        /// All possible child glycan box combinations derived from this glycan box.
        /// </summary>
        /// <remarks>
        /// Boxes from <see cref="BuildOGlycanBoxes(int, bool, double)"/> and <see cref="BuildNOGlycanBoxes(int, bool, double)"/>
        /// build these on first read. An N+O search makes hundreds of thousands of boxes with about 16 child boxes each, and
        /// the localization graph reads them only for boxes that match a precursor mass. Building them all up front was most
        /// of the box-building time and kept millions of objects alive for the whole search. Every read returns the same
        /// array, which <see cref="LocalizationGraph"/> relies on. A box made directly with a constructor has none, and an
        /// assigned array is returned as assigned.
        /// </remarks>
        public GlycanBox[] ChildGlycanBoxes
        {
            get
            {
                if (_childGlycanBoxes == null && _childBoxBuilder != ChildBoxBuilder.None)
                {
                    // Racing first reads may each build the array; only one is published and every read returns it.
                    LazyInitializer.EnsureInitialized(ref _childGlycanBoxes, BuildChildGlycanBoxes);
                }
                return _childGlycanBoxes;
            }
            set
            {
                _childGlycanBoxes = value;
            }
        }

        private GlycanBox[] _childGlycanBoxes;

        /// <summary>
        /// Which child builder a box made by the box builders uses on first read of <see cref="ChildGlycanBoxes"/>.
        /// </summary>
        private enum ChildBoxBuilder : byte
        {
            None,
            OGlycan,
            NOGlycan,
        }

        private ChildBoxBuilder _childBoxBuilder;

        /// <summary>
        /// True once this box holds its child boxes, whether built on read or assigned.
        /// </summary>
        internal bool HasBuiltChildGlycanBoxes => Volatile.Read(ref _childGlycanBoxes) != null;

        private GlycanBox[] BuildChildGlycanBoxes()
        {
            return _childBoxBuilder == ChildBoxBuilder.OGlycan
                ? BuildChildOGlycanBoxes(NumberOfMods, ModIds, TargetDecoy).ToArray()
                : BulidChildNOBoxes(NumberOfMods, ModIds, TargetDecoy).ToArray();
        }

        /// <summary>
        /// The global collection of all possible O-glycan boxes.
        /// </summary>
        public static GlycanBox[] OGlycanBoxes { get; set; }

        public static GlycanBox[] NOGlycanBoxes { get; set; }

        /// <summary>
        /// The summed glycan composition for this box, where each element represents the count of a specific monosaccharide type.
        /// </summary>
        public byte[] Kind { get; private set; }

            //TO DO: Decoy O-glycan can be created, but the results need to be reasoned.
            //public static int[] SugarShift = new int[]{ -16205282, -20307937, -29109542, -14605791, -30709033, -15005282, -36513219, -40615874, 16205282, 20307937, 29109542, 14605791, 30709033, 15005282, 36513219, 40615874 };
            private readonly static int[] SugarShift = new int[] //still unclear about the shift...
            {
                7103710, 10300920, 11502690, 12904260, 14706840, 5702150, 13705890, 12809500, 11308410, 13104050,
                11404290, 9705280, 12805860, 15610110, 8703200, 10104770, 9906840, 18607930, 16306330,
                -7103710, -10300920, -11502690, -12904260, -14706840, -5702150, -13705890, -12809500, -11308410, -13104050,
                -11404290, -9705280, -12805860, -15610110, -8703200, -10104770, -9906840, -18607930, -16306330,

            };

        /// <summary>
        /// Use the glycan from database to create all possible combination glycan set into GlycanBox. 
        /// </summary>
        /// <param name="maxNum"> The maxNum is maximum glycans allowed on one peptides </param>
        /// <returns> The glycanBox collection, glycanBox[]</returns>
        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum)
        {
            return BuildOGlycanBoxes(maxNum, false);
        }
        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum, bool buildDecoy)
        {
            return BuildOGlycanBoxes(maxNum, buildDecoy, double.MaxValue);
        }

        /// <summary>
        /// Default for <see cref="BuildOGlycanBoxes(int, bool, double)"/> and <see cref="BuildNOGlycanBoxes(int, bool, double)"/> as used by the search, in Da.
        /// </summary>
        public const double DefaultMaximumGlycanBoxMass = 4000;

        /// <param name="maxBoxMass"> Boxes heavier than this (Da) are skipped before their child boxes are built. </param>
        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum, bool buildDecoy, double maxBoxMass)
        {

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalOGlycans.Length), i))
                {
                    GlycanBox glycanBox = new GlycanBox(idCombine.ToArray());
                    if (glycanBox.Mass > maxBoxMass)
                    {
                        continue;
                    }
                    glycanBox.TargetDecoy = true;
                    glycanBox._childBoxBuilder = ChildBoxBuilder.OGlycan;

                    yield return glycanBox;

                    if (buildDecoy)
                    {
                        GlycanBox glycanBox_decoy = new GlycanBox(idCombine.ToArray(),false); // decoy glycanBox
                        glycanBox_decoy.TargetDecoy = false;
                        glycanBox_decoy._childBoxBuilder = ChildBoxBuilder.OGlycan;
                        yield return glycanBox_decoy;
                    }
                }
            }
        }

        /// <summary>
        /// Generate all possible child/fragment box of the specific glycanBox. The childBoxes is uesd for LocalizationGraph.
        /// </summary>
        /// <param name="maxNum"></param>
        /// <param name="glycanIds"> The glycanBox, ex. [0,0,1] means glycan0 + glycan0 + glycan1 </param>
        /// <param name="targetDecoy"></param>
        /// <returns> The ChildBox collection, ChildBox[] </returns>
        public static IEnumerable<GlycanBox> BuildChildOGlycanBoxes(int maxNum, int[] glycanIds, bool targetDecoy = true)
        {
            yield return new GlycanBox(new int[0], targetDecoy);
            HashSet<int[]> seen = new HashSet<int[]>(IdSequenceComparer.Instance);
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombs(Enumerable.Range(0, maxNum), i)) //get all combinations of glycans on the peptide, ex. we have three glycosite and three glycan maybe on that (A,B,C)
                {                                                                           //the combination of glycans on the peptide can be (A),(A+B),(A+C),(B+C),(A+B+C) totally six
                    int[] ids = new int[i];
                    int n = 0;
                    foreach (var id in idCombine)
                    {
                        ids[n++] = glycanIds[id];
                    }

                    if (seen.Add(ids))
                    {
                        GlycanBox glycanBox = new GlycanBox(ids, targetDecoy);

                        yield return glycanBox;
                    }

                }
            }
        }

        public static IEnumerable<GlycanBox> BuildNOGlycanBoxes(int maxOGlycanNum, bool buildDecoy)
        {
            return BuildNOGlycanBoxes(maxOGlycanNum, buildDecoy, double.MaxValue);
        }

        /// <param name="maxBoxMass"> Boxes heavier than this (Da) are skipped before their child boxes are built. </param>
        public static IEnumerable<GlycanBox> BuildNOGlycanBoxes(int maxOGlycanNum, bool buildDecoy, double maxBoxMass)
        {
            int[] oGlycansIds;

            foreach (var box in AllNGlycanOnlyBox(buildDecoy, maxBoxMass))
            {
                yield return box;
            }

            for (int i = 1; i <= maxOGlycanNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalOGlycans.Length), i))
                {
                    oGlycansIds = idCombine.ToArray();
                    // Most O/N pairings exceed the cap once the O-glycans alone are heavy, so screen on the summed
                    // glycan masses (with a 1 Da margin) before allocating a box. The exact test below decides.
                    double oGlycanMass = oGlycansIds.Sum(id => (double)GlobalOGlycans[id].Mass) / 1E5;
                    if (oGlycanMass > maxBoxMass + 1)
                    {
                        continue;
                    }
                    List<int> keptNGlycanIds = new List<int>();
                    for (int j = 0; j < GlobalNGlycans.Count + 1; j++)
                    {
                        // the index for N-glycan will be start from -1, -2, -3...
                        // nglycanId = 0 means no glycan on the peptide.
                        int nglycanId = - j;
                        if (nglycanId != 0 && oGlycanMass + (double)GlobalNGlycans[nglycanId].Mass / 1E5 > maxBoxMass + 1)
                        {
                            continue;
                        }
                        GlycanBox glycanBox = new GlycanBox(oGlycansIds, nglycanId, true);
                        if (glycanBox.Mass > maxBoxMass)
                        {
                            continue;
                        }
                        keptNGlycanIds.Add(nglycanId);
                        glycanBox.TargetDecoy = true;
                        glycanBox._childBoxBuilder = ChildBoxBuilder.NOGlycan;
                        yield return glycanBox;
                    }

                    if (buildDecoy)
                    {
                        // A decoy is built only where its target composition survived the mass cap.
                        foreach (int nglycanId in keptNGlycanIds)
                        {
                            GlycanBox glycanBox_decoy = new GlycanBox(oGlycansIds, nglycanId,false); // decoy glycanBox
                            glycanBox_decoy.TargetDecoy = false;
                            glycanBox_decoy._childBoxBuilder = ChildBoxBuilder.NOGlycan;
                            yield return glycanBox_decoy;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Generate all possible glycan boxes that contain only one N-glycan.
        /// </summary>
        /// <param name="buildDecoy"></param>
        /// <returns></returns>
        private static IEnumerable<GlycanBox> AllNGlycanOnlyBox(bool buildDecoy, double maxBoxMass)
        {
            List<int> keptNGlycanIds = new List<int>();
            // first consider the case when there is no N-glycan in the box.
            for (int j = 1; j < GlobalNGlycans.Count + 1; j++)
            {
                // the index for N-glycan will be start from -1, -2, -3...
                // nglycanId = 0 means no glycan on the peptide.
                int nglycanId = -j;
                GlycanBox glycanBox = new GlycanBox(null, nglycanId, true);
                if (glycanBox.Mass > maxBoxMass)
                {
                    continue;
                }
                keptNGlycanIds.Add(nglycanId);
                glycanBox.TargetDecoy = true;
                glycanBox._childBoxBuilder = ChildBoxBuilder.NOGlycan;
                yield return glycanBox;
            }
            if (buildDecoy)
            {
                foreach (int nglycanId in keptNGlycanIds)
                {
                    GlycanBox glycanBox_decoy = new GlycanBox(null, nglycanId, false); // decoy glycanBox
                    glycanBox_decoy.TargetDecoy = false;
                    glycanBox_decoy._childBoxBuilder = ChildBoxBuilder.NOGlycan;
                    yield return glycanBox_decoy;
                }
            }
        }

        public static IEnumerable<GlycanBox> BulidChildNOBoxes(int maxNum, int[] glycanIds, bool isTarget = true)
        {
            yield return new GlycanBox(new int[0], isTarget);
            HashSet<int[]> seen = new HashSet<int[]>(IdSequenceComparer.Instance);

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombs(Enumerable.Range(0, maxNum), i)) //get all combinations of glycans on the peptide, ex. we have three glycosite and three glycan maybe on that (A,B,C)
                {                                                                           //the combination of glycans on the peptide can be (A),(A+B),(A+C),(B+C),(A+B+C) totally six
                    int[] ids = new int[i];
                    int n = 0;
                    int oGlycanCount = 0;
                    foreach (var id in idCombine)
                    {
                        ids[n] = glycanIds[id];
                        if (ids[n] > -1)
                        {
                            oGlycanCount++;
                        }
                        n++;
                    }

                    if (seen.Add(ids))
                    {
                        int[] oGlycanIds = new int[oGlycanCount];
                        // If there is no N-glycan on the peptide, the nGlycanid = 0; otherwise it is the first negative id.
                        int nGlycanid = 0;
                        int o = 0;
                        foreach (int id in ids)
                        {
                            if (id > -1)
                            {
                                oGlycanIds[o++] = id;
                            }
                            else if (nGlycanid == 0)
                            {
                                nGlycanid = id;
                            }
                        }
                        GlycanBox glycanBox = new GlycanBox(oGlycanIds,nGlycanid, isTarget);

                        yield return glycanBox;
                    }

                }
            }
        }

        /// <summary>
        /// Constructor of GlycanBox.
        /// </summary>
        /// <param name="ids"> The glycanBox composition, each number represent one glycan index in the database</param>
        /// <param name="targetDecoy"></param>
        public GlycanBox(int[] ids, bool Istarget = true):base(ids)
        {
            byte[] kind = new byte[Glycan.KindCapacity];
            foreach (var id in ModIds) //ModIds is the same as ids.
            {
                for (int i = 0; i < kind.Length; i++)   
                {
                    kind[i] += GlobalOGlycans[id].Kind[i]; //kind is the sum of all glycan Kind in the Box.
                }
            }
            Kind = kind;

            if (Istarget)
            {
                Mass = (double)Glycan.GetMass(Kind) / 1E5;
            }
            else
            {
                Random random = new Random();
                int shiftInd = random.Next(SugarShift.Length);
                Mass = (double)(Glycan.GetMass(Kind) + SugarShift[shiftInd]) / 1E5;
            }
        }

        public GlycanBox(int[] oGlycanIds, int nGlycanIds, bool Istarget = true) : base(oGlycanIds, nGlycanIds)
        {
            OGlycanIds = oGlycanIds;
            NGlycanId = nGlycanIds;
            byte[] kind = new byte[Glycan.KindCapacity];

            if (!oGlycanIds.IsNullOrEmpty())
            {
                foreach (var id in oGlycanIds) //ModIds is the same as ids.
                {
                    for (int i = 0; i < kind.Length; i++)
                    {
                        kind[i] += GlobalOGlycans[id].Kind[i]; //kind is the sum of all glycan Kind in the Box.
                    }
                }
            }

            if (nGlycanIds != 0)
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    kind[i] += GlobalNGlycans[nGlycanIds].Kind[i];
                }
            }
            Kind = kind;

            if (Istarget)
            {
                Mass = (double)Glycan.GetMass(Kind) / 1E5;
            }
            else
            {
                Random random = new Random();
                int shiftInd = random.Next(SugarShift.Length);
                Mass = (double)(Glycan.GetMass(Kind) + SugarShift[shiftInd]) / 1E5;
            }
        }

        /// <summary>
        /// Equality of glycan id sequences by content and order, as comparing their comma-joined strings did.
        /// </summary>
        private sealed class IdSequenceComparer : IEqualityComparer<int[]>
        {
            public static readonly IdSequenceComparer Instance = new IdSequenceComparer();

            public bool Equals(int[] x, int[] y)
            {
                return ReferenceEquals(x, y) || (x != null && y != null && x.AsSpan().SequenceEqual(y));
            }

            public int GetHashCode(int[] ids)
            {
                var hash = new HashCode();
                foreach (int id in ids)
                {
                    hash.Add(id);
                }
                return hash.ToHashCode();
            }
        }

        /// <summary>
        /// The motif (e.g. "S", "T", "Nxs", "Nxt") of a glycan id. Non-negative ids are O-glycans, negative ids are N-glycans.
        /// </summary>
        internal static string MotifOf(int modId)
        {
            return modId >= 0 ? GlobalOGlycans[modId].Target.ToString() : GlobalNGlycans[modId].Target.ToString();
        }

        private MotifCount _motifCount;
        private GlycanBoxLocalizationCache _localizationCache;

        /// <summary>
        /// How many glycans in this box need each motif. Built on first use; a box is shared by all search threads.
        /// </summary>
        internal MotifCount GetMotifCount()
        {
            return LazyInitializer.EnsureInitialized(ref _motifCount, () => new MotifCount(ModIds));
        }

        /// <summary>
        /// Everything the localization graph derives from this box's child boxes alone. Built on first use.
        /// </summary>
        internal GlycanBoxLocalizationCache GetLocalizationCache()
        {
            return LazyInitializer.EnsureInitialized(ref _localizationCache, () => new GlycanBoxLocalizationCache(this));
        }

        public string GlycanIdString // the composition of glycanBox. Example: [1,2,3] means glycan1 + glycan2 + glycan3 are on the peptide.
        {
            get
            {
                return string.Join(",", ModIds.Select(p => p.ToString()));
            }
        }
    }
}
