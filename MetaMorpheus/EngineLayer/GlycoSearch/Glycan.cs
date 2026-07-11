using System.Collections.Generic;
using System.Linq;
using Chemistry;
using System;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using Omics.Modifications;

namespace EngineLayer
{
    public class GlycanIon
    {
        public GlycanIon(string ionStruct, int ionMass, byte[] ionKind, int lossIonMass)
        {
            IonStruct = ionStruct;     // Always set null, deprecated.
            IonMass = ionMass;
            IonKind = ionKind;
            LossIonMass = lossIonMass; // Neutral loss mass = Glycan.Mass - IonMass
        }
        public string IonStruct { get; set; }
        public int IonMass { get; set; }
        public int LossIonMass { get; set; }
        public byte[] IonKind { get; set; }
    }

    /// <summary>
    /// Glycan represents a glycan modification, which can be O-glycan or N-glycan.
    /// Included information like glycan structure, mass, kind, ions, and type.
    /// </summary>
    public class Glycan :  Modification
    {
        public Glycan(string struc, int mass, byte[] kind, List<GlycanIon> ions, bool decoy, string motif, GlycanType type = GlycanType.O_glycan) 
            : base( _monoisotopicMass: mass / 1E5, _locationRestriction: "Anywhere.") // Divide the mass by 1E5 to convert it to monoisotopic mass in Daltons (Da).
        {
            // Glycan Properties
            Struc = struc;
            Mass = mass; // Glycan mass is stored as an integer scaled by 1e5 to improve performance. Divide by 1e5 to obtain the monoisotopic mass in Daltons (Da).
            Kind = kind;
            Ions = ions;
            Decoy = decoy;
            Type = type;

            // Modification Properties
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            // Generate the neural loss and diagnostic ions for O_glycan.
            if (type == GlycanType.O_glycan)
            {
                ModificationType = "O-linked glycosylation"; // Set the modification type.
                if (Ions != null)
                {
                    List<double> lossMasses = Ions.Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList();
                    neutralLosses.Add(DissociationType.HCD, lossMasses);
                    neutralLosses.Add(DissociationType.CID, lossMasses);
                    neutralLosses.Add(DissociationType.EThcD, lossMasses);
                }
            }

            // Generate the neural loss and diagnostic ions for N_glycan.
            else if (type == GlycanType.N_glycan)
            {
                ModificationType = "N-linked glycosylation"; // Set the modification type.
                if (Ions != null)
                {
                    List<double> lossMasses = Ions.Where(p=>p.IonMass < 57000000).Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList();
                    neutralLosses.Add(DissociationType.HCD, lossMasses);
                    neutralLosses.Add(DissociationType.CID, lossMasses);
                    neutralLosses.Add(DissociationType.EThcD, lossMasses);
                }
            }

            Dictionary<DissociationType, List<double>> diagnosticIons = new Dictionary<DissociationType, List<double>>();
            diagnosticIons.Add(DissociationType.HCD, GlycanDiagnosticIons.Select(p => (double)p / 1E5).ToList()); // Divided by 1E5 to convert mass(int) to monoisotopic mass(double).
            diagnosticIons.Add(DissociationType.CID, GlycanDiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.EThcD, GlycanDiagnosticIons.Select(p => (double)p / 1E5).ToList());
            ModificationMotif.TryGetMotif(motif, out ModificationMotif finalMotif); //TO DO: only one motif can be write here.
            var id = Glycan.GetKindString(Kind);

            OriginalId = id; // Set the original ID to the glycan kind string, which is unique for each glycan.
            Target = finalMotif; // Set the target motif for the modification.
            NeutralLosses = neutralLosses; // Set the neutral losses for the modification.
            base.DiagnosticIons = diagnosticIons; // Set the diagnostic ions for the modification.

            if (OriginalId != null)
            {
                IdWithMotif = OriginalId + " on " + Target.ToString();
                OriginalId = OriginalId;
            }
            else
                OriginalId = OriginalId;
        }

        /// <summary>
        /// In this constructor, we will generate the glycan only by the glycan kind and type.
        /// So there is no ions information, and the diagnostic ions will not be generated.
        /// </summary>
        /// <param name="kind"></param>
        /// <param name="type"></param>
        public Glycan(byte[] kind, string motif, GlycanType type)
            : this(null, GetMass(kind), kind, null, false, motif, type)
        {
        }

        /// <summary>
        /// Glycan ID, which is the index of glycan in the glycan database.
        /// </summary>
        public int GlyId { get; set; }
        /// <summary>
        /// Glycan structure string representing the glycan structure and linkage. Example: (N(H(A))(N(H(A))(F)))
        /// </summary>
        public string Struc { get; private set; }
        /// <summary>
        /// Glycan mass, stored as an integer scaled by 1e5.
        /// </summary>
        public int Mass { get; private set; }
        /// <summary>
        /// Type of glycan (N-glycan, O-glycan).
        /// </summary>
        public GlycanType Type;


        /// <summary>
        /// Glycan composition array. Each number corresponds to one type of monosaccharide (order matches Glycan.CharMassDic).
        /// </summary>
        public byte[] Kind { get; private set; }
        /// <summary>
        /// Glycan composition string. Example: H2N2A2F1.
        /// </summary>
        public string Composition
        {
            get
            {
                return Glycan.GetKindString(Kind);
            }
        }
        /// <summary>
        /// List of glycan fragment ions.
        /// </summary>
        public List<GlycanIon> Ions { get; set; }
        /// <summary>
        /// Indicates whether the glycan is a decoy.
        /// </summary>
        public bool Decoy { get; private set; }

        /// <summary>
        /// Set of diagnostic ion masses (B ions) for the glycan, used for glycopeptide identification.
        /// </summary>
        public HashSet<int> GlycanDiagnosticIons
        {
            get
            {
                HashSet<int> diagnosticIons = new HashSet<int>();
                if (Kind[0] >= 1) //if we have Hexose(the number more than one), then we have the corresponding diagonsitic ions as below.
                {
                    diagnosticIons.Add(10902895 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(11503951 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(16306064 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[1] >= 1) // if we have HexNAc(the number more than one), then we have the corresponding diagonsitic ions as below.
                {
                    diagnosticIons.Add(12605550 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(13805550 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(14406607 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(16806607 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(18607663 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(20408720 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[1] >= 1 && Kind[0] >= 1) // if we have HexNAc and Hexose, then we have the corresponding diagonsitic ions as below.
                {
                    diagnosticIons.Add(36614002 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[2] >= 1) //If we have NeuNAc, then we have the corresponding diagonsitic ions as below.
                {
                    diagnosticIons.Add(27409268 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(29210324 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[3] >= 1) //If we have NeuNGc, then we have the corresponding diagonsitic ions as below.
                {
                    diagnosticIons.Add(29008759 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(30809816 - hydrogenAtomMonoisotopicMass);
                }
                // Custom-monosaccharide diagnostic ions: emitted at the value supplied by the user
                // in MonosaccharidesCustom.tsv (no hydrogen-mass offset is applied; users provide
                // observed m/z directly).
                foreach (var kv in _customDiagnosticIonsByIndex)
                {
                    if (kv.Key < Kind.Length && Kind[kv.Key] >= 1)
                    {
                        foreach (int ion in kv.Value)
                        {
                            diagnosticIons.Add(ion);
                        }
                    }
                }
                return diagnosticIons;
            }
        }

        #region Glycan information

        private static readonly int hydrogenAtomMonoisotopicMass =  Convert.ToInt32(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 1E5);


        // Master ordered list of monosaccharide slots. Position in this list IS the Kind[] index.
        // The 11 built-in entries are seeded here; LoadCustomMonosaccharides appends additional
        // entries at indices 11+. Built-in indices (0..10) are referenced by name in places like
        // GlycanDiagnosticIons and OGlycanCompositionCombinationChildIons filter rules, so existing
        // built-in semantics are unchanged by adding customs.
        //
        // H: C6O5H10 Hexose, N: C8O5NH13 HexNAc, A: C11O8NH17 Neu5Ac, G: C11H17NO9 Neu5Gc,
        // F: C6O4H10 Fucose, P: PO3H Phosphate, S: SO3H Sulfo, Y: Na Sodium, C: Acetyl for Neu5Ac,
        // X: C5H10O5 Xylose, K: Kdn
        private static readonly List<(string CanonicalName, char Code, int MassScaled)> _builtInKindEntries =
            new List<(string, char, int)>
        {
            ("Hex",     'H', 16205282),
            ("HexNAc",  'N', 20307937),
            ("NeuAc",   'A', 29109542),
            ("NeuGc",   'G', 30709033),
            ("Fuc",     'F', 14605791),
            ("Phospho", 'P',  7996633),
            ("Sulfo",   'S',  7995681),
            ("Na",      'Y',  2298977),
            ("Ac",      'C',  4201056),
            ("Xylose",  'X', 15005282),
            ("Kdn",     'K', 25006897),
        };

        // Active list (built-ins + any loaded customs). Length defines Kind[] capacity.
        private static List<(string CanonicalName, char Code, int MassScaled)> _kindEntries =
            new List<(string, char, int)>(_builtInKindEntries);

        // Custom monosaccharides' diagnostic ions (int, scaled by 1e5), keyed by Kind[] index.
        // Emitted by GlycanDiagnosticIons whenever the glycan has count >= 1 at the custom index.
        private static Dictionary<int, int[]> _customDiagnosticIonsByIndex = new Dictionary<int, int[]>();

        /// <summary>
        /// Number of monosaccharide slots in the Kind[] array (built-ins + any registered customs).
        /// </summary>
        public static int KindCapacity => _kindEntries.Count;

        /// <summary>
        /// Dictionary mapping monosaccharide character codes to their integer mass (scaled by 1e5).
        /// Rebuilt whenever a custom monosaccharide is registered.
        /// </summary>
        public static Dictionary<char, int> CharMassDic { get; private set; } = BuildCharMassDic();

        /// <summary>
        /// Dictionary mapping monosaccharide names to their character code and index in the Kind array.
        /// Rebuilt whenever a custom monosaccharide is registered. "dHex" is permanently aliased to "Fuc".
        /// </summary>
        public static Dictionary<string, Tuple<char, int>> NameCharDic { get; private set; } = BuildNameCharDic();

        private static Dictionary<char, int> BuildCharMassDic()
        {
            var dic = new Dictionary<char, int>(_kindEntries.Count);
            foreach (var e in _kindEntries)
            {
                dic[e.Code] = e.MassScaled;
            }
            return dic;
        }

        private static Dictionary<string, Tuple<char, int>> BuildNameCharDic()
        {
            var dic = new Dictionary<string, Tuple<char, int>>(_kindEntries.Count + 1);
            for (int i = 0; i < _kindEntries.Count; i++)
            {
                var e = _kindEntries[i];
                dic[e.CanonicalName] = Tuple.Create(e.Code, i);
            }
            // Preserved built-in alias: dHex -> Fuc (same slot, same code).
            if (dic.ContainsKey("Fuc"))
            {
                dic["dHex"] = dic["Fuc"];
            }
            return dic;
        }

        /// <summary>
        /// Register a custom monosaccharide. Called by GlycanDatabase.LoadCustomMonosaccharides
        /// at startup. Appends a new slot at the end of Kind[] (index = KindCapacity-after-add - 1).
        /// </summary>
        /// <param name="name">Unique name, must not collide with a built-in or already-registered name.</param>
        /// <param name="code">Unique single ASCII letter, must not collide with a built-in or already-registered code.</param>
        /// <param name="massScaled">Monoisotopic mass in Da multiplied by 1e5 (e.g. 176.03209 Da -> 17603209).</param>
        /// <param name="diagnosticIonsScaled">Optional diagnostic ion m/z values (scaled by 1e5); null or empty for none.</param>
        public static void RegisterCustomMonosaccharide(string name, char code, int massScaled, int[] diagnosticIonsScaled)
        {
            if (string.IsNullOrWhiteSpace(name))
                throw new ArgumentException("Monosaccharide name must be non-empty.", nameof(name));
            if (NameCharDic.ContainsKey(name))
                throw new ArgumentException($"Monosaccharide name '{name}' already exists (built-in or previously-registered).", nameof(name));
            if (CharMassDic.ContainsKey(code))
                throw new ArgumentException($"Single-char code '{code}' already exists (built-in or previously-registered).", nameof(code));
            if (!((code >= 'A' && code <= 'Z') || (code >= 'a' && code <= 'z')))
                throw new ArgumentException($"Single-char code '{code}' must be an ASCII letter.", nameof(code));

            int newIndex = _kindEntries.Count;
            _kindEntries.Add((name, code, massScaled));
            if (diagnosticIonsScaled != null && diagnosticIonsScaled.Length > 0)
            {
                _customDiagnosticIonsByIndex[newIndex] = diagnosticIonsScaled;
            }
            CharMassDic = BuildCharMassDic();
            NameCharDic = BuildNameCharDic();
        }

        /// <summary>
        /// Removes all registered custom monosaccharides, restoring the built-in 11. Intended for tests.
        /// </summary>
        public static void ResetCustomMonosaccharides()
        {
            _kindEntries = new List<(string, char, int)>(_builtInKindEntries);
            _customDiagnosticIonsByIndex = new Dictionary<int, int[]>();
            CharMassDic = BuildCharMassDic();
            NameCharDic = BuildNameCharDic();
        }

        /// <summary>
        /// Set of common oxonium ion masses (int, scaled by 1e5) used for initial glycopeptide peak filtering.
        /// </summary>
        public readonly static HashSet<int> CommonOxoniumIons = new HashSet<int>
        {13805550, 16806607, 18607663, 20408720, 36614002 };

        /// <summary>
        /// Array of all oxonium ion masses (int, scaled by 1e5) used for building oxonium intensity lists.
        /// </summary>
        public readonly static int[] AllOxoniumIons = new int[]
        {10902895, 11503951, 12605550, 12703952, 13805550, 14406607, 16306064, 16806607, 18607663, 20408720, 27409268, 29008759, 29210324, 30809816, 36614002, 65723544, 67323035};

        /// <summary>
        /// Dictionary mapping N-glycan core ion indices to their monoisotopic mass (double).
        /// </summary>
        public readonly static Dictionary<int, double> TrimannosylCores = new Dictionary<int, double>
        {
            //Each of the mass represent as a N-Glycan core. 
            { 0, 0}, //Y0
            { 83, 83.038194}, //Y*
            { 203, 203.079373}, //Y1
            { 406, 406.158745}, //Y2
            { 568, 568.211568}, //Y3
            { 730, 730.264392}, //Y4
            { 892, 892.317215}, //Y5
            { 349, 349.137281}, //Y2F
            { 552, 552.216654}  //Y3F
        };

        #endregion

        #region Glycan Structure manipulation

        //There are two ways to represent a glycan in string
        //Composition: HexNAc(2)Hex(5)NeuAc(1)NeuGc(1)Fuc(1)Phospho(1)Sulfo(1)Na(1)Ac(1)Xylose(1),
        //Struct(Linkage): (N(H(A))(N(H(A))(F)))

        /// <summary>
        /// Only for Gdb. The method generate a glycan object by reading the glycan structure string from database.
        /// </summary>
        /// <param name="theGlycanStruct"> structrue string ex. (N(H(A))(N(H(A))(F)))</param>
        /// <param name="id"></param>
        /// <param name="isOglycan"></param>
        /// <returns> Glycan Object </returns>
        public static List<Glycan> Struct2Glycan(string theGlycanStruct, int id, bool isOglycan = false)
        {
            Node node = Struct2Node(theGlycanStruct);              // String to tree structure.
            List<Node> nodeIons = GetAllChildrenCombination(node); // Get all possible fragmentation & neutralLoss of a glycan.
            int mass = Glycan.GetMass(theGlycanStruct);            // Get glycan mass.
            byte[] kind = Glycan.GetKind(theGlycanStruct);         // Get glycan composition array, EX. [2, 5, 1, 1, 1, 1, 1, 1, 1, 1].
            List<GlycanIon> glycanIons = new List<GlycanIon>();
            HashSet<double> ionMasses = new HashSet<double>();
            foreach (var aNodeIon in nodeIons)
            {
                var ionMass = Glycan.GetMass(Node2Struct(aNodeIon)); // Get the ionMass
                if (!ionMasses.Contains(ionMass) && ionMass != mass) // Avoid duplicate ions with the same mass. Ex. N(H)N and N(N(H)) have the same ionMass.
                {                                                    // We also avoid the ionMass equals to the glycan mass. Because we won't assume the whole glycan is a fragment ion.
                    ionMasses.Add(ionMass);
                    var ionKind = Glycan.GetKind(Node2Struct(aNodeIon));
                    var lossIonMass = GetIonLossMass(kind, ionKind);
                    GlycanIon glycanIon = new GlycanIon(null, ionMass, ionKind, lossIonMass);
                    glycanIons.Add(glycanIon);
                }
            }
            if (!isOglycan)
            {
                glycanIons.Add(new GlycanIon(null, 8303819, new byte[KindCapacity], mass - 8303819)); //Cross-ring mass
            }
            glycanIons.Add(new GlycanIon(null, 0, kind, mass)); //That is Y0 ion. The whole glycan dropped from the glycopeptide. Like a netural loss.

            List<Glycan> glycans = new List<Glycan>();
            if (isOglycan) //Because we will generate two o-Glycan with different motifs
            {
                GlycanType glycanType = GlycanType.O_glycan;
                Glycan Oglycan_S = new Glycan(theGlycanStruct, mass, kind, glycanIons.OrderBy(p => p.IonMass).ToList(), false, "S", glycanType);
                Oglycan_S.GlyId = id;
                Glycan Oglycan_T = new Glycan(theGlycanStruct, mass, kind, glycanIons.OrderBy(p => p.IonMass).ToList(), false, "T", glycanType);
                Oglycan_T.GlyId = id+1;

                glycans.Add(Oglycan_S);
                glycans.Add(Oglycan_T);

                return glycans;
            }
            else 
            {
                GlycanType glycanType = GlycanType.N_glycan;
                Glycan N_glycan_Nxs = new Glycan(theGlycanStruct, mass, kind, glycanIons.OrderBy(p => p.IonMass).ToList(), false, "Nxs", glycanType);
                N_glycan_Nxs.GlyId = id;
                Glycan N_glycan_Nxt = new Glycan(theGlycanStruct, mass, kind, glycanIons.OrderBy(p => p.IonMass).ToList(), false, "Nxt", glycanType);
                N_glycan_Nxt.GlyId = id+1;

                glycans.Add(N_glycan_Nxs);
                glycans.Add(N_glycan_Nxt);
                return glycans;
            }

        }


        /// <summary>
        /// Convert the glycan structure string to tree format
        /// </summary>
        /// <param name="theGlycanStruct"> linkage inforamtion ex. (N(H))</param>
        /// <returns> glycan tree node ex. Current Nonde = Node(N, 0), left Child = Node(H, 1)</returns>
        public static Node Struct2Node(string theGlycanStruct)
        {
            int level = 0;
            Node curr = new Node(theGlycanStruct[1], level);     // The first character is always '(', so the second character is the root of the tree. In this case of (N(H)), N is the root.
            for (int i = 2; i < theGlycanStruct.Length - 1; i++) // Try to extract the following characters.
            {
                if (theGlycanStruct[i] == '(')                   // Skip the '(' character.
                {
                    continue;
                }
                if (theGlycanStruct[i] == ')')                   // When we meet a ')', we need to go back to the parent node.
                {
                    curr = curr.Father;
                    level--;
                }
                else         // While meeting a character, we need to decide where to put it in the tree. (putting priority: left -> right side -> middle)
                {
                    level++; // Move to the level.(Deeper/Child level)
                    if (curr.LeftChild == null)
                    {
                        curr.LeftChild = new Node(theGlycanStruct[i], level);
                        curr.LeftChild.Father = curr;
                        curr = curr.LeftChild;
                    }
                    else if (curr.RightChild == null)
                    {
                        curr.RightChild = new Node(theGlycanStruct[i], level);
                        curr.RightChild.Father = curr;
                        curr = curr.RightChild;
                    }
                    else if (curr.MiddleChild == null)
                    {
                        curr.MiddleChild = curr.LeftChild;
                        curr.LeftChild = new Node(theGlycanStruct[i], level);
                        curr.LeftChild.Father = curr;
                        curr = curr.LeftChild;
                    }
                }
            }
            return curr; 

        }

        
        /// <summary>
        /// Generate all possible fragments(subset) of a glycan. The fragments are also represented by a Node.
        /// </summary>
        /// <param name="node"></param>
        /// <returns> The all combination of the Glycan fragment. Presented by Node </returns>
        private static List<Node> GetAllChildrenCombination(Node node)
        {
            List<Node> nodes = new List<Node>();
            var curr = node;
            if (curr.LeftChild == null && curr.RightChild == null)
            {
                nodes.Add(curr);
            }
            else
            {
                List<Node> l = GetAllChildrenCombination(curr.LeftChild);
                nodes.Add(new Node(curr.Value));
                if (curr.RightChild != null)
                {
                    List<Node> r = GetAllChildrenCombination(curr.RightChild);
                    if (curr.MiddleChild != null)
                    {
                        List<Node> m = GetAllChildrenCombination(curr.MiddleChild);

                        foreach (var lNode in l)
                        {
                            var c = new Node(curr.Value);
                            c.LeftChild = lNode;
                            nodes.Add(c);
                        }
                        foreach (var rNode in r)
                        {
                            var c = new Node(curr.Value);
                            c.RightChild = rNode;
                            nodes.Add(c);
                        }
                        foreach (var mNode in m)
                        {
                            var c = new Node(curr.Value);
                            c.MiddleChild = mNode;
                            nodes.Add(c);
                        }
                        foreach (var lNode in l)
                        {
                            foreach (var rNode in r)
                            {
                                var c = new Node(curr.Value);
                                c.LeftChild = lNode;
                                c.RightChild = rNode;
                                nodes.Add(c);
                            }
                        }

                        foreach (var rNode in r)
                        {
                            foreach (var mNode in m)
                            {
                                var c = new Node(curr.Value);
                                c.RightChild = rNode;
                                c.MiddleChild = mNode;
                                nodes.Add(c);
                            }
                        }

                        foreach (var lNode in l)
                        {
                            foreach (var mNode in m)
                            {
                                var c = new Node(curr.Value);
                                c.LeftChild = lNode;
                                c.MiddleChild = mNode;
                                nodes.Add(c);
                            }
                        }

                        foreach (var lNode in l)
                        {
                            foreach (var rNode in r)
                            {
                                foreach (var mNode in m)
                                {
                                    var c = new Node(curr.Value);
                                    c.LeftChild = lNode;
                                    c.RightChild = rNode;
                                    c.MiddleChild = mNode;
                                    nodes.Add(c);
                                }
                            }
                        }
                    }
                    else
                    {
                        foreach (var lNode in l)
                        {
                            var c = new Node(curr.Value);
                            c.LeftChild = lNode;
                            nodes.Add(c);
                        }
                        foreach (var rNode in r)
                        {
                            var c = new Node(curr.Value);
                            c.RightChild = rNode;
                            nodes.Add(c);
                        }
                        foreach (var lNode in l)
                        {
                            foreach (var rNode in r)
                            {
                                var c = new Node(curr.Value);
                                c.LeftChild = lNode;
                                c.RightChild = rNode;
                                nodes.Add(c);
                            }
                        }
                    }
                }
                else
                {
                    foreach (var lNode in l)
                    {
                        var c = new Node(curr.Value);
                        c.LeftChild = lNode;
                        nodes.Add(c);
                    }
                }
            }

            return nodes;
        }

        //Node structure to string structure.
        // input: Node(N, 0) -> left Child = Node(H, 1), output: (N(H))
        private static string Node2Struct(Node node)
        {
            string output = "";
            if (node != null)
            {
                output += "(" + node.Value + Node2Struct(node.MiddleChild) + Node2Struct(node.LeftChild) + Node2Struct(node.RightChild) + ")";
            }
            return output;
        }

        /// <summary>
        /// Calculate the mass difference of two glycan kind.
        /// </summary>
        /// <param name="Kind"> Composition of the glycan </param>
        /// <param name="ionKind"> Composition of the glycanIon </param>
        /// <returns> Mass different between the glycan and its glycanIon </returns>
        public static int GetIonLossMass(byte[] Kind, byte[] ionKind)
        {
            byte[] lossKind = new byte[Kind.Length];
            for (int i = 0; i < Kind.Length; i++)
            {
                lossKind[i] = (byte)(Kind[i] - ionKind[i]);
            }
            return Glycan.GetMass(lossKind);
        }

        #endregion

        #region Transfer information
        /// <summary>
        /// Get glycan mass by glycan structure string
        /// </summary>
        /// <param name="structure"> ex.(N(H(A))(N(H(A))(F))) </param>
        /// <returns> The glycan Mass </returns>
        private static int GetMass(string structure)
        {
            int y = 0;
            foreach (var entry in _kindEntries)
            {
                y += entry.MassScaled * structure.Count(p => p == entry.Code);
            }
            return y;
        }

        /// <summary>
        /// Get glycan mass by glycan composition
        /// </summary>
        /// <param name="kind"> [2, 2, 2, 0, 1, 0, 0, 0, 0, 0, 0, ...] </param>
        /// <returns> The glycan mass </returns>
        public static int GetMass(byte[] kind)
        {
            int mass = 0;
            int upper = Math.Min(kind.Length, _kindEntries.Count);
            for (int i = 0; i < upper; i++)
            {
                mass += _kindEntries[i].MassScaled * kind[i];
            }
            return mass;
        }


        /// <summary>
        /// Get glycan composition by the structure string
        /// </summary>
        /// <param name="structure"> structure format : (N(H(A))(N(H(A))(F))) </param>
        /// <returns> The kind array, length == KindCapacity. </returns>
        public static byte[] GetKind(string structure)
        {
            var kind = new byte[_kindEntries.Count];
            for (int i = 0; i < _kindEntries.Count; i++)
            {
                kind[i] = Convert.ToByte(structure.Count(p => p == _kindEntries[i].Code));
            }
            return kind;
        }


        /// <summary>
        /// Get glycan composition text from the glycan kind[]. Slots with zero count are omitted.
        /// </summary>
        /// <param name="Kind"> ex. [2, 2, 2, 0, 1, 0, 0, 0, 0, 0, 0, ...] </param>
        /// <returns> The composition text ex. H2N2A2F1 </returns>
        public static string GetKindString(byte[] Kind)
        {
            var sb = new System.Text.StringBuilder();
            int upper = Math.Min(Kind.Length, _kindEntries.Count);
            for (int i = 0; i < upper; i++)
            {
                if (Kind[i] != 0)
                {
                    sb.Append(_kindEntries[i].Code).Append(Kind[i]);
                }
            }
            return sb.ToString();
        }

        #endregion


        #region Combination or Permutation functions not directly related to glycan, use carefully these function don't deal duplicate elements.


        public static IEnumerable<IEnumerable<T>> GetKCombs<T>(IEnumerable<T> list, int length) where T : IComparable
        {
            if (length == 1) return list.Select(t => new T[] { t });  // Return the list of the single element.
            return GetKCombs(list, length - 1).SelectMany(t => list.Where(o => o.CompareTo(t.Last()) > 0), (t1, t2) => t1.Concat(new T[] { t2 }));
        }

        public static IEnumerable<IEnumerable<T>> GetKCombsWithRept<T>(IEnumerable<T> list, int length) where T : IComparable
        {
            if (length == 1) return list.Select(t => new T[] { t });
            return GetKCombsWithRept(list, length - 1).SelectMany(t => list.Where(o => o.CompareTo(t.Last()) >= 0), (t1, t2) => t1.Concat(new T[] { t2 }));
        }

        public static IEnumerable<IEnumerable<T>> GetPermutations<T>(IEnumerable<T> list, int length)
        {
            if (length == 1)
            {
                return list.Select(t => new T[] { t }); 
            }
            return GetPermutations(list, length - 1).SelectMany(t => list.Where(o => !t.Contains(o)), (t1, t2) => t1.Concat(new T[] { t2 }));
        }

        public static IEnumerable<IEnumerable<T>> GetPermutationsWithRept<T>(IEnumerable<T> list, int length)
        {
            if (length == 1) return list.Select(t => new T[] { t });
            return GetPermutationsWithRept(list, length - 1).SelectMany(t => list, (t1, t2) => t1.Concat(new T[] { t2 }));
        }

        #endregion

        #region Functions are not used now, could be useful in the future.      

        /// <summary>
        /// Test the equality of two glycan objects. Including the glycan mass and the glycan ions should be totally indentical.
        /// </summary>
        /// <param name="glycan1"></param>
        /// <param name="glycan2"></param>
        /// <returns></returns>
        public static bool Equals(Glycan glycan1, Glycan glycan2)
        {
            if (glycan1.Mass == glycan2.Mass)
            {
                if (glycan1.Ions.Count() == glycan2.Ions.Count())
                {

                    for (int i = 0; i < glycan1.Ions.Count(); i++)
                    {
                        if (glycan1.Ions[i].IonMass != glycan2.Ions[i].IonMass)
                        {
                            return false;
                        }
                    }
                    return true;
                }
            }
            return false;
        }

        public static Glycan[] BuildTargetDecoyGlycans(IEnumerable<Glycan> glycans) //Build target-decoy glycans for testing.
        {
            List<Glycan> allGlycans = new List<Glycan>();

            Random random = new Random();
            foreach (var aGlycan in glycans)
            {
                allGlycans.Add(aGlycan);
                List<GlycanIon> glycanIons = new List<GlycanIon>();
                foreach (var ion in aGlycan.Ions)
                {
                    var value = random.Next(100000, 3000000); //Based on pGlyco [1, 30] and GlycoPAT [-50, 50].
                    GlycanIon glycanIon = new GlycanIon(null, ion.IonMass + value, ion.IonKind, ion.LossIonMass - value);
                    glycanIons.Add(glycanIon);
                }
                var DecoyGlycan = new Glycan(aGlycan.Struc, aGlycan.Mass, aGlycan.Kind, glycanIons, true, aGlycan.Target.ToString(), aGlycan.Type);
                DecoyGlycan.GlyId = aGlycan.GlyId;
                allGlycans.Add(DecoyGlycan);
            }
            return allGlycans.OrderBy(p => p.Mass).ToArray();
        }

        #endregion

    }    

}