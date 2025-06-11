using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using EngineLayer.GlycoSearch;
using Proteomics;
using MassSpectrometry;
using Omics.Modifications;
using UsefulProteomicsDatabases.Generated;

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

    public class Glycan :  Modification
    {
        public Glycan(string struc, int mass, byte[] kind, List<GlycanIon> ions, bool decoy, string motif, GlycanType type = GlycanType.O_glycan) 
            : base( _monoisotopicMass: mass / 1E5, _locationRestriction: "Anywhere.")
        {
            Struc = struc;
            Mass = mass;
            Kind = kind;
            Ions = ions;
            Decoy = decoy;

            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            // Generate the neural loss and diagnostic ions for O_glycan.
            if (type == GlycanType.O_glycan && Ions != null)
            {
                List<double> lossMasses = Ions.Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList();
                neutralLosses.Add(DissociationType.HCD, lossMasses);
                neutralLosses.Add(DissociationType.CID, lossMasses);
                neutralLosses.Add(DissociationType.EThcD, lossMasses);
                ModificationType = "O-Glycosylation"; // Set the modification type to N-Glycosylation.
            }

            // Generate the neural loss and diagnostic ions for N_glycan.
            else if (type == GlycanType.N_glycan && Ions != null)
            {
                List<double> lossMasses = Ions.Where(p => p.IonMass < 57000000).Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList(); //570 is a cutoff for glycan ion size 2N1H, which will generate fragment ions. 
                neutralLosses.Add(DissociationType.HCD, lossMasses);
                neutralLosses.Add(DissociationType.CID, lossMasses);
                neutralLosses.Add(DissociationType.EThcD, lossMasses);
                ModificationType = "N-Glycosylation"; // Set the modification type to N-Glycosylation.
                
            }

            Dictionary<DissociationType, List<double>> diagnosticIons = new Dictionary<DissociationType, List<double>>();
            diagnosticIons.Add(DissociationType.HCD, DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.CID, DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.EThcD, DiagnosticIons.Select(p => (double)p / 1E5).ToList());


            ModificationMotif.TryGetMotif(motif, out ModificationMotif finalMotif); //TO DO: only one motif can be write here.
            var id = Glycan.GetKindString(Kind);

            OriginalId = id; // Set the original ID to the glycan kind string, which is unique for each glycan.
            Target = finalMotif; // Set the target motif for the modification.
            NeutralLosses = neutralLosses; // Set the neutral losses for the modification.
            base.DiagnosticIons = diagnosticIons; // Set the diagnostic ions for the modification.
        }

        public Glycan(byte[] kind)
        {
            Kind = kind;
            Mass = GetMass(kind);
        }

        public int GlyId { get; set; }            // Glycan ID, which is the index of glycan in the glycan database.
        public string Struc { get; private set; } // Glycan structure string represented the glycan structure and linkage. Ex. (N(H(A))(N(H(A))(F)))
        public int Mass { get; private set; }

        public enum Type; // GlycanType N_glycan, O_glycan, Undefined;


        public byte[] Kind { get; private set; }  // Glycans are composed of several types of mono suagr. In Kind, each number correspond to one type (corresponded order as Glycan.CharMassDic).
        public string Composition                 // Glycan composition string. Ex. H2N2A2F1.
        {
            get
            {
                return Glycan.GetKindString(Kind);
            }
        }
        public List<GlycanIon> Ions { get; set; }
        public bool Decoy { get; private set; }

        public HashSet<int> DiagnosticIons // B ions (the sugar fragment dropped from the glycopeptide), used for the N-glycan. There are more ions to set...
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
                return diagnosticIons;
            }

        }

        #region Glycan information

        private static readonly int hydrogenAtomMonoisotopicMass =  Convert.ToInt32(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 1E5);


        //Glycan mass dictionary
        //H: C6O5H10 Hexose, N: C8O5NH13 HexNAc, A: C11O8NH17 Neu5Ac, G: C11H17NO9 Neu5Gc, F: C6O4H10 Fucose, 
        //P: PO3H Phosphate, S: SO3H Sulfo, Y: Na Sodium, C:Acetyl for Neu5Ac
        //X: C5H10O5 Xylose
        private readonly static Dictionary<char, int> CharMassDic = new Dictionary<char, int> {
            { 'H', 16205282 },
            { 'N', 20307937 },
            { 'A', 29109542 },
            { 'G', 30709033 },
            { 'F', 14605791 },
            { 'P', 7996633 },
            { 'S', 7995681 },
            { 'Y', 2298977 },
            { 'C',  4201056 },
            { 'X', 15005282 },
            { 'K', 25006897 },
        };

        // The corresponding index for sugar and Kind.
        public readonly static Dictionary<string, Tuple<char, int>> NameCharDic = new Dictionary<string, Tuple<char, int>>
        {
            {"Hex", new Tuple<char, int>('H', 0) },
            {"HexNAc", new Tuple<char, int>('N', 1) },
            {"NeuAc", new Tuple<char, int>('A', 2) },
            {"NeuGc", new Tuple<char, int>('G', 3) },
            {"Fuc",  new Tuple<char, int>('F', 4)},
            {"Phospho", new Tuple<char, int>('P', 5)},
            {"Sulfo", new Tuple<char, int>('S', 6) },
            {"Na", new Tuple<char, int>('Y', 7) },
            {"Ac", new Tuple<char, int>('C', 8) },
            {"Xylose", new Tuple<char, int>('X', 9) },
            {"Kdn", new Tuple<char,int>('K',10)}
        };

        //The same ion as we describe above in the diagnostic ions. That just for the initial filtering for glycopeptide peaks. Not used now.
        public readonly static HashSet<int> CommonOxoniumIons = new HashSet<int> 
        {13805550, 16806607, 18607663, 20408720, 36614002 };

        //The same ion as we describe above in the diagnostic ions. Used for building the oxoniumIntensity list.
        public readonly static int[] AllOxoniumIons = new int[] 
        {10902895, 11503951, 12605550, 12703952, 13805550, 14406607, 16306064, 16806607, 18607663, 20408720, 27409268, 29008759, 29210324, 30809816, 36614002, 65723544, 67323035};

        //TrimannosylCore. Only useful for N-Glyco peptides.
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
        public static Glycan Struct2Glycan(string theGlycanStruct, int id, bool isOglycan = false)
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
                glycanIons.Add(new GlycanIon(null, 8303819, new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, mass - 8303819)); //Cross-ring mass
            }
            glycanIons.Add(new GlycanIon(null, 0, kind, mass)); //That is Y0 ion. The whole glycan dropped from the glycopeptide. Like a netural loss.

            Glycan glycan = new Glycan(theGlycanStruct, mass, kind, glycanIons.OrderBy(p => p.IonMass).ToList(), false);
            glycan.GlyId = id;
            return glycan;
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
            int y = CharMassDic['H'] * structure.Count(p => p == 'H') +
                CharMassDic['N'] * structure.Count(p => p == 'N') +
                CharMassDic['A'] * structure.Count(p => p == 'A') +
                CharMassDic['G'] * structure.Count(p => p == 'G') +
                CharMassDic['F'] * structure.Count(p => p == 'F') +
                CharMassDic['P'] * structure.Count(p => p == 'P') +
                CharMassDic['S'] * structure.Count(p => p == 'S') +
                CharMassDic['Y'] * structure.Count(p => p == 'Y') +
                CharMassDic['C'] * structure.Count(p => p == 'C') +
                CharMassDic['X'] * structure.Count(p => p == 'X') +
                CharMassDic['K'] * structure.Count(p => p == 'K')
                ;
            return y;
        }

        /// <summary>
        /// Get glycan mass by glycan composition
        /// </summary>
        /// <param name="kind"> [2, 2, 2, 0, 1, 0, 0, 0, 0, 0] </param>
        /// <returns> The glycan mass </returns>
        public static int GetMass(byte[] kind) 
        {
            int mass = CharMassDic['H'] * kind[0] +
            CharMassDic['N'] * kind[1] +
            CharMassDic['A'] * kind[2] +
            CharMassDic['G'] * kind[3] +
            CharMassDic['F'] * kind[4] +
            CharMassDic['P'] * kind[5] +
            CharMassDic['S'] * kind[6] +
            CharMassDic['Y'] * kind[7] +
            CharMassDic['C'] * kind[8] +
            CharMassDic['X'] * kind[9] +
            CharMassDic['K'] * kind[10]
            ;

            return mass;
        }


        /// <summary>
        /// Get glycan composition by the structure string
        /// </summary>
        /// <param name="structure"> structure format : (N(H(A))(N(H(A))(F))) </param>
        /// <returns> The kind List ex [2, 2, 2, 0, 1, 0, 0, 0, 0, 0].</returns>
        public static byte[] GetKind(string structure) 
        {
            var kind = new byte[] 
            { Convert.ToByte(structure.Count(p => p == 'H')),
                Convert.ToByte(structure.Count(p => p == 'N')),
                Convert.ToByte(structure.Count(p => p == 'A')),
                Convert.ToByte(structure.Count(p => p == 'G')),
                Convert.ToByte(structure.Count(p => p == 'F')),
                Convert.ToByte(structure.Count(p => p == 'P')),
                Convert.ToByte(structure.Count(p => p == 'S')),
                Convert.ToByte(structure.Count(p => p == 'Y')),
                Convert.ToByte(structure.Count(p => p == 'C')),
                Convert.ToByte(structure.Count(p => p == 'X')),
                Convert.ToByte(structure.Count(p => p == 'K'))
            };
            return kind;
        }


        /// <summary>
        /// Get glycan composition text from the glycan kind[].
        /// </summary>
        /// <param name="Kind"> ex. [2, 2, 2, 0, 1, 0, 0, 0, 0, 0] </param>
        /// <returns> The composition text ex. H2N2A2F1 </returns>
        public static string GetKindString(byte[] Kind)
        {
            string H = Kind[0]==0 ? "" : "H" + Kind[0].ToString();
            string N = Kind[1] == 0 ? "" : "N" + Kind[1].ToString();
            string A = Kind[2] == 0 ? "" : "A" + Kind[2].ToString();
            string G = Kind[3] == 0 ? "" : "G" + Kind[3].ToString();
            string F = Kind[4] == 0 ? "" : "F" + Kind[4].ToString();
            string P = Kind[5] == 0 ? "" : "P" + Kind[5].ToString();
            string S = Kind[6] == 0 ? "" : "S" + Kind[6].ToString();
            string Y = Kind[7] == 0 ? "" : "Y" + Kind[7].ToString();
            string C = Kind[8] == 0 ? "" : "C" + Kind[8].ToString();
            string X = Kind[9] == 0 ? "" : "X" + Kind[9].ToString();
            string K = Kind[10] == 0 ? "" : "K" + Kind[10].ToString();
            string kindString = H + N + A + G + F + P + S + Y + C + X + K;
            return kindString;
        }

        #endregion

        //TO THINK: Is it reasonable to transfer Glycan to Modification the first time Glycan is read in? Which could save time.
        //Use glycan index and modification index to reduce space.

        /// <summary>
        /// Input the N-glycan object, and transfer it to the modification object.
        /// </summary>
        /// <param name="glycan"></param>
        /// <returns></returns>
        public static Modification NGlycanToModification(Glycan glycan)
        {
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            if (glycan.Ions!=null)
            {
                List<double> lossMasses = glycan.Ions.Where(p => p.IonMass < 57000000).Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList(); //570 is a cutoff for glycan ion size 2N1H, which will generate fragment ions. 
                neutralLosses.Add(DissociationType.HCD, lossMasses);
                neutralLosses.Add(DissociationType.CID, lossMasses);
                neutralLosses.Add(DissociationType.EThcD, lossMasses);
            }

            Dictionary<DissociationType, List<double>> diagnosticIons = new Dictionary<DissociationType, List<double>>();
            diagnosticIons.Add(DissociationType.HCD, glycan.DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.CID, glycan.DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.EThcD, glycan.DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            ModificationMotif.TryGetMotif("N", out ModificationMotif finalMotif); //TO DO: only one motif can be write here.
            var id = Glycan.GetKindString(glycan.Kind);
            Modification modification = new Modification(
                _originalId: id,
                _modificationType: "N-Glycosylation",
                _monoisotopicMass: (double)glycan.Mass / 1E5,
                _locationRestriction: "Anywhere.",
                _target: finalMotif,
                _neutralLosses: neutralLosses,
                _diagnosticIons: diagnosticIons
            );
            return modification;
        }

        /// <summary>
        /// Input the O-glycan object, and transfer it to the modification object.
        /// </summary>
        /// <param name="glycan"></param>
        /// <returns> The modification object </returns>
        public static Modification OGlycanToModification(Glycan glycan) //try to transfer the glycan object to modification object.
        {
            //TO THINK: what the neutralLoss for O-Glyco?
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();

            if (glycan.Ions!=null)
            {
                List<double> lossMasses = glycan.Ions.Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList();
                neutralLosses.Add(DissociationType.HCD, lossMasses);
                neutralLosses.Add(DissociationType.CID, lossMasses);
                neutralLosses.Add(DissociationType.EThcD, lossMasses);
            }

            Dictionary<DissociationType, List<double>> diagnosticIons = new Dictionary<DissociationType, List<double>>();
            diagnosticIons.Add(DissociationType.HCD, glycan.DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.CID, glycan.DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            diagnosticIons.Add(DissociationType.EThcD, glycan.DiagnosticIons.Select(p => (double)p / 1E5).ToList());
            ModificationMotif.TryGetMotif("X", out ModificationMotif finalMotif); //TO DO: only one motif can be write here.

            var id = Glycan.GetKindString(glycan.Kind);
            Modification modification = new Modification(
                _originalId: id,
                _modificationType: "O-Glycosylation",
                _monoisotopicMass: (double)glycan.Mass / 1E5,
                _locationRestriction: "Anywhere.",
                _target: finalMotif,
                _neutralLosses: neutralLosses,
                _diagnosticIons: diagnosticIons
            );
            return modification;
        }

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
                var aDecoyGlycan = new Glycan(aGlycan.Struc, aGlycan.Mass, aGlycan.Kind, glycanIons, true);
                aDecoyGlycan.GlyId = aGlycan.GlyId;
                allGlycans.Add(aDecoyGlycan);
            }
            return allGlycans.OrderBy(p => p.Mass).ToArray();
        }

        #endregion

    }    

}