using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;
using System.Dynamic;
using System.Reflection.Metadata;
using System.Text;
namespace EngineLayer
{
    public class Monosaccharide
    {
        public Monosaccharide(byte id, string name, char symbol, int mass)
        {
            Id = id;
            Name = name;
            Symbol = symbol;
            Mass = mass;
        }
        public byte Id { get; }
        public string Name { get; }
        public char Symbol { get; }
        public int Mass { get; }

        public static Monosaccharide[] LoadMonosaccharide(string filePath)
        {
            List<Monosaccharide> sugars = new List<Monosaccharide>();
            using (StreamReader lines = new StreamReader(filePath))
            {
                bool firstLine = true;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine();

                    if (firstLine)
                    {
                        //Skip the first line
                        firstLine = false;
                        continue;
                    }

                    var xs = line.Split('\t');

                    Monosaccharide mono = new Monosaccharide(byte.Parse(xs[0].Trim()), xs[1].Trim(), char.Parse(xs[2].Trim()), int.Parse(xs[3].Trim()));
                    sugars.Add(mono);
                }
            }
         
            return sugars.ToArray();
        }
        public static Dictionary<char, int> GetCharMassDic(Monosaccharide[] monosaccharides)
        {
            Dictionary<char, int> sugarMass = new Dictionary<char, int>();
            foreach (var mono in monosaccharides)
            {
                sugarMass.Add(mono.Symbol, mono.Mass);
            }
            return sugarMass;
        }
        public static Dictionary<string, byte> GetNameIdDic(Monosaccharide[] monosaccharides)
        {
            Dictionary<string, byte> nameId = new Dictionary<string, byte>();
            foreach (var mono in monosaccharides)
            {
                nameId.Add(mono.Name, mono.Id);
            }
            return nameId;
        }
    }
    public class GlycanIon
    {
        public GlycanIon(string ionStruct, int ionMass, byte[] ionKind, int lossIonMass)
        {
            IonStruct = ionStruct;
            IonMass = ionMass;
            IonKind = ionKind;
            LossIonMass = lossIonMass;
        }
        public string IonStruct { get; set; }
        public int IonMass { get; set; }
        public int LossIonMass { get; set; }//Glycan.Mass - IonMass
        public byte[] IonKind { get; set; }
    }

    public class Glycan
    {
        public static readonly int CrossRingMass = 8303819;

        public Glycan(string struc, int mass, byte[] kind, List<GlycanIon> ions, bool decoy)
        {
            Struc = struc;
            Mass = mass;
            Kind = kind;
            Ions = ions;
            Decoy = decoy;
        }

        public Glycan(byte[] kind)
        {
            Kind = kind;
            Mass = GetMass(kind);
        }

        public int GlyId { get; set; }
        public string Struc { get; private set; }
        public int Mass { get; private set; }

        //Glycans are composed of several different types of mono saccharides. In Kind, each number correspond to one type of mono saccharide in the same order as Glycan.CharMassDic. 
        public byte[] Kind { get; private set; }
        public string Composition
        {
            get
            {
                return Glycan.GetKindString(Kind);
            }
        }
        public List<GlycanIon> Ions { get; set; }
        public bool Decoy { get; private set; }

        public HashSet<int> DiagnosticIons
        {
            get
            {
                HashSet<int> diagnosticIons = new HashSet<int>();
                if (Kind[0] >= 1)
                {
                    diagnosticIons.Add(10902895 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(11503951 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(16306064 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[1] >= 1)
                {
                    diagnosticIons.Add(12605550 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(13805550 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(14406607 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(16806607 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(18607663 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(20408720 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[1] >= 1 && Kind[0] >= 1)
                {
                    diagnosticIons.Add(36614002 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[2] >= 1)
                {
                    diagnosticIons.Add(27409268 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(29210324 - hydrogenAtomMonoisotopicMass);
                }
                if (Kind[3] >= 1)
                {
                    diagnosticIons.Add(29008759 - hydrogenAtomMonoisotopicMass);
                    diagnosticIons.Add(30809816 - hydrogenAtomMonoisotopicMass);
                }
                return diagnosticIons;
            }

        }

        #region Glycan information

        private static readonly int hydrogenAtomMonoisotopicMass = Convert.ToInt32(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 1E5);

        public static int SugarLength
        {
            get
            {
                return GlobalVariables.Monosaccharides.Length;
            }
        }

        //Glycan mass dictionary
        //H: C6O5H10 Hexose, N: C8O5NH13 HexNAc, A: C11O8NH17 Neu5Ac, G: C11H17NO9 Neu5Gc, F: C6O4H10 Fucose, 
        //P: PO3H Phosphate, S: SO3H Sulfo, Y: Na Sodium, C:Acetyl for Neu5Ac
        //X: C5H10O5 Xylose
        //Related code include: GetMass, GetKind, GetKindString, GlycanBox constructor, search byte[].
        private static Dictionary<char, int> CharMassDic
        {
            get
            {
                return Monosaccharide.GetCharMassDic(GlobalVariables.Monosaccharides);
            }
        }
      
        public static Dictionary<string, byte> NameIdDic
        {
            get
            {
                return Monosaccharide.GetNameIdDic(GlobalVariables.Monosaccharides);
            }
        }

        public readonly static HashSet<int> CommonOxoniumIons = new HashSet<int>
        {13805550, 16806607, 18607663, 20408720, 36614002 };

        public readonly static int[] AllOxoniumIons = new int[]
        {10902895, 11503951, 12605550, 12703952, 13805550, 14406607, 16306064, 16806607, 18607663, 20408720, 27409268, 29008759, 29210324, 30809816, 36614002, 65723544, 67323035};

        //TrimannosylCore is only useful for N-Glyco peptides.
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

        //There are two ways to represent a glycan in string, one only combination, the other structure.
        //The method generate a glycan by read in a glycan structure string from database.
        public static Glycan Struct2Glycan(string theGlycanStruct, int id, bool ToGenerateIons = true, bool isOglycan = false)
        {
            if (!ToGenerateIons)
            {
                var _kind = GetKind(theGlycanStruct);

                return new Glycan(_kind);
            }

            Node node = Struct2Node(theGlycanStruct);
            List<Node> nodeIons = GetAllChildrenCombination(node);
            int mass = Glycan.GetMass(theGlycanStruct);
            byte[] kind = Glycan.GetKind(theGlycanStruct);
            List<GlycanIon> glycanIons = new List<GlycanIon>();
            HashSet<double> ionMasses = new HashSet<double>();
            foreach (var aNodeIon in nodeIons)
            {
                var ionMass = Glycan.GetMass(Node2Struct(aNodeIon));
                if (!ionMasses.Contains(ionMass) && ionMass != mass)
                {
                    ionMasses.Add(ionMass);
                    var ionKind = Glycan.GetKind(Node2Struct(aNodeIon));
                    var lossIonMass = GetIonLossMass(kind, ionKind);
                    GlycanIon glycanIon = new GlycanIon(null, ionMass, ionKind, lossIonMass);
                    glycanIons.Add(glycanIon);
                }
            }
            if (!isOglycan)
            {
                byte[] aIonKind = new byte[SugarLength];
                glycanIons.Add(new GlycanIon(null, CrossRingMass, aIonKind, mass - CrossRingMass)); //Cross-ring mass
            }
            byte[] zeroKind = new byte[SugarLength];
            glycanIons.Add(new GlycanIon(null, 0, zeroKind, mass));

            Glycan glycan = new Glycan(theGlycanStruct, mass, kind, glycanIons.OrderBy(p => p.IonMass).ToList(), false);
            glycan.GlyId = id;
            return glycan;
        }

        //Glycan are represented in tree structures composed of Node. The function here is to transfer a string into connected Node.
        public static Node Struct2Node(string theGlycanStruct)
        {
            int level = 0;
            Node curr = new Node(theGlycanStruct[1], level);
            for (int i = 2; i < theGlycanStruct.Length - 1; i++)
            {
                if (theGlycanStruct[i] == '(')
                {
                    continue;
                }
                if (theGlycanStruct[i] == ')')
                {
                    curr = curr.Father;
                    level--;
                }
                else
                {
                    level++;
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

        //The function is to generate all possible fragmentation/neutral loss of a glycan, which is a subset of glycan. 
        //Node is tree structured glycan. subset of glycans are also represented by Node.
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
        private static string Node2Struct(Node node)
        {
            string output = "";
            if (node != null)
            {
                output += "(" + node.Value + Node2Struct(node.MiddleChild) + Node2Struct(node.LeftChild) + Node2Struct(node.RightChild) + ")";
            }
            return output;
        }

        //kind are compositions of glycan. The function here is to generate mass difference of two glycan.
        public static int GetIonLossMass(byte[] Kind, byte[] ionKind)
        {
            byte[] lossKind = new byte[SugarLength];
            for (int i = 0; i < Kind.Length; i++)
            {
                lossKind[i] = (byte)(Kind[i] - ionKind[i]);
            }
            return Glycan.GetMass(lossKind);
        }

        #endregion

        #region Transfer information

        public static int GetMass(string structure)
        {
            int y = 0;
            for (int i = 0; i < SugarLength; i++)
            {
                y += CharMassDic[GlobalVariables.Monosaccharides[i].Symbol] * structure.Count(p => p == GlobalVariables.Monosaccharides[i].Symbol);
            }
            return y;
        }

        public static int GetMass(byte[] kind)
        {
            int mass = 0;

            for (int i = 0; i < SugarLength; i++)
            {  
                mass += CharMassDic[GlobalVariables.Monosaccharides[i].Symbol] * kind[i];
            }
            return mass;
        }

        public static byte[] GetKind(string structure)
        {
            byte[] kind = new byte[SugarLength];

            for (int i = 0; i < SugarLength; i++)
            {
                kind[i] = Convert.ToByte(structure.Count(p => p == GlobalVariables.Monosaccharides[i].Symbol));
            }

            return kind;
        }

        public static string GetKindString(byte[] Kind)
        {
            StringBuilder kindString = new StringBuilder();
            for (int i = 0; i < SugarLength; i++)
            {
                kindString.Append(Kind[i] == 0 ? "" : GlobalVariables.Monosaccharides[i].Symbol + Kind[i].ToString());
            }

            return kindString.ToString();
        }

        public static double GetNGlycanLocalMass(int countOfNGlycan)
        {
            byte[] kind = new byte[SugarLength];
            kind[0] = (byte)countOfNGlycan;

            int mass = GetMass(kind);
            return (double)mass / 1E5;
        }

        #endregion

        //TO THINK: Is it reasonable to transfer Glycan to Modification the first time Glycan is read in? Which could save time.
        //Use glycan index and modification index to reduce space.
        public static Modification NGlycanToModification(Glycan glycan)
        {
            Dictionary<DissociationType, List<double>> neutralLosses = new Dictionary<DissociationType, List<double>>();
            if (glycan.Ions!=null)
            {
                List<double> lossMasses = glycan.Ions.Where(p => p.IonMass == CrossRingMass || p.IonMass == 20307937).Select(p => (double)p.LossIonMass / 1E5).OrderBy(p => p).ToList(); //570 is a cutoff for glycan ion size 2N1H, which will generate fragment ions. 
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
                _featureType: "Nxs/t",
                _neutralLosses: neutralLosses,
                _diagnosticIons: diagnosticIons
            );
            return modification;
        }

        //This is a simple version of _NGlycanToModification for GlobalVariables read N-glycan
        public static Modification _NGlycanToModification(Glycan glycan)
        {
            ModificationMotif.TryGetMotif("N", out ModificationMotif finalMotif); //TO DO: only one motif can be write here.
            var id = Glycan.GetKindString(glycan.Kind);
            Modification modification = new Modification(
                _originalId: id,
                _modificationType: "N-Glycosylation",
                _monoisotopicMass: (double)glycan.Mass / 1E5,
                _locationRestriction: "Anywhere.",
                _target: finalMotif,
                _featureType: "Nxs/t"
            );
            return modification;
        }

        public static Modification OGlycanToModification(Glycan glycan)
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
                _featureType: "S/T",
                _neutralLosses: neutralLosses,
                _diagnosticIons: diagnosticIons
            );
            return modification;
        }

        #region Combination or Permutation functions not directly related to glycan, use carefully these function don't deal duplicate elements.

        public static IEnumerable<IEnumerable<T>> GetKCombs<T>(IEnumerable<T> list, int length) where T : IComparable
        {
            if (length == 1) return list.Select(t => new T[] { t });
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

        public static Glycan[] BuildTargetDecoyGlycans(IEnumerable<Glycan> glycans)
        {
            List<Glycan> allGlycans = new List<Glycan>();

            Random random = new Random();
            foreach (var aGlycan in glycans)
            {
                //Add original glycan
                allGlycans.Add(aGlycan);

                //Generate decoy glycan by shiftting glycanIons
                List<GlycanIon> glycanIons = new List<GlycanIon>();
                foreach (var ion in aGlycan.Ions)
                {
                    var value = random.Next(100000, 3000000); //Based on pGlyco [1, 30] and GlycoPAT [-50, 50].
                    GlycanIon glycanIon = new GlycanIon(null, ion.IonMass + value, ion.IonKind, ion.LossIonMass - value);
                    glycanIons.Add(glycanIon);
                }
                var aDecoyGlycan = new Glycan(aGlycan.Struc, aGlycan.Mass, aGlycan.Kind, glycanIons, true);

                ////Generate decoy glycan by shiftting N-glycan mass
                //var aDecoyGlycan = new Glycan(aGlycan.Struc, aGlycan.Mass + 2000000, aGlycan.Kind, aGlycan.Ions, true);

                aDecoyGlycan.GlyId = aGlycan.GlyId;
                allGlycans.Add(aDecoyGlycan);
            }
            return allGlycans.OrderBy(p => p.Mass).ToArray();
        }

        #endregion

    }    

}