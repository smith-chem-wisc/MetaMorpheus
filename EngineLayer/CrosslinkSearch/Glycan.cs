using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;

namespace EngineLayer
{
    public class Glycan
    {
        private static readonly int hydrogenAtomMonoisotopicMass =  Convert.ToInt32(PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 1E5);

        //H: C6O5H10, N: C8O5NH13, A: C11O8NH17, G: C11H17NO9, F: C6O4H10, X: C5H10O5, K: C9H16O9, P: PO3H, S: SO3H, R: C6H10O7
        private static Dictionary<char, int> CharMassDic = new Dictionary<char, int>() { { 'H', 16205282 }, { 'N', 20307937 }, { 'A', 29109542 }, { 'G', 30709033 }, { 'F', 14605791 }, {'X', 15005282 }, {'K', 26807943 }, {'P', 7996633 }, {'S', 8096464 }, {'R', 19404265 } };

        public static HashSet<int> oxoniumIons = new HashSet<int>()
        {13805550, 16806607, 18607663, 20408720, 36614002 };

        public static int[] allOxoniumIons = new int[]
        {10902895, 11503951, 12605550, 12703952, 13805550, 14406607, 16306064, 16806607, 18607663, 20408720, 27409268, 29008759, 29210324, 30809816, 36614002, 65723544, 67323035};

        public static Dictionary<int, double> TrimannosylCores = new Dictionary<int, double>()
        {
            //HashSet {83038194, 20307937, 40615875, 56821157, 73026439, 89231722, 34913728, 55221665};
            { 0, 0},
            { 83, 83.038194},
            { 203, 203.079373},
            { 406, 406.158745},
            { 568, 568.211568},
            { 730, 730.264392},
            { 892, 892.317215},
            { 349, 349.137281},
            { 552, 552.216654}
            //{ "Y0", 0},
            //{ "Y*", 83.038194},
            //{ "Y1", 203.079373},
            //{ "Y2", 406.158745},
            //{ "Y3", 568.211568},
            //{ "Y4", 730.264392},
            //{ "Y5", 892.317215},
            //{ "Y2F", 349.137281},
            //{ "Y3F", 552.216654}
        };

        public Glycan(string struc, int mass, byte[] kind, List<GlycanIon> ions, bool decoy)
        {
            Struc = struc;
            Mass = mass;
            Kind = kind;
            Ions = ions;
            Decoy = decoy;
        }
        
        public int GlyId { get; set; }
        public int GlyType { get; set; }
        public string Struc { get; set; }
        public int Mass { get; set; }
        public byte[] Kind { get; set; }
        public List<GlycanIon> Ions { get; set; }
        public bool Decoy { get; set; }
        
        public static Node Struct2Node(string theGlycanStruct)
        {
            int level = 0;
            Node curr = new Node(theGlycanStruct[1], level);
            for (int i = 2; i < theGlycanStruct.Length - 1; i++)
            {

                if (theGlycanStruct[i] != null)
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
                        else if(curr.RightChild == null)
                        {
                            curr.RightChild = new Node(theGlycanStruct[i], level);
                            curr.RightChild.Father = curr;
                            curr = curr.RightChild;
                        }
                        else if(curr.MiddleChild == null)
                        {
                            curr.MiddleChild = curr.LeftChild;
                            curr.LeftChild = new Node(theGlycanStruct[i], level);
                            curr.LeftChild.Father = curr;
                            curr = curr.LeftChild;
                        }
                    }

                }
            }
            return curr;
        }    

        private static string Node2Struct(Node node)
        {
            string output = "";
            if (node != null)
            {
                output += "(" + node.Value + Node2Struct(node.MiddleChild) + Node2Struct(node.LeftChild) + Node2Struct(node.RightChild) + ")";
            }
            return output;
        }

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

        private static int GetMass(string structure)
        {
            int y = CharMassDic['H'] * structure.Count(p => p == 'H') +
                CharMassDic['N'] * structure.Count(p => p == 'N') +
                CharMassDic['A'] * structure.Count(p => p == 'A') +
                CharMassDic['G'] * structure.Count(p => p == 'G') +
                CharMassDic['F'] * structure.Count(p => p == 'F');
            return y;
        }

        public static byte[] GetKind(string structure)
        {
            byte[] kind = new byte[] { Convert.ToByte(structure.Count(p => p == 'H')), Convert.ToByte(structure.Count(p => p == 'N')), Convert.ToByte(structure.Count(p => p == 'A')) , Convert.ToByte(structure.Count(p => p == 'G')) , Convert.ToByte(structure.Count(p => p == 'F')) };
            return kind;
        }

        public static int GetIonLossMass(byte[] Kind, byte[] ionKind)
        {
            byte[] lossKind = new byte[Kind.Length];
            for (int i = 0; i < Kind.Length; i++)
            {
                lossKind[i] = (byte)(Kind[i] - ionKind[i]);
            }
            return GetMass(lossKind);
        }

        //Find glycans from structured glycan database
        public static List<Glycan> GetAllIonMassFromKind(byte[] kind, Dictionary<string, List<Glycan>> groupedGlycans)
        {
            var kindKey = GetKindString(kind);
            List<Glycan> glycans = new List<Glycan>();

            groupedGlycans.TryGetValue(kindKey, out glycans);

            if (glycans==null)
            {
                //if not in the structured glycan database, find a smaller one.
                bool notFound = true;
                while (notFound)
                {
                    var childKinds = BuildChildKindKey(kind);
                    foreach (var child in childKinds)
                    {
                        var key = GetKindString(child);
                        if (groupedGlycans.TryGetValue(key, out glycans))
                        {
                            notFound = false;
                            break;
                        }
                    }

                    if (notFound == true)
                    {
                        glycans = GetAllIonMassFromKind(childKinds[0], groupedGlycans);
                        notFound = false;
                    }
                }
            }

            return glycans;
        }

        private static List<byte[]> BuildChildKindKey(byte[] kind)
        {
            List<byte[]> childKinds = new List<byte[]>();
            for (int i = kind.Length -1; i >=0; i--)
            {             
                if (kind[i] >= 1)
                {
                    var childKind = new byte[kind.Length];
                    Array.Copy(kind, childKind, kind.Length);
                    childKind[i]--;
                    childKinds.Add(childKind);
                }           
            }
            return childKinds;
        }

        public static IEnumerable<Glycan> LoadKindGlycan(string filePath, IEnumerable<Glycan> NGlycans)
        {
            var groupedGlycans = NGlycans.GroupBy(p => GetKindString(p.Kind)).ToDictionary(p => p.Key, p => p.ToList());

            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine();

                    byte[] kind = new byte[10] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  };
                    var x = line.Split('(', ')');
                    int i = 0;
                    while (i < x.Length - 1)
                    {
                        switch (x[i])
                        {
                            case "Hex":
                                kind[0] = byte.Parse(x[i + 1]);
                                break;
                            case "HexNAc":
                                kind[1] = byte.Parse(x[i + 1]);
                                break;
                            case "NeuAc":                            
                                kind[2] = byte.Parse(x[i + 1]);
                                break;
                            case "NeuGc":
                                kind[3] = byte.Parse(x[i + 1]);
                                break;
                            case "Fuc":
                                kind[4] = byte.Parse(x[i + 1]);
                                break;
                            case "Xyl":
                                kind[5] = byte.Parse(x[i + 1]);
                                break;
                            case "KND":
                                kind[6] = byte.Parse(x[i + 1]);
                                break;
                            case "Phosphate":
                                kind[7] = byte.Parse(x[i + 1]);
                                break;
                            case "Sulfate":
                                kind[8] = byte.Parse(x[i + 1]);
                                break;
                            case "HexA":
                                kind[9] = byte.Parse(x[i + 1]);
                                break;
                            default:
                                break;
                        }
                        i = i + 2;
                    }                  
                    var mass = GetMass(kind);
                    
                    var glycans = GetAllIonMassFromKind(kind, groupedGlycans);

                    var glycan = new Glycan(glycans.First().Struc, mass, kind, glycans.First().Ions, true);
                    glycan.GlyId = id++; 
                    yield return glycan; 
                }
            }
        }

        public static bool DistingushGlycans(Glycan glycan1, Glycan glycan2)
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

        public static int GetMass(byte[] kind)
        {
            int mass = 0;
            if (kind.Length <= 5)
            {
                mass = CharMassDic['H'] * kind[0] +
                CharMassDic['N'] * kind[1] +
                CharMassDic['A'] * kind[2] +
                CharMassDic['G'] * kind[3] +
                CharMassDic['F'] * kind[4];
            }
            else
            {
                mass = CharMassDic['H'] * kind[0] +
                CharMassDic['N'] * kind[1] +
                CharMassDic['A'] * kind[2] +
                CharMassDic['G'] * kind[3] +
                CharMassDic['F'] * kind[4] +
                CharMassDic['X'] * kind[5] +
                CharMassDic['K'] * kind[6] +
                CharMassDic['P'] * kind[7] +
                CharMassDic['S'] * kind[8] +
                CharMassDic['R'] * kind[9];
            }
            
            return mass;
        }

        public HashSet<int> GetDiagnosticIons()
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

        public static Glycan Struct2Glycan(string theGlycanStruct, int id)
        {
            Node node = Struct2Node(theGlycanStruct);
            List<Node> nodeIons = GetAllChildrenCombination(node);
            int mass = GetMass(theGlycanStruct);
            byte[] kind = GetKind(theGlycanStruct);
            List<GlycanIon> glycanIons = new List<GlycanIon>();
            HashSet<double> ionMasses = new HashSet<double>();
            foreach (var aNodeIon in nodeIons)
            {
                var ionMass = GetMass(Node2Struct(aNodeIon));
                if (!ionMasses.Contains(ionMass) && ionMass != mass)
                {
                    ionMasses.Add(ionMass);
                    var ionKind = GetKind(Node2Struct(aNodeIon));
                    var lossIonMass = GetIonLossMass(kind, ionKind);
                    GlycanIon glycanIon = new GlycanIon(null, ionMass, ionKind, lossIonMass);
                    glycanIons.Add(glycanIon);
                }
            }
            glycanIons.Add(new GlycanIon(null, 8303819, new byte[] { 0, 0, 0, 0, 0 }, mass - 8303819)); //Cross-ring mass
            glycanIons = glycanIons.OrderBy(p => p.IonMass).ToList();
            //glycanIons.RemoveAt(glycanIons.Count - 1);

            Glycan glycan = new Glycan(theGlycanStruct, mass, kind, glycanIons, false);
            glycan.GlyId = id;
            return glycan;
        }

        public static IEnumerable<Glycan> LoadGlycan(string filePath)
        {
            using (StreamReader glycans = new StreamReader(filePath))
            {
                int id = 1;
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();
                    yield return Struct2Glycan(line, id++);
                }
            }
        }

        public static Glycan[] BuildTargetDecoyGlycans(IEnumerable<Glycan> glycans)
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
            return allGlycans.OrderBy(p=>p.Mass).ToArray();
        }

        public static string GetKindString(byte[] Kind)
        {
            string H =  "H" + Kind[0].ToString();
            string N =  "N" + Kind[1].ToString();
            string A =  "A" + Kind[2].ToString();
            string G =  "G" + Kind[3].ToString();
            string F =  "F" + Kind[4].ToString();
            string kindString = H + N + A + G + F;
            return kindString;
        }

        public static string GetKindString(string structure)
        {
            string H =  "H" + structure.Count(p => p == 'H').ToString();
            string N =  "N" + structure.Count(p => p == 'N').ToString();
            string A =  "A" + structure.Count(p => p == 'A').ToString();
            string G =  "G" + structure.Count(p => p == 'G').ToString();
            string F =  "F" + structure.Count(p => p == 'F').ToString();
            string kindString = H + N + A + G + F;
            return kindString;
        }

        public static IEnumerable<GlycanBox> BuildGlycanBoxes(IEnumerable<Glycan> glycans, int maxNum)
        {

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in GetKCombsWithRept(Enumerable.Range(0, glycans.Count()), i))
                {
                    GlycanBox glycanBox = new GlycanBox();
                    foreach (var id in idCombine)
                    {
                        glycanBox.glycans.Add(glycans.ElementAt(id));

                    }
                    yield return glycanBox;
                }
            }

        }

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
            if (length == 1) return list.Select(t => new T[] { t });
            return GetPermutations(list, length - 1).SelectMany(t => list.Where(o => !t.Contains(o)), (t1, t2) => t1.Concat(new T[] { t2 }));
        }

        public static IEnumerable<IEnumerable<T>> GetPermutationsWithRept<T>(IEnumerable<T> list, int length)
        {
            if (length == 1) return list.Select(t => new T[] { t });
            return GetPermutationsWithRept(list, length - 1).SelectMany(t => list, (t1, t2) => t1.Concat(new T[] { t2 }));
        }

        //Functions are not used now.      

        private bool SameComponentGlycan(Glycan glycan)
        {
            return this.Kind == glycan.Kind;
        }

        private static Node GobackRootNode(Node node)
        {
            while (node.Father != null)
            {
                node = node.Father;
            }
            return node;
        }

        private static SortedSet<int> GetAllChildrenMass(Node node)
        {
            SortedSet<int> masses = new SortedSet<int>();
            var allC = GetAllChildrenCombination(node);
            foreach (var aC in allC)
            {
                masses.Add(GetMass(Node2Struct(aC)));
            }
            return masses;
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

    public class GlycanBox
    {
        public int Mass {
            get
            {
                return Glycan.GetMass(Kind);
            }
        }

        public bool IsNGlycan;

        public List<Glycan> glycans { get; set; } = new List<Glycan>();

        public List<GlycanIon> CommonGlycanIons
        {
            get
            {
                //TO DO: A combination of glycanIons need to be considered.
                List<GlycanIon> glycanIons = new List<GlycanIon>();
                HashSet<double> ionMasses = new HashSet<double>();
                foreach (var ion in glycans.SelectMany(p=>p.Ions))
                {
                    if (!ionMasses.Contains(ion.IonMass))
                    {
                        ionMasses.Add(ion.IonMass);
                        var lossIonMass = Glycan.GetIonLossMass(Kind, ion.IonKind);
                        GlycanIon glycanIon = new GlycanIon(null, ion.IonMass, ion.IonKind, lossIonMass);
                        glycanIons.Add(glycanIon);
                    }
                }

                return glycanIons;
            }
        }

        public int NumberOfGlycans { get { return glycans.Count; } }

        public byte[] Kind
        {
            get
            {          
                if (IsNGlycan)
                {
                    return glycans.First().Kind;
                }
                else
                {
                    byte[] kind = new byte[5] { 0, 0, 0, 0, 0 };
                    foreach (var aglycan in glycans)
                    {
                        for (int i = 0; i < 5; i++)
                        {
                            kind[i] += aglycan.Kind[i];
                        }
                    }
                    return kind;
                }                                 
            }
        }

    }

}