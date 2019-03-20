using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;

namespace EngineLayer
{
    public class Glycan
    {
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        //H: C6O5H10, N: C8O5NH13, A: C11O8NH17, G: C11H17NO9, F: C6O4H10
        private static Dictionary<char, double> CharMassDic = new Dictionary<char, double>() { { 'H', 162.052823 }, { 'N', 203.079373 }, { 'A', 291.095417 }, { 'G', 307.090331 }, { 'F', 146.057909 } };

        public static Dictionary<int, double> oxoniumIons = new Dictionary<int, double>()
        {
            { 126, 126.055 },
            { 138, 138.055 },
            { 144, 144.065 },
            { 168, 168.066 },
            { 186, 186.076 },
            { 204, 204.087 },
            { 366, 366.140 },
            { 274, 274.092 },
            { 292, 292.103 }
        };

        public static Dictionary<int, double> TrimannosylCores = new Dictionary<int, double>()
        {
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

        public Glycan(string struc, double mass, byte[] kind, List<GlycanIon> ions)
        {
            Struc = struc;
            Mass = mass;
            Kind = kind;
            Ions = ions;
        }
        public int GlyId { get; set; }
        public int GlyType { get; set; }
        public string Struc { get; set; }
        public double Mass { get; set; }
        public byte[] Kind { get; set; }
        public List<GlycanIon> Ions { get; set; }
        
        private static Node Struct2Node(string theGlycanStruct)
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
                        curr = curr.father;
                        level--;
                    }
                    else
                    {
                        level++;
                        if (curr.lChild == null)
                        {
                            curr.lChild = new Node(theGlycanStruct[i], level);
                            curr.lChild.father = curr;
                            curr = curr.lChild;
                        }
                        else
                        {
                            curr.rChild = new Node(theGlycanStruct[i], level);
                            curr.rChild.father = curr;
                            curr = curr.rChild;
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
                output += "(" + node.value + Node2Struct(node.lChild) + Node2Struct(node.rChild) + ")";
            }
            return output;
        }

        private static List<Node> GetAllChildrenCombination(Node node)
        {
            List<Node> nodes = new List<Node>();
            var curr = node;
            if (curr.lChild == null && curr.rChild == null)
            {
                nodes.Add(curr);
            }
            else
            {
                List<Node> l = GetAllChildrenCombination(curr.lChild);
                nodes.Add(new Node(curr.value));
                if (curr.rChild != null)
                {
                    List<Node> r = GetAllChildrenCombination(curr.rChild);
                    foreach (var lNode in l)
                    {
                        var c = new Node(curr.value);
                        c.lChild = lNode;
                        nodes.Add(c);
                    }
                    foreach (var rNode in r)
                    {
                        var c = new Node(curr.value);
                        c.rChild = rNode;
                        nodes.Add(c);
                    }
                    foreach (var lNode in l)
                    {
                        foreach (var rNode in r)
                        {
                            var c = new Node(curr.value);
                            c.lChild = lNode;
                            c.rChild = rNode;
                            nodes.Add(c);
                        }
                    }
                }
                else
                {
                    foreach (var lNode in l)
                    {
                        var c = new Node(curr.value);
                        c.lChild = lNode;
                        nodes.Add(c);
                    }
                }
            }

            return nodes; 
        }

        private static double GetMass(string structure)
        {
            double y = CharMassDic['H'] * structure.Count(p => p == 'H') +
                CharMassDic['N'] * structure.Count(p => p == 'N') +
                CharMassDic['A'] * structure.Count(p => p == 'A') +
                CharMassDic['G'] * structure.Count(p => p == 'G') +
                CharMassDic['F'] * structure.Count(p => p == 'F');
            return y;
        }

        private static byte[] GetKind(string structure)
        {
            byte[] kind = new byte[] { Convert.ToByte(structure.Count(p => p == 'H')), Convert.ToByte(structure.Count(p => p == 'N')), Convert.ToByte(structure.Count(p => p == 'A')) , Convert.ToByte(structure.Count(p => p == 'G')) , Convert.ToByte(structure.Count(p => p == 'F')) };
            return kind;
        }

        private static SortedSet<double> GetAllChildrenMass(Node node)
        {
            SortedSet<double> masses = new SortedSet<double>();
            var allC = GetAllChildrenCombination(node);
            foreach (var aC in allC)
            {
                masses.Add(GetMass(Node2Struct(aC)));
            }
            return masses;
        }

        public static double GetMass(byte[] kind)
        {
            double y = CharMassDic['H'] * kind[0] +
                CharMassDic['N'] * kind[1] +
                CharMassDic['A'] * kind[2] +
                CharMassDic['G'] * kind[3] +
                CharMassDic['F'] * kind[4];
            return y;
        }

        public Dictionary<int, double> GetDiagnosticIons()
        {
            Dictionary<int, double> diagnosticIons = new Dictionary<int, double>();
            if (Kind[1] >= 1)
            {
                diagnosticIons.Add(126, 126.055 - hydrogenAtomMonoisotopicMass);
                diagnosticIons.Add(138, 138.055 - hydrogenAtomMonoisotopicMass);
                diagnosticIons.Add(144, 144.065 - hydrogenAtomMonoisotopicMass);
                diagnosticIons.Add(168, 168.066 - hydrogenAtomMonoisotopicMass);
                diagnosticIons.Add(186, 186.076 - hydrogenAtomMonoisotopicMass);
                diagnosticIons.Add(204, 204.087 - hydrogenAtomMonoisotopicMass);
            }
            if (Kind[1] >= 1 && Kind[0] >= 1)
            {
                diagnosticIons.Add(366, 366.140 - hydrogenAtomMonoisotopicMass);
            }
            if (Kind[2] >= 1)
            {
                diagnosticIons.Add(274, 274.092 - hydrogenAtomMonoisotopicMass);
                diagnosticIons.Add(292, 292.103 - hydrogenAtomMonoisotopicMass);
            }
            return diagnosticIons;
        }

        public static Glycan Struct2Glycan(string theGlycanStruct, int id)
        {
            Node node = Struct2Node(theGlycanStruct);
            List<Node> nodeIons = GetAllChildrenCombination(node);
            double mass = GetMass(theGlycanStruct);
            byte[] kind = GetKind(theGlycanStruct);
            List<GlycanIon> glycanIons = new List<GlycanIon>();
            HashSet<double> ionMasses = new HashSet<double>();
            foreach (var aNodeIon in nodeIons)
            {
                var ionMass = GetMass(Node2Struct(aNodeIon));
                if (!ionMasses.Contains(ionMass))
                {
                    ionMasses.Add(ionMass);
                    var ionKind = GetKind(Node2Struct(aNodeIon));
                    GlycanIon glycanIon = new GlycanIon(0, ionMass, ionKind);
                    glycanIons.Add(glycanIon);
                }
            }
            var halfIonKind = new byte[] { 0, 0, 0, 0, 0 };
            glycanIons.Add(new GlycanIon(0, 83.038194, halfIonKind)); //Cross-ring mass
            glycanIons = glycanIons.OrderBy(p => p.IonMass).ToList();
            //glycanIons.RemoveAt(glycanIons.Count - 1);

            Glycan glycan = new Glycan(theGlycanStruct, mass, kind, glycanIons);
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

        public static string GetKindString(byte[] Kind)
        {
            string H = (Kind[0] > 0) ? "H" + Kind[0].ToString() : "";
            string N = (Kind[1] > 0) ? "N" + Kind[1].ToString() : "";
            string A = (Kind[2] > 0) ? "A" + Kind[2].ToString() : "";
            string G = (Kind[3] > 0) ? "G" + Kind[3].ToString() : "";
            string F = (Kind[4] > 0) ? "F" + Kind[4].ToString() : "";
            string kindString = H + N + A + G + F;
            if (kindString == "")
            {
                kindString = "@";
            }
            return kindString;
        }

        public static string GetKindString(string structure)
        {
            string H = (structure.Count(p => p == 'H') > 0) ? "H" + structure.Count(p => p == 'H').ToString() : "";
            string N = (structure.Count(p => p == 'N') > 0) ? "N" + structure.Count(p => p == 'N').ToString() : "";
            string A = (structure.Count(p => p == 'A') > 0) ? "A" + structure.Count(p => p == 'A').ToString() : "";
            string G = (structure.Count(p => p == 'G') > 0) ? "G" + structure.Count(p => p == 'G').ToString() : "";
            string F = (structure.Count(p => p == 'F') > 0) ? "F" + structure.Count(p => p == 'F').ToString() : "";
            string kindString = H + N + A + G + F;
            if (kindString == "")
            {
                kindString = "@";
            }
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
            while (node.father != null)
            {
                node = node.father;
            }
            return node;
        }

    }

    public class GlycanIon
    {
        public GlycanIon(int ionStruct, double ionMass, byte[] ionKind)
        {
            IonStruct = ionStruct;
            IonMass = ionMass;
            IonKind = ionKind;
        }
        public int IonStruct { get; set; }
        public double IonMass { get; set; }
        public byte[] IonKind { get; set; }
    }

    public class GlycanBox
    {
        public double Mass {
            get
            {
                return Glycan.GetMass(Kind);
            }
        }

        public List<Glycan> glycans { get; set; } = new List<Glycan>();

        public List<GlycanIon> CommonGlycanIons { get; set; }

        public int NumberOfGlycans { get { return glycans.Count; } }

        public byte[] Kind
        {
            get
            {
                byte[] kind = new byte[5] { 0,0,0,0,0};
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