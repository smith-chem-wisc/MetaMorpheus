using System.Collections.Generic;

namespace EngineLayer
{
    public class Glycan
    {
        public Glycan(int glyId, int glyType, string struc, double mass, int[] kind, List<GlycanIon> ions )
        {
            GlyId = glyId;
            GlyType = glyType;
            Struc = struc;
            Mass = mass;
            Kind = kind;
            Ions = ions;
        }
        public int GlyId { get; set; }
        public int GlyType { get; set; }
        public string Struc { get; set; }
        public double Mass { get; set; }
        public int[] Kind { get; set; }
        public List<GlycanIon> Ions {get; set;}
        public Dictionary<int, double> GetDiagnosticIons()
        {
            Dictionary<int, double> diagnosticIons = new Dictionary<int, double>();
            if (Kind[1] >= 1)
            {
                diagnosticIons.Add(126, 126.055);
                diagnosticIons.Add(138, 138.055);
                diagnosticIons.Add(144, 144.065);
                diagnosticIons.Add(168, 168.066);
                diagnosticIons.Add(186, 186.076);
                diagnosticIons.Add(204, 204.087);
            }
            if (Kind[1]>=1 && Kind[0]>=1)
            {
                diagnosticIons.Add(366, 366.140);
            }
            if (Kind[2]>=1)
            {
                diagnosticIons.Add(274, 274.092);
                diagnosticIons.Add(292, 292.103);
            }
            return diagnosticIons;
        }
    }
    public class GlycanIon
    {
        public GlycanIon(int ionStruct, double ionMass, int[] ionKind)
        {
            IonStruct = ionStruct;
            IonMass = ionMass;
            IonKind = ionKind;
        }
        public int IonStruct { get; set; }
        public double IonMass { get; set; }
        public int[] IonKind { get; set; }
    }
}
