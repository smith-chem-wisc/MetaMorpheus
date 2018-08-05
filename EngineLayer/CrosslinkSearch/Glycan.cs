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
