using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class TheoreticalProducts
    {
        public List<Product> products { get; set; } = new List<Product>();

        public int position { get; set; }

        public List<int> positions { get; set; }
    }
}
