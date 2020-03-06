using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using Proteomics.Fragmentation;

namespace EngineLayer
{
    public class ModBox
    {
        public ModBox(int[] ids)
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = true;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass { get; set; }

        public double DecoyMass { get; set; }

        public bool TargetDecoy { get; set; }

    }
}
