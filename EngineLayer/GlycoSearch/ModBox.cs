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
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public virtual double Mass { get; set; }

    }
}
