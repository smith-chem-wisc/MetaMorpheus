using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class LibraryTargetDecoyDictionary
    {
        public LibraryTargetDecoyDictionary()
        {

        }
        public Dictionary<PeptideWithSetModifications, PeptideWithSetModifications> TargetDecoyPairs;

    }
}
