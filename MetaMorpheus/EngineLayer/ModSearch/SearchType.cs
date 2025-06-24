using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ModSearch
{
    public enum SearchType
    {
        O_GlycanSearch, //consider O-glycan modification
        N_GlycanSearch, // consider N-glycan modification
        RegularModSearch, //Only consider regular modification
        MixedModSearch, // consider all modifications, both O-glycan, n-glycan and regular modifications
    }
}
