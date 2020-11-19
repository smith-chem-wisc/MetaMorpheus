using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public enum LibraryDecoyType
    {
        /// <summary>
        /// Generate no decoy
        /// </summary>
        None,

        SwapPrecurser,

        /// <summary>
        /// Reverse the protein sequence, possibly keeping the initiating methionine in place
        /// </summary>
        Reverse,

        /// <summary>
        /// No clue...
        /// </summary>
        Slide,

        /// <summary>
        /// Generate decoy by:
        /// 1. simulating proteolytic digesiton (if any)
        /// 2. shuffling the resulting peptides, keeping the cleavage site and possibly initiating methionine in place,
        /// 3. concatenating them back together into
        /// </summary>
        Shuffle,

        /// <summary>
        /// No clue... not implemented
        /// </summary>
        Random
    }
}
