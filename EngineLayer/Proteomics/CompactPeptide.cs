using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide
    {

        #region Public Fields

        public readonly byte[] BaseSequence;
        public readonly ushort varMod1Type;
        public readonly ushort varMod1Loc;
        public readonly ushort varMod2Type;
        public readonly ushort varMod2Loc;
        public readonly ushort varMod3Type;
        public readonly ushort varMod3Loc;
        public double MonoisotopicMassIncludingFixedMods;

        #endregion Public Fields

        #region Public Constructors

        public CompactPeptide(PeptideWithSetModifications yyy, Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            varMod1Type = 0;
            varMod1Loc = 0;
            varMod2Type = 0;
            varMod2Loc = 0;
            varMod3Type = 0;
            varMod3Loc = 0;
            foreach (var kvpp in yyy.allModsOneIsNterminus)
            {
                var ok = modsDictionary[kvpp.Value];
                if (ok > 0)
                {
                    ushort oneBasedLoc = (ushort)kvpp.Key;
                    var mod = kvpp.Value;
                    if (varMod1Type == 0)
                    {
                        varMod1Type = ok;
                        varMod1Loc = oneBasedLoc;
                    }
                    else if (varMod2Type == 0)
                    {
                        varMod2Type = ok;
                        varMod2Loc = oneBasedLoc;
                    }
                    else
                    {
                        varMod3Type = ok;
                        varMod3Loc = oneBasedLoc;
                    }
                }
            }
            MonoisotopicMassIncludingFixedMods = yyy.MonoisotopicMass;

            BaseSequence = yyy.BaseSequence.Select(b => (byte)b).ToArray();
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Equals(object obj)
        {
            var cp = obj as CompactPeptide;
            if (cp == null)
                return false;
            return (BaseSequence.SequenceEqual(cp.BaseSequence) &&
                    varMod1Type == cp.varMod1Type &&
                    varMod1Loc == cp.varMod1Loc &&
                    varMod2Type == cp.varMod2Type &&
                    varMod2Loc == cp.varMod2Loc &&
                    varMod3Type == cp.varMod3Type &&
                    varMod3Loc == cp.varMod3Loc);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                var result = 0;
                foreach (byte b in BaseSequence)
                    result = (result * 31) ^ b;
                return result;
            }
        }

        #endregion Public Methods

    }
}