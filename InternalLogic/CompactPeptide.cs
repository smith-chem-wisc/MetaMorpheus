using OldInternalLogic;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    [Serializable]
    public class CompactPeptide
    {

        #region Public Fields

        public readonly byte[] BaseSequence;
        public readonly ushort varMod1Type;
        public readonly byte varMod1Loc;
        public readonly ushort varMod2Type;
        public readonly byte varMod2Loc;
        public readonly ushort varMod3Type;
        public readonly byte varMod3Loc;
        public float MonoisotopicMass;

        #endregion Public Fields

        #region Public Constructors

        public CompactPeptide(PeptideWithSetModifications yyy, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            varMod1Type = 0;
            varMod1Loc = 0;
            varMod2Type = 0;
            varMod2Loc = 0;
            varMod3Type = 0;
            varMod3Loc = 0;
            foreach (var kvpp in yyy.twoBasedVariableAndLocalizeableModificationss)
            {
                int twoBasedLoc = kvpp.Key;
                var mod = kvpp.Value;
                // Variable
                if (varMod1Type == 0)
                {
                    // Set first
                    if (variableModifications.Contains(mod))
                        varMod1Type = (ushort)(variableModifications.IndexOf(mod) + 1);
                    else
                        varMod1Type = (ushort)(32767 + localizeableModifications.IndexOf(mod) + 1);
                    varMod1Loc = (byte)twoBasedLoc;
                }
                else if (varMod2Type == 0)
                {
                    // Set second
                    if (variableModifications.Contains(mod))
                        varMod2Type = (ushort)(variableModifications.IndexOf(mod) + 1);
                    else
                        varMod2Type = (ushort)(32767 + localizeableModifications.IndexOf(mod) + 1);
                    varMod2Loc = (byte)twoBasedLoc;
                }
                else
                {
                    // Set third
                    if (variableModifications.Contains(mod))
                        varMod3Type = (ushort)(variableModifications.IndexOf(mod) + 1);
                    else
                        varMod3Type = (ushort)(32767 + localizeableModifications.IndexOf(mod) + 1);
                    varMod3Loc = (byte)twoBasedLoc;
                }
            }

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