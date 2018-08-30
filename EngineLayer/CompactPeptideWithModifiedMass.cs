using System;
using System.Collections.Generic;
using System.Text;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer
{
    //[Serializable]
    //public class CompactPeptideWithModifiedMass : CompactPeptideBase
    //{

    //    public CompactPeptideWithModifiedMass(CompactPeptideBase cp, double monoisotopicMassIncludingFixedMods)
    //    {
    //        CTerminalMasses = cp.CTerminalMasses;
    //        NTerminalMasses = cp.NTerminalMasses;
    //        MonoisotopicMassIncludingFixedMods = cp.MonoisotopicMassIncludingFixedMods;
    //        ModifiedMass = monoisotopicMassIncludingFixedMods;
    //    }

    //    public double ModifiedMass { get; set; }

    //    public void SwapMonoisotopicMassWithModifiedMass()
    //    {
    //        double tempDouble = MonoisotopicMassIncludingFixedMods;
    //        MonoisotopicMassIncludingFixedMods = ModifiedMass;
    //        ModifiedMass = tempDouble;
    //    }

    //    public void CropTerminalMasses(FragmentationTerminus terminusType)
    //    {
    //        List<double> tempList = new List<double>();
    //        double[] masses = terminusType == FragmentationTerminus.N ? NTerminalMasses : CTerminalMasses;
    //        for (int i = 0; i < masses.Length; i++)
    //        {
    //            if (masses[i] < MonoisotopicMassIncludingFixedMods)
    //            {
    //                tempList.Add(masses[i]);
    //            }
    //            else if (terminusType == FragmentationTerminus.N)
    //            {
    //                NTerminalMasses = tempList.ToArray();
    //                break;
    //            }
    //            else
    //            {
    //                CTerminalMasses = tempList.ToArray();
    //                break;
    //            }
    //        }
    //    }
    //}

    [Serializable]
    public class PeptideWithSetModifications_ModifiedMass
    {
        public double[] CterminalMasses { get; private set; }
        public double[] NterminalMasses { get; private set; }
        public double PeptideWithSetModifications_MonoisotopicMass { get; private set; }
        public double ModifiedMass { get; private set; }


        public PeptideWithSetModifications_ModifiedMass(PeptideWithSetModifications p, double newMonoMass, MassSpectrometry.DissociationType d)
        {
            var c = p.Fragment(d, FragmentationTerminus.C);
            List<double> l = new List<double>();
            foreach (Product prod in c)
            {
                l.Add(prod.NeutralMass);
            }
            this.CterminalMasses = l.ToArray();

            c = p.Fragment(d, FragmentationTerminus.N);
            l = new List<double>();
            foreach (Product prod in c)
            {
                l.Add(prod.NeutralMass);
            }
            this.NterminalMasses = l.ToArray();
            this.PeptideWithSetModifications_MonoisotopicMass = p.MonoisotopicMass;
            this.ModifiedMass = newMonoMass;
        }

        public void SwapMonoisotopicMassWithModifiedMass()
        {
            double tempDouble = PeptideWithSetModifications_MonoisotopicMass;
            PeptideWithSetModifications_MonoisotopicMass = ModifiedMass;
            ModifiedMass = tempDouble;
        }

        public void CropTerminalMasses(FragmentationTerminus terminusType)
        {
            List<double> tempList = new List<double>();
            double[] masses = terminusType == FragmentationTerminus.N ? NterminalMasses : CterminalMasses;
            for (int i = 0; i < masses.Length; i++)
            {
                if (masses[i] < PeptideWithSetModifications_MonoisotopicMass)
                {
                    tempList.Add(masses[i]);
                }
                else if (terminusType == FragmentationTerminus.N)
                {
                    NterminalMasses = tempList.ToArray();
                    break;
                }
                else
                {
                    CterminalMasses = tempList.ToArray();
                    break;
                }
            }
        }

    }

}
