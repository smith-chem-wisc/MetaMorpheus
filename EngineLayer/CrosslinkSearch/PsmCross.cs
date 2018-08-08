using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : PeptideSpectralMatch
    {
        public CompactPeptide compactPeptide;

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams) : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams)
        {
            compactPeptide = theBestPeptide;
        }

        public double BestScore { get; set; } //For the current psmCross
        public int ModPosition { get; set; }
        public PsmCross BetaPsmCross { get; set; }
        public double DScore { get; set; }
        public CrosslinkInfo CrosslinkInfo { get; set; }
        public Glycan Glycan {get; set;}
        public List<MatchedFragmentIon> MatchedIons { get; set; }

        private List<TheoreticalFragmentIon> GetTheoreticalFragmentIons(List<ProductType> productTypes)
        {
            List<TheoreticalFragmentIon> theoreticalFragmentIons = new List<TheoreticalFragmentIon>();
            foreach (var pt in productTypes)
            {
                int ionNumberAdd = 1;
                if (pt == ProductType.BnoB1ions)
                {
                    // first generated b ion is b2, not b1, if we're skipping b1 ions
                    ionNumberAdd++;
                }
                List<ProductType> temp = new List<ProductType> { pt };
                var productMasses = compactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(temp);

                for (int i = 0; i < productMasses.Length; i++)
                {
                    theoreticalFragmentIons.Add(new TheoreticalFragmentIon(productMasses[i], double.NaN, 1, pt, i + ionNumberAdd));
                }
            }


            return theoreticalFragmentIons;
        }

        public Dictionary<int, List<TheoreticalFragmentIon>> XlGetTheoreticalFramentIons(List<ProductType> productTypes, bool Charge_2_3, CrosslinkerTypeClass crosslinker,List<int> modPos, double modMass)
        {
            Dictionary<int, List<TheoreticalFragmentIon>> AllTheoreticalFragmentIonsLists = new Dictionary<int, List<TheoreticalFragmentIon>>();

            List<TheoreticalFragmentIon> baseTheoreticalFragmentIons = GetTheoreticalFragmentIons(productTypes);

            foreach (var iPos in modPos)
            {
                List<TheoreticalFragmentIon> currentIons = new List<TheoreticalFragmentIon>();

                foreach (var iIon in baseTheoreticalFragmentIons)
                {
                    var iType = iIon.ProductType;
                    switch (iType)
                    {
                        case ProductType.BnoB1ions:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else
                            {
                                currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                if (crosslinker.Cleavable)
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassShort, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassLong, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                }
                            }
                            break;
                        case ProductType.C:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            break;
                        case ProductType.Y:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else
                            {
                                currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                if (crosslinker.Cleavable)
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassShort, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassLong, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                }
                            }
                            break;
                        case ProductType.Zdot:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            break;
                    }
                }

                if (Charge_2_3)
                {
                    var length = currentIons.Count;
                    for (int i = 0; i < length; i++)
                    {
                        currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 2, currentIons[i].ProductType, currentIons[i].IonNumber));
                        currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 3, currentIons[i].ProductType, currentIons[i].IonNumber));
                    }
                }

                AllTheoreticalFragmentIonsLists.Add(iPos, currentIons);
            }

            return AllTheoreticalFragmentIonsLists;
        }

        //TO DO: the second ModPostion jPos is not recorded. 
        public Dictionary<int, List<TheoreticalFragmentIon>> XlLoopGetTheoreticalFramentIons(List<ProductType> productTypes, bool Charge_2_3, CrosslinkerTypeClass crosslinker, List<int> modPos, double modMass)
        {
            Dictionary<int, List<TheoreticalFragmentIon>> AllTheoreticalFragmentIonsLists = new Dictionary<int, List<TheoreticalFragmentIon>>();

            List<TheoreticalFragmentIon> baseTheoreticalFragmentIons = GetTheoreticalFragmentIons(productTypes);

            if (modPos.Count() >= 2)
            {
                for (int iPos = 0; iPos < modPos.Count() - 1; iPos++)
                {
                    for (int jPos = iPos + 1; jPos < modPos.Count(); jPos++)
                    {
                        List<TheoreticalFragmentIon> currentIons = new List<TheoreticalFragmentIon>();

                        foreach (var iIon in baseTheoreticalFragmentIons)
                        {
                            var iType = iIon.ProductType;
                            switch (iType)
                            {
                                case ProductType.BnoB1ions:
                                    if (iIon.IonNumber < iPos)
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if(iIon.IonNumber >= jPos)
                                    {
                                        currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                    }
                                    break;
                                case ProductType.C:
                                    if (iIon.IonNumber < iPos)
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if(iIon.IonNumber >= jPos) { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                                    break;
                                case ProductType.Y:
                                    if (iIon.IonNumber < iPos)
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if(iIon.IonNumber >= jPos)
                                    {
                                        currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));          
                                    }
                                    break;
                                case ProductType.Zdot:
                                    if (iIon.IonNumber < iPos)
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if (iIon.IonNumber >= jPos) { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                                    break;
                            }
                        }

                        if (Charge_2_3)
                        {
                            var length = currentIons.Count;
                            for (int i = 0; i < length; i++)
                            {
                                currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 2, currentIons[i].ProductType, currentIons[i].IonNumber));
                                currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 3, currentIons[i].ProductType, currentIons[i].IonNumber));
                            }
                        }
                        AllTheoreticalFragmentIonsLists.Add(iPos, currentIons);
                    }
                }
            }

            return AllTheoreticalFragmentIonsLists;
        }

        public static int[] GenerateIntensityRanks(double[] experimental_intensities)
        {
            var y = experimental_intensities.ToArray();
            var x = Enumerable.Range(1, y.Length).OrderBy(p => p).ToArray();
            Array.Sort(y, x);
            var experimental_intensities_rank = Enumerable.Range(1, y.Length).OrderByDescending(p => p).ToArray();
            Array.Sort(x, experimental_intensities_rank);
            return experimental_intensities_rank;
        }

        public List<int> XlPosCal(string crosslinkerModSites)
        {
            Tolerance tolerance = new PpmTolerance(1);
            List<int> xlpos = new List<int>();
            foreach (char item in crosslinkerModSites)
            {
                if (tolerance.Within(compactPeptide.NTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(0);
                }
                for (int i = 1; i < compactPeptide.NTerminalMasses.Length; i++)
                {
                    if (tolerance.Within(compactPeptide.NTerminalMasses[i] - compactPeptide.NTerminalMasses[i - 1], Residue.GetResidue(item).MonoisotopicMass))
                    {
                        xlpos.Add(i);
                    }
                }
                if (tolerance.Within(compactPeptide.CTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(compactPeptide.NTerminalMasses.Length);
                }
            }
            xlpos.Sort();
            return xlpos;
        }

        public Dictionary<int, List<TheoreticalFragmentIon>> GlyGetTheoreticalFramentIons(List<ProductType> productTypes, bool Charge_2_3, List<int> modPos)
        {
            var diagnosticIons = Glycan.GetDiagnosticIons();

            Dictionary<int, List<TheoreticalFragmentIon>> AllTheoreticalFragmentIonsLists = new Dictionary<int, List<TheoreticalFragmentIon>>();

            List<TheoreticalFragmentIon> baseTheoreticalFragmentIons = GetTheoreticalFragmentIons(productTypes);

            double modMass = 0;
            if (Glycan.Ions.Count != 0)
            {
                modMass = Glycan.Ions.First().IonMass;
            }

            foreach (var iPos in modPos)
            {
                List<TheoreticalFragmentIon> currentIons = new List<TheoreticalFragmentIon>();
                foreach (var ion in diagnosticIons)
                {
                    currentIons.Add(new TheoreticalFragmentIon(ion.Value - hydrogenAtomMonoisotopicMass, double.NaN, 1, ProductType.None, ion.Key));
                }
                
                foreach (var iIon in baseTheoreticalFragmentIons)
                {
                    var iType = iIon.ProductType;
                    switch (iType)
                    {
                        case ProductType.BnoB1ions:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else
                            {
                                currentIons.Add(iIon);
                                if (Glycan.Ions.Count != 0)
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                }
                            }
                            break;
                        case ProductType.C:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + Glycan.Mass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            break;
                        case ProductType.Y:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else
                            {
                                currentIons.Add(iIon);
                                if (Glycan.Ions.Count != 0)
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));

                                }
                            }
                            break;
                        case ProductType.Zdot:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + Glycan.Mass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            break;
                    }
                }

                for (int i = 0; i < Glycan.Ions.Count; i++)
                {
                    currentIons.Add(new TheoreticalFragmentIon(Glycan.Ions[i].IonMass + compactPeptide.MonoisotopicMassIncludingFixedMods, double.NaN, 1, ProductType.X, i));
                }

                if (Charge_2_3)
                {
                    var length = currentIons.Count;
                    for (int i = 0; i < length; i++)
                    {
                        currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 2, currentIons[i].ProductType, currentIons[i].IonNumber));
                        currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 3, currentIons[i].ProductType, currentIons[i].IonNumber));
                    }
                }

                AllTheoreticalFragmentIonsLists.Add(iPos, currentIons);
            }

            return AllTheoreticalFragmentIonsLists;
        }

        //To do: Optimize motif
        public List<int> GlyPosCal(CompactPeptide compactPeptide, string modSites)
        {
            Tolerance tolerance = new PpmTolerance(1);
            List<int> xlpos = new List<int>();
            foreach (char item in modSites)
            {
                if (tolerance.Within(compactPeptide.NTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(0);
                }
                for (int i = 1; i < compactPeptide.NTerminalMasses.Length; i++)
                {
                    if (tolerance.Within(compactPeptide.NTerminalMasses[i] - compactPeptide.NTerminalMasses[i - 1], Residue.GetResidue(item).MonoisotopicMass))
                    {
                        xlpos.Add(i);
                    }
                }
                if (tolerance.Within(compactPeptide.CTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(compactPeptide.NTerminalMasses.Length);
                }
            }
            xlpos.Sort();
            return xlpos;
        }

        public void GetBestMatch(Ms2ScanWithSpecificMass theScan, List<ProductType> productTypes, bool Charge_2_3, CommonParameters commonParameters, List<int> modPos)
        {
            var pmmhList = GlyGetTheoreticalFramentIons(productTypes, Charge_2_3, modPos);

            BestScore = 0;
            foreach (var pmmh in pmmhList)
            {
                var matchedIons = MetaMorpheusEngine.MatchFragmentIons( theScan.TheScan.MassSpectrum, pmmh.Value, commonParameters);
                var score = MetaMorpheusEngine.CalculatePeptideScore(theScan.TheScan, matchedIons, 0);
                if (score > BestScore)
                {
                    BestScore = score;
                    MatchedIons = matchedIons;
                    ModPosition = pmmh.Key;
                }
            }
        }
    }
}