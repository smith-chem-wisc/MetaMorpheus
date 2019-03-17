using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkedPeptide
    {
        public static List<Tuple<double, List<Product>>> XlGetTheoreticalFragments(DissociationType dissociationType, Crosslinker crosslinker,
            int crosslinkerPosition, double otherPeptideMass, PeptideWithSetModifications peptide)
        {
            List<double> massesToLocalize = new List<double>();
            if (crosslinker.Cleavable)
            {
                massesToLocalize.Add(crosslinker.CleaveMassShort);
                massesToLocalize.Add(crosslinker.CleaveMassLong);
            }

            //The cleavage of crosslinker depend on the crosslinker structure, dissociation type or energy. In theory it won't always break.
            massesToLocalize.Add(crosslinker.TotalMass + otherPeptideMass);

            List<Tuple<double, List<Product>>> massToLocalizeWithTheoreticalProducts = new List<Tuple<double, List<Product>>>();
            
            foreach (double massToLocalize in massesToLocalize)
            {
                HashSet<double> masses = new HashSet<double>();
                List<Product> theoreticalProducts = new List<Product>();
                Dictionary<int, Modification> testMods = new Dictionary<int, Modification> { { crosslinkerPosition + 1, new Modification(_monoisotopicMass: massToLocalize) } };

                foreach (var mod in peptide.AllModsOneIsNterminus)
                {
                    testMods.Add(mod.Key, mod.Value);
                }

                var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                    peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

                // add fragmentation ions for this crosslinker position guess
                foreach (var fragment in testPeptide.Fragment(dissociationType, FragmentationTerminus.Both))
                {
                    if (!masses.Contains(fragment.NeutralMass))
                    {
                        theoreticalProducts.Add(fragment);
                        masses.Add(fragment.NeutralMass);
                    }
                }

                // add signature ions
                if (crosslinker.Cleavable)
                {
                    theoreticalProducts.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.None, peptide.MonoisotopicMass + massToLocalize,
                        peptide.Length, peptide.Length), 0));
                }

                massToLocalizeWithTheoreticalProducts.Add(new Tuple<double, List<Product>>(massToLocalize, theoreticalProducts));
            }

            return massToLocalizeWithTheoreticalProducts;

        }

        public static Dictionary<Tuple<int, int>, List<Product>> XlLoopGetTheoreticalFragments(DissociationType dissociationType, Modification loopMass,
            List<int> modPos, PeptideWithSetModifications peptide)
        {
            Dictionary<Tuple<int, int>, List<Product>> AllTheoreticalFragmentsLists = new Dictionary<Tuple<int, int>, List<Product>>();
            var originalFragments = peptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();

            foreach (int position1 in modPos)
            {
                foreach (int position2 in modPos)
                {
                    if (position2 <= position1)
                    {
                        continue;
                    }

                    // add N and C terminal fragments that do not contain the loop
                    Tuple<int, int> loopPositions = new Tuple<int, int>(position1, position2);
                    List<Product> loopFragments = originalFragments
                        .Where(p => p.TerminusFragment.Terminus == FragmentationTerminus.N && p.TerminusFragment.AminoAcidPosition < position1
                        || p.TerminusFragment.Terminus == FragmentationTerminus.C && p.TerminusFragment.AminoAcidPosition > position2).ToList();

                    // add N-terminal fragments containing the loop
                    Dictionary<int, Modification> modDict = new Dictionary<int, Modification>();
                    if (peptide.AllModsOneIsNterminus.Any())
                    {
                        double combinedModMass = loopMass.MonoisotopicMass.Value + peptide.AllModsOneIsNterminus.Where(v => v.Key <= position2 + 1).Sum(p => p.Value.MonoisotopicMass.Value);
                        Modification combined = new Modification(_monoisotopicMass: combinedModMass);
                        modDict.Add(position1 + 1, combined);

                        foreach (var mod in peptide.AllModsOneIsNterminus.Where(m => m.Key > position2 + 1))
                        {
                            modDict.Add(mod.Key, mod.Value);
                        }
                    }
                    else
                    {
                        modDict.Add(position1 + 1, loopMass);
                    }
                    PeptideWithSetModifications peptideWithLoop = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams,
                        peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory,
                        peptide.PeptideDescription, peptide.MissedCleavages, modDict, peptide.NumFixedMods);
                    loopFragments.AddRange(peptideWithLoop.Fragment(dissociationType, FragmentationTerminus.Both)
                        .Where(p => p.TerminusFragment.Terminus == FragmentationTerminus.N && p.TerminusFragment.AminoAcidPosition >= position2));

                    // add C-terminal fragments containing the loop
                    modDict.Clear();
                    if (peptide.AllModsOneIsNterminus.Any())
                    {
                        double combinedModMass = loopMass.MonoisotopicMass.Value + peptide.AllModsOneIsNterminus.Where(v => v.Key >= position1 + 1).Sum(p => p.Value.MonoisotopicMass.Value);
                        Modification combined = new Modification(_monoisotopicMass: combinedModMass);
                        modDict.Add(position2 + 1, combined);

                        foreach (var mod in peptide.AllModsOneIsNterminus.Where(m => m.Key < position1 + 1))
                        {
                            modDict.Add(mod.Key, mod.Value);
                        }
                    }
                    else
                    {
                        modDict.Add(position2 + 1, loopMass);
                    }
                    peptideWithLoop = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams,
                        peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory,
                        peptide.PeptideDescription, peptide.MissedCleavages, modDict, peptide.NumFixedMods);
                    loopFragments.AddRange(peptideWithLoop.Fragment(dissociationType, FragmentationTerminus.Both)
                        .Where(p => p.TerminusFragment.Terminus == FragmentationTerminus.C && p.TerminusFragment.AminoAcidPosition <= position1));

                    AllTheoreticalFragmentsLists.Add(loopPositions, loopFragments);
                }
            }

            return AllTheoreticalFragmentsLists;
        }

        //TO DO: The situation can be very complicated. Currently only MS2(h), MS2(h)-MS2(h), MS2(h)-MS3(l) have been considered.
        public static List<Product> XLAllowedProducts(MsDataScan TheScan, List<Tuple<double, List<Product>>> massToLocalizeWithProducts, CommonParameters commonParameters, Crosslinker crosslinker, PeptideWithSetModifications peptide)
        {
            List<Product> allowedProducts = new List<Product>();

            if (TheScan.MsnOrder == 3 && crosslinker.Cleavable)
            {
                var MS3ScanPrecursorFilter = new SinglePpmAroundZeroSearchMode(commonParameters.ProductMassTolerance.Value);
                foreach (var aMassToLocalizeWithProduct in massToLocalizeWithProducts)
                {
                    //This filter for MS3 low res scan may not working perfect. If signature ions are not selected.
                    if ((TheScan.SelectedIonChargeStateGuess.HasValue && MS3ScanPrecursorFilter.Accepts(aMassToLocalizeWithProduct.Item1 + peptide.MonoisotopicMass, TheScan.SelectedIonChargeStateGuess.Value * (TheScan.SelectedIonMZ.Value - 1.0073) )>=0)
                        || (MS3ScanPrecursorFilter.Accepts(aMassToLocalizeWithProduct.Item1 + peptide.MonoisotopicMass,  TheScan.SelectedIonMZ.Value - 1.0073) >= 0)
                        || (MS3ScanPrecursorFilter.Accepts(aMassToLocalizeWithProduct.Item1 + peptide.MonoisotopicMass, 2*(TheScan.SelectedIonMZ.Value - 1.0073)) >= 0))
                    {
                        foreach (var aProduct in aMassToLocalizeWithProduct.Item2)
                        {
                            if (!TheScan.DissociationType.HasValue || AllowedProductByDessociationType(aProduct, TheScan.DissociationType.Value))
                            {
                                allowedProducts.Add(aProduct);
                            }
                        }
                    }
                }
            }
            else
            {
                HashSet<double> masses = new HashSet<double>();
                foreach (var aProduct in massToLocalizeWithProducts.SelectMany(p=>p.Item2).ToArray())
                {
                    //!TheScan.DissociationType.HasValue or DissociationType.Unknown is used in XLtest
                    if (!TheScan.DissociationType.HasValue || AllowedProductByDessociationType(aProduct, TheScan.DissociationType.Value))
                    {
                        if (!masses.Contains(aProduct.NeutralMass))
                        {
                            allowedProducts.Add(aProduct);
                        }
                       
                    }
                }
            }
            
            return allowedProducts;
        }

        public static List<Product> SingleAllowedProducts(MsDataScan TheScan, List<Product> products)
        {
            List<Product> allowedProducts = new List<Product>();
            foreach (var aProduct in products)
            {
                if (!TheScan.DissociationType.HasValue || AllowedProductByDessociationType(aProduct, TheScan.DissociationType.Value))
                {
                    allowedProducts.Add(aProduct);
                }
            }
            return allowedProducts;
        }

        public static bool AllowedProductByDessociationType(Product product, DissociationType dissociationType)
        {
            if (product.ProductType == ProductType.M || dissociationType == DissociationType.Unknown) 
            {
                return true;
            }
            if (dissociationType == DissociationType.CID || dissociationType == DissociationType.HCD || dissociationType == DissociationType.EThcD)
            {
                if (product.ProductType == ProductType.b || product.ProductType == ProductType.y)
                {
                    return true;
                }
            }
            if (dissociationType == DissociationType.ETD || dissociationType == DissociationType.EThcD)
            {
                if (product.ProductType == ProductType.c || product.ProductType == ProductType.zDot)
                {
                    return true;
                }
            }
            return false;
        }
    }
}
