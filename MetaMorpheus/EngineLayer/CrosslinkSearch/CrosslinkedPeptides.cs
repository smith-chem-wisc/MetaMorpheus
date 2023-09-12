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
        public static IEnumerable<Tuple<int, List<Product>>> XlGetTheoreticalFragments(DissociationType dissociationType, Crosslinker crosslinker,
            List<int> possibleCrosslinkerPositions, double otherPeptideMass, PeptideWithSetModifications peptide)
        {
            List<double> massesToLocalize = new List<double>();
            if (crosslinker.Cleavable && crosslinker.CleaveDissociationTypes.Contains(dissociationType))
            {
                massesToLocalize.Add(crosslinker.CleaveMassShort);
                massesToLocalize.Add(crosslinker.CleaveMassLong);
            }
            else
            {
                massesToLocalize.Add(crosslinker.TotalMass + otherPeptideMass);
            }

            var fragments = new List<Product>();
            HashSet<double> masses = new HashSet<double>();

            foreach (int crosslinkerPosition in possibleCrosslinkerPositions)
            {
                List<Product> theoreticalProducts = new List<Product>(); //need a new one each time to pass as a reference, don't clear
                masses.Clear();

                foreach (double massToLocalize in massesToLocalize)
                {
                    Dictionary<int, Modification> testMods = new Dictionary<int, Modification> { { crosslinkerPosition + 1, new Modification(_monoisotopicMass: massToLocalize) } };

                    foreach (var mod in peptide.AllModsOneIsNterminus)
                    {
                        testMods.Add(mod.Key, mod.Value);
                    }

                    var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                    peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);
                    
                    testPeptide.Fragment(dissociationType, FragmentationTerminus.Both, fragments);

                    // add fragmentation ions for this crosslinker position guess
                    foreach (var fragment in fragments)
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
                        theoreticalProducts.Add(new Product(ProductType.M, FragmentationTerminus.None, peptide.MonoisotopicMass + massToLocalize,
                            peptide.Length, peptide.Length, 0));
                    }
                }

                yield return new Tuple<int, List<Product>>(crosslinkerPosition, theoreticalProducts);
            }
        }

        public static Dictionary<Tuple<int, int>, List<Product>> XlLoopGetTheoreticalFragments(DissociationType dissociationType, Modification loopMass,
            List<int> modPos, PeptideWithSetModifications peptide)
        {
            Dictionary<Tuple<int, int>, List<Product>> AllTheoreticalFragmentsLists = new Dictionary<Tuple<int, int>, List<Product>>();
            var originalFragments = new List<Product>();
            peptide.Fragment(dissociationType, FragmentationTerminus.Both, originalFragments);
            var loopProducts = new List<Product>();

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
                        .Where(p => p.Terminus == FragmentationTerminus.N && p.AminoAcidPosition < position1
                        || p.Terminus == FragmentationTerminus.C && p.AminoAcidPosition > position2).ToList();

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
                    
                    peptideWithLoop.Fragment(dissociationType, FragmentationTerminus.Both, loopProducts);
                    loopFragments.AddRange(loopProducts.Where(p => p.Terminus == FragmentationTerminus.N && p.AminoAcidPosition >= position2));

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

                    peptideWithLoop.Fragment(dissociationType, FragmentationTerminus.Both, loopProducts);
                    loopFragments.AddRange(
                        loopProducts.Where(p => p.Terminus == FragmentationTerminus.C && p.AminoAcidPosition <= position1));

                    AllTheoreticalFragmentsLists.Add(loopPositions, loopFragments);
                }
            }

            return AllTheoreticalFragmentsLists;
        }
    }
}
