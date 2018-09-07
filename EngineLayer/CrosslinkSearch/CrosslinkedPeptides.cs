using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkedPeptide
    {
        public static Dictionary<List<int>, List<Product>> XlGetTheoreticalFragmentIons(DissociationType dissociationType, bool Charge_2_3, Crosslinker crosslinker, List<int> possiblePositionsOfMod, double modMass, PeptideWithSetModifications peptide)
        {
            Dictionary<List<int>, List<Product>> AllTheoreticalFragmentsLists = new Dictionary<List<int>, List<Product>>();

            List<Product> baseTheoreticalFragmentIons = peptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();

            foreach (var crosslinkerPosition in possiblePositionsOfMod)
            {
                List<Product> currentProducts = new List<Product>();

                // add signature ions (molecular ion plus short or long piece of the crosslinker)
                if (crosslinker.Cleavable)
                {
                    currentProducts.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.None, peptide.MonoisotopicMass + crosslinker.CleaveMassShort, peptide.Length, peptide.Length), 0));
                    currentProducts.Add(new Product(ProductType.M, new NeutralTerminusFragment(FragmentationTerminus.None, peptide.MonoisotopicMass + crosslinker.CleaveMassLong, peptide.Length, peptide.Length), 0));
                }

                foreach (var product in baseTheoreticalFragmentIons)
                {
                    var productType = product.ProductType;

                    //TODO: use dissociation type and not ProductType to figure out product masses because UVPD could generate C ions and fragment the cleavable crosslinker, etc
                    switch (productType)
                    {
                        case ProductType.b:
                            if (product.TerminusFragment.AminoAcidPosition < crosslinkerPosition)
                            {
                                currentProducts.Add(product);
                            }
                            else
                            {
                                if (crosslinker.Cleavable)
                                {
                                    // cleavable crosslinker fragment mass; the mass is the fragment mass plus either the short or long piece of the fragmented crosslinker
                                    currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.CleaveMassShort, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                    currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.CleaveMassLong, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                }
                                else
                                {
                                    // non cleavable crosslinker mass; modMass is the mass of the thing it's crosslinked to (e.g., the other peptide, RNA, etc)
                                    currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + modMass + crosslinker.TotalMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                }
                            }
                            break;
                        case ProductType.c:
                            if (product.TerminusFragment.AminoAcidPosition < crosslinkerPosition)
                            {
                                currentProducts.Add(product);
                            }
                            else
                            {
                                currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + modMass + crosslinker.TotalMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                            }
                            break;
                        case ProductType.y:
                            if (product.TerminusFragment.AminoAcidPosition < crosslinkerPosition)
                            {
                                currentProducts.Add(product);
                            }
                            else
                            {
                                if (crosslinker.Cleavable)
                                {
                                    currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.CleaveMassShort, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                    currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.CleaveMassLong, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                }
                                else
                                {
                                    currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + modMass + crosslinker.TotalMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                }
                            }
                            break;
                        case ProductType.zPlusOne:
                            if (product.TerminusFragment.AminoAcidPosition < crosslinkerPosition)
                            {
                                currentProducts.Add(product);
                            }
                            else
                            {
                                currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + modMass + crosslinker.TotalMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                            }
                            break;
                        default:
                            break;
                    }
                }

                AllTheoreticalFragmentsLists.Add(new List<int> { crosslinkerPosition }, currentProducts);
            }

            return AllTheoreticalFragmentsLists;
        }

        //TO DO: the second ModPostion jPos is not recorded. 
        public static Dictionary<List<int>, List<Product>> XlLoopGetTheoreticalFragmentIons(DissociationType dissociationType, bool Charge_2_3, Crosslinker crosslinker, List<int> modPos, double modMass, PeptideWithSetModifications peptide)
        {
            Dictionary<List<int>, List<Product>> AllTheoreticalFragmentsLists = new Dictionary<List<int>, List<Product>>();

            List<Product> baseTheoreticalFragmentIons = peptide.Fragment(dissociationType, FragmentationTerminus.Both).ToList();

            if (modPos.Count >= 2)
            {
                for (int iPos = 0; iPos < modPos.Count - 1; iPos++)
                {
                    for (int jPos = iPos + 1; jPos < modPos.Count; jPos++)
                    {
                        List<Product> currentProducts = new List<Product>();

                        foreach (var product in baseTheoreticalFragmentIons)
                        {
                            var iType = product.ProductType;
                            switch (iType)
                            {
                                case ProductType.b:
                                    if (product.TerminusFragment.FragmentNumber < modPos[iPos])
                                    {
                                        currentProducts.Add(product);
                                    }
                                    else if (product.TerminusFragment.FragmentNumber >= modPos[jPos])
                                    {
                                        currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.LoopMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                    }
                                    break;
                                case ProductType.c:
                                    if (product.TerminusFragment.FragmentNumber < modPos[iPos])
                                    {
                                        currentProducts.Add(product);
                                    }
                                    else if (product.TerminusFragment.FragmentNumber >= modPos[jPos])
                                    {
                                        currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.LoopMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                    }
                                    break;
                                case ProductType.y:
                                    if (product.TerminusFragment.FragmentNumber < peptide.Length - modPos[jPos] + 1)
                                    {
                                        currentProducts.Add(product);
                                    }
                                    else if (product.TerminusFragment.FragmentNumber >= peptide.Length - modPos[iPos] + 1)
                                    {
                                        currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.LoopMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));
                                    }
                                    break;
                                case ProductType.zPlusOne:
                                    if (product.TerminusFragment.FragmentNumber < peptide.Length - modPos[jPos] + 1)
                                    {
                                        currentProducts.Add(product);
                                    }
                                    else if (product.TerminusFragment.FragmentNumber >= peptide.Length - modPos[iPos] + 1)
                                    {
                                        currentProducts.Add(new Product(product.ProductType, new NeutralTerminusFragment(product.TerminusFragment.Terminus, product.NeutralMass + crosslinker.LoopMass, product.TerminusFragment.FragmentNumber, product.TerminusFragment.AminoAcidPosition), 0));                                        
                                    }
                                    break;
                                default:
                                    break;
                            }
                        }

                        AllTheoreticalFragmentsLists.Add(new List<int> { modPos[iPos], modPos[jPos] }, currentProducts);
                    }
                }
            }

            return AllTheoreticalFragmentsLists;
        }
    }
}
