using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Test
{
    [TestFixture]
    class MultiProteaseParsimonyTest
    {
        [Test]
        public static void MultiProteaseCompactPeptideMatchingTest()
        {
            string[] sequences = { "AB--------",   // 1: contains unique
                                   "--C-------",   // 2: one hit wonder
                                   "---D---HHH--", // 3: subset
                                   "-B-D---HHH--", // 4: D should go to 4, not 3 (3 is subset)
                                   "-B--E-----",   // 5: subsumable
                                   "----EFG---",   // 6: indistinguishable from 8 (J will not be a "detected" PSM)
                                   "-----F----",   // 7: lone pep shared w/ decoy
                                   "--------I-",   // 8: I should go to 9, not 8
                                   "-B------I-",   // 9: I should go to 9, not 8
                                   "----EFG--J"    // 10: indistinguishable from 6 (J will not be a "detected" PSM)
                                   };

            IEnumerable<string> sequencesInducingCleavage = new List<string> { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "-" };
            var protease = new Protease("test1", sequencesInducingCleavage, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var protease2 = new Protease("test2", sequencesInducingCleavage, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var protease3 = new Protease("test3", sequencesInducingCleavage, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var peptideList = new HashSet<PeptideWithSetModifications>();
            GlobalVariables.ProteaseDictionary.Add(protease.Name, protease);
            GlobalVariables.ProteaseDictionary.Add(protease2.Name, protease2);
            GlobalVariables.ProteaseDictionary.Add(protease3.Name, protease3);
            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, MinPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, MinPeptideLength: 1);
            DigestionParams digestionParams3 = new DigestionParams(protease: protease3.Name, MinPeptideLength: 1);

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "A": peptideList.Add(peptide); break;
                        case "B": peptideList.Add(peptide); break;
                        case "C": peptideList.Add(peptide); break;
                        case "D": peptideList.Add(peptide); break;
                        case "E": peptideList.Add(peptide); break;
                        case "F": peptideList.Add(peptide); break;
                        case "G": peptideList.Add(peptide); break;
                        case "H": peptideList.Add(peptide); break;
                        case "I": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];
            HashSet<PeptideWithSetModifications>[] virtualPeptideSets = new HashSet<PeptideWithSetModifications>[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
            }

            // creates protein list
            for (int i = 0; i < virtualPeptideSets.Length; i++)
            {
                virtualPeptideSets[i] = new HashSet<PeptideWithSetModifications>();

                foreach (var virtualPeptide in peptideList)
                {
                    string peptideBaseSequence = string.Join("", peptideList.ElementAt(i).BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    if (virtualPeptide.BaseSequence.Contains(peptideBaseSequence))
                    {
                        virtualPeptideSets[i].Add(virtualPeptide);
                    }
                }
            }

            // populates initial peptide-virtualpeptide dictionary
            for (int i = 0; i < peptides.Length; i++)
            {
                if (!dictionary.ContainsKey(peptides[i]))
                {
                    dictionary.Add(peptides[i], virtualPeptideSets[i]);
                }
            }


            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var kvp in dictionary)
            {
                foreach (var peptide in kvp.Value)
                {
                    switch (peptide.BaseSequence)
                    {
                        case "A": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;
                        case "B": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;
                        case "C": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams3)); break;
                        case "D": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams3)); break;
                        case "E": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));break;
                        case "F": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams3)); break;
                        case "G": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;
                        case "H": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;
                        case "I": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams3)); break;
                    }
                }
            }

            
            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            List<DigestionParams> digestionParamsList = new List<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            digestionParamsList.Add(digestionParams3);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass>{ mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
                new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, false, null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var proteaseSortedCompactPeptidesToProteinPeptidesMatching = results.ProteaseSortedCompactPeptideToProteinPeptideMatching;

            Assert.AreEqual(3, proteaseSortedCompactPeptidesToProteinPeptidesMatching.Count);
            Assert.AreEqual(2, proteaseSortedCompactPeptidesToProteinPeptidesMatching[protease].Count);
            Assert.AreEqual(3, proteaseSortedCompactPeptidesToProteinPeptidesMatching[protease2].Count);
            Assert.AreEqual(4, proteaseSortedCompactPeptidesToProteinPeptidesMatching[protease3].Count);
        }
        
        [Test]
        public static void MultiParsimonyEngineTest()
        {

            string[] sequences = { "AB--------",   // 1: contains unique
                                   "--C-------",   // 2: one hit wonder
                                   "---D---HHH--", // 3: subset
                                   "-B-D---HHH--", // 4: D should go to 4, not 3 (3 is subset)
                                   "-B--E-----",   // 5: subsumable
                                   "----EFG---",   // 6: indistinguishable from 8 (J will not be a "detected" PSM)
                                   "-----F----",   // 7: lone pep shared w/ decoy
                                   "--------I-",   // 8: I should go to 9, not 8
                                   "-B------I-",   // 9: I should go to 9, not 8
                                   "----EFG--J"    // 10: indistinguishable from 6 (J will not be a "detected" PSM)
                                   };

            var proteaseSortedCompactPeptidesToProteinPeptidesMatching = new Dictionary<Protease, Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>> ();
            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(proteaseSortedCompactPeptidesToProteinPeptidesMatching, false, null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();
            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;
        }
    }
}
 