using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class MultiProteaseParsimonyTest
    {
        [Test]
        public static void MultiProteaseUniqueTest()
        {
            string[] sequences = {
                "-XYZ--ABC",
                "-XYZ-EFGABC",
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("-", TerminusType.C), new Tuple<string, TerminusType>("Z", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("test1", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("test2", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFGABC": peptideList.Add(peptide); break;
                        case "XYZ": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "-XYZ-EFG": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[2] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[5] });
            dictionary.Add(peptides[3], new HashSet<PeptideWithSetModifications> { PWSM[3] });
            dictionary.Add(peptides[4], new HashSet<PeptideWithSetModifications> { PWSM[4] });

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
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }

                        case "EFGABC": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "XYZ": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams)); break;

                        case "-XYZ-EFG": psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2)); break;
                    }
                }
            }

            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Oxidation of M", "Common Variable", motif2, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            Assert.AreEqual(4, CompactPeptidesToProteinPeptidesMatching.Count);

            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.Count);
            Assert.AreEqual("XYZ", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("XYZ", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(1).BaseSequence);

            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.Count);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(1).BaseSequence);

            Assert.AreEqual(1, CompactPeptidesToProteinPeptidesMatching.ElementAt(2).Value.Count);
            Assert.AreEqual("EFGABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(2).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("2", CompactPeptidesToProteinPeptidesMatching.ElementAt(2).Value.ElementAt(0).Protein.Accession);

            Assert.AreEqual(1, CompactPeptidesToProteinPeptidesMatching.ElementAt(3).Value.Count);
            Assert.AreEqual("-XYZ-EFG", CompactPeptidesToProteinPeptidesMatching.ElementAt(3).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("2", CompactPeptidesToProteinPeptidesMatching.ElementAt(3).Value.ElementAt(0).Protein.Accession);

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);

            if (proteinGroups.ElementAt(0).ProteinGroupName == "1")
            {
                Assert.AreEqual(2, proteinGroups.ElementAt(0).AllPeptides.Count);
                Assert.AreEqual(1, proteinGroups.ElementAt(0).UniquePeptides.Count);
                Assert.AreEqual("XYZ", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).BaseSequence);

                Assert.AreEqual(4, proteinGroups.ElementAt(1).AllPeptides.Count);
                Assert.AreEqual(3, proteinGroups.ElementAt(1).UniquePeptides.Count);
                Assert.AreEqual("XYZ", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).BaseSequence);
                Assert.AreEqual("test2", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
                Assert.AreEqual("EFGABC", proteinGroups.ElementAt(1).AllPeptides.ElementAt(2).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(1).AllPeptides.ElementAt(2).DigestionParams.Protease.Name);
                Assert.AreEqual("-XYZ-EFG", proteinGroups.ElementAt(1).AllPeptides.ElementAt(3).BaseSequence);
                Assert.AreEqual("test2", proteinGroups.ElementAt(1).AllPeptides.ElementAt(3).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).BaseSequence);
                Assert.AreEqual("EFGABC", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(1).BaseSequence);
                Assert.AreEqual("-XYZ-EFG", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(2).BaseSequence);
            }
            else
            {
                Assert.AreEqual(2, proteinGroups.ElementAt(1).AllPeptides.Count);
                Assert.AreEqual(1, proteinGroups.ElementAt(1).UniquePeptides.Count);
                Assert.AreEqual("XYZ", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).BaseSequence);

                Assert.AreEqual(4, proteinGroups.ElementAt(0).AllPeptides.Count);
                Assert.AreEqual(3, proteinGroups.ElementAt(0).UniquePeptides.Count);
                Assert.AreEqual("XYZ", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).BaseSequence);
                Assert.AreEqual("test2", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
                Assert.AreEqual("EFGABC", proteinGroups.ElementAt(0).AllPeptides.ElementAt(2).BaseSequence);
                Assert.AreEqual("test1", proteinGroups.ElementAt(0).AllPeptides.ElementAt(2).DigestionParams.Protease.Name);
                Assert.AreEqual("-XYZ-EFG", proteinGroups.ElementAt(0).AllPeptides.ElementAt(3).BaseSequence);
                Assert.AreEqual("test2", proteinGroups.ElementAt(0).AllPeptides.ElementAt(3).DigestionParams.Protease.Name);
                Assert.AreEqual("ABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).BaseSequence);
                Assert.AreEqual("EFGABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(1).BaseSequence);
                Assert.AreEqual("-XYZ-EFG", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(2).BaseSequence);
            }
        }

        [Test]
        public static void MultiProteaseIndistiguishableTest()
        {
            string[] sequences = {
                "ABCEFG",
                "EFGABC",
            };

            List<Tuple<string, TerminusType>> sequencesInducingCleavage = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("C", TerminusType.C) };
            List<Tuple<string, TerminusType>> sequencesInducingCleavage2 = new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) };

            var protease = new Protease("testA", sequencesInducingCleavage, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testB", sequencesInducingCleavage2, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            var peptideList = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<ModificationWithMass>(), new List<ModificationWithMass>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
            }

            // creates the initial dictionary of "peptide" and "virtual peptide" matches
            var dictionary = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            CompactPeptide[] peptides = new CompactPeptide[peptideList.Count];

            PeptideWithSetModifications[] PWSM = new PeptideWithSetModifications[peptideList.Count];

            // creates peptide list
            for (int i = 0; i < peptideList.Count; i++)
            {
                peptides[i] = new CompactPeptide(peptideList.ElementAt(i), TerminusType.None);
                PWSM[i] = peptideList.ElementAt(i);
            }

            dictionary.Add(peptides[0], new HashSet<PeptideWithSetModifications> { PWSM[0], PWSM[3] });
            dictionary.Add(peptides[1], new HashSet<PeptideWithSetModifications> { PWSM[1], PWSM[2] });

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
                        case "ABC":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }

                        case "EFG":
                            if (peptide.DigestionParams == digestionParams)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams));
                                break;
                            }
                            if (peptide.DigestionParams == digestionParams2)
                            {
                                psms.Add(new PeptideSpectralMatch(peptide.CompactPeptide(TerminusType.None), 0, 10, 0, scan, digestionParams2));
                                break;
                            }
                            else { break; }
                    }
                }
            }

            List<ProductType> IonTypes = new List<ProductType>();
            ProductType BnoB1ions = ProductType.BnoB1ions;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(BnoB1ions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            List<ModificationWithMass> modVarList = new List<ModificationWithMass> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            List<ModificationWithMass> modFixedList = new List<ModificationWithMass> { mod };
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
            new SequencesToActualProteinPeptidesEngine(psms, p, modFixedList, modVarList, IonTypes, digestionParamsList, true, new CommonParameters(), null);
            var results = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var CompactPeptidesToProteinPeptidesMatching = results.CompactPeptideToProteinPeptideMatching;
            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.Count);

            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.Count);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("ABC", CompactPeptidesToProteinPeptidesMatching.ElementAt(0).Value.ElementAt(1).BaseSequence);

            Assert.AreEqual(2, CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.Count);
            Assert.AreEqual("EFG", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(0).BaseSequence);
            Assert.AreEqual("EFG", CompactPeptidesToProteinPeptidesMatching.ElementAt(1).Value.ElementAt(1).BaseSequence);

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(CompactPeptidesToProteinPeptidesMatching, false, new CommonParameters(), null);
            var proteinAnalysisResults = (ProteinParsimonyResults)ppe.Run();

            List<ProteinGroup> proteinGroups = proteinAnalysisResults.ProteinGroups;
            Assert.AreEqual(2, proteinGroups.Count);

            Assert.AreEqual(2, proteinGroups.ElementAt(0).AllPeptides.Count);
            Assert.AreEqual(2, proteinGroups.ElementAt(0).UniquePeptides.Count);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("testA", proteinGroups.ElementAt(0).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).BaseSequence);
            Assert.AreEqual("testA", proteinGroups.ElementAt(0).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(0).UniquePeptides.ElementAt(1).BaseSequence);

            Assert.AreEqual(2, proteinGroups.ElementAt(1).AllPeptides.Count);
            Assert.AreEqual(2, proteinGroups.ElementAt(1).UniquePeptides.Count);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("testB", proteinGroups.ElementAt(1).AllPeptides.ElementAt(0).DigestionParams.Protease.Name);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).BaseSequence);
            Assert.AreEqual("testB", proteinGroups.ElementAt(1).AllPeptides.ElementAt(1).DigestionParams.Protease.Name);
            Assert.AreEqual("ABC", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(0).BaseSequence);
            Assert.AreEqual("EFG", proteinGroups.ElementAt(1).UniquePeptides.ElementAt(1).BaseSequence);
        }
    }
}