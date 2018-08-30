using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Text;
using EngineLayer;
using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using System.Linq;
using IO.MzML;
using Proteomics.ProteolyticDigestion;

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

            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("-", FragmentationTerminus.C), new Tuple<string, FragmentationTerminus>("Z", FragmentationTerminus.C) };
            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage2 = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) };

            var protease1 = new Protease("test1", sequencesInducingCleavage, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease1.Name, protease1);
            var protease2 = new Protease("test2", sequencesInducingCleavage2, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);

            //a hashset of peptideWithSetModifications from all proteases
            var pwsmList = new HashSet<PeptideWithSetModifications>();

            var proteinList = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                proteinList.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            DigestionParams digestionParams = new DigestionParams(protease: protease1.Name, minPeptideLength: 1);
            DigestionParams digestionParams2 = new DigestionParams(protease: protease2.Name, minPeptideLength: 1);

            foreach (var protein in proteinList)
            {
                //using protease #1
                foreach (var pwsm in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (pwsm.BaseSequence)
                    {
                        case "ABC": pwsmList.Add(pwsm); break;
                        case "EFGABC": pwsmList.Add(pwsm); break;
                        case "XYZ": pwsmList.Add(pwsm); break;
                    }
                }
                //using protease #2
                foreach (var pwsm in protein.Digest(digestionParams2, new List<Modification>(), new List<Modification>()))
                {
                    switch (pwsm.BaseSequence)
                    {
                        case "ABC": pwsmList.Add(pwsm); break;
                        case "-XYZ-EFG": pwsmList.Add(pwsm); break;
                    }
                }
            }
                        
            
            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            // dummy scan
            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            foreach (var peptide in pwsmList)
            {
                switch (peptide.BaseSequence)
                {
                    case "ABC":
                        if (peptide.DigestionParams == digestionParams)
                        {
                            psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                            break;
                        }
                        if (peptide.DigestionParams == digestionParams2)
                        {
                            psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>()));
                            break;
                        }
                        else { break; }

                    case "EFGABC": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;

                    case "XYZ": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;

                    case "-XYZ-EFG": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>())); break;
                }
            }

            List<ProductType> IonTypes = new List<ProductType>();
            ProductType Bions = ProductType.B;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(Bions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_id: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            Modification mod2 = new Modification(_id: "Oxidation of M", _modificationType: "Common Variable", _target: motif2, _locationRestriction: "Anyhwere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
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

            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("C", FragmentationTerminus.C) };
            List<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage2 = new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) };

            var protease = new Protease("testA", sequencesInducingCleavage, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var protease2 = new Protease("testB", sequencesInducingCleavage2, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
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
            //generates HashSet of PeptidesWithSetMods
            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
                foreach (var peptide in protein.Digest(digestionParams2, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {

                        case "ABC": peptideList.Add(peptide); break;
                        case "EFG": peptideList.Add(peptide); break;
                    }
                }
            }

           
            // builds psm list to match to peptides
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File");

            // creates psms for specific PWSM
            foreach (var PWSM in peptideList)
            {
                if (PWSM.DigestionParams == digestionParams)
                {
                    psms.Add(new PeptideSpectralMatch(PWSM, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                }
                if (PWSM.DigestionParams == digestionParams2)
                {
                    psms.Add(new PeptideSpectralMatch(PWSM, 0, 10, 0, scan, digestionParams2, new List<MatchedFragmentIon>()));
                }
            }           
            

            List<ProductType> IonTypes = new List<ProductType>();
            ProductType Bions = ProductType.B;
            ProductType Yions = ProductType.Y;
            IonTypes.Add(Bions);
            IonTypes.Add(Yions);

            HashSet<DigestionParams> digestionParamsList = new HashSet<DigestionParams>();
            digestionParamsList.Add(digestionParams);
            digestionParamsList.Add(digestionParams2);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod = new Modification(_id: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            List<Modification> modVarList = new List<Modification> { mod };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif2);
            List<Modification> modFixedList = new List<Modification> { mod };

            ProteinParsimonyEngine ppe = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null);
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