﻿using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class RobTest
    {
        [Test]
        public static void TestParsimony()
        {
            // creates some proteins to test parsimony with
            string[] proteinSequences = {
                                   "AB--------",   // 1: contains unique
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

            var proteins = new List<Protein>();

            for (int i = 0; i < proteinSequences.Length; i++)
            {
                proteins.Add(new Protein(proteinSequences[i], (i + 1).ToString()));
            }
            proteins.Add(new Protein("-----F----*", "D1", isDecoy: true));
            proteins.Add(new Protein("-----F----**", "C1", isContaminant: true));
            proteins.Add(new Protein("----E----**", "C2", isContaminant: true));

            // create the protease
            IEnumerable<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage = new List<Tuple<string, FragmentationTerminus>>
            { new Tuple<string, FragmentationTerminus>("A", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("B", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("C", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("D", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("E", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("F", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("H", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("I", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("J", FragmentationTerminus.C),
                new Tuple<string, FragmentationTerminus>("-", FragmentationTerminus.C) };

            var protease = new Protease("test", sequencesInducingCleavage, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1);

            // digest the proteins
            var peptides = new HashSet<PeptideWithSetModifications>();
            foreach (Protein protein in proteins)
            {
                foreach (PeptideWithSetModifications peptide in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                {
                    switch (peptide.BaseSequence)
                    {
                        case "A": peptides.Add(peptide); break;
                        case "B": peptides.Add(peptide); break;
                        case "C": peptides.Add(peptide); break;
                        case "D": peptides.Add(peptide); break;
                        case "E": peptides.Add(peptide); break;
                        case "F": peptides.Add(peptide); break;
                        case "G": peptides.Add(peptide); break;
                        case "H": peptides.Add(peptide); break;
                        case "I": peptides.Add(peptide); break;
                    }
                }
            }

            // create PSMs for the peptides
            Dictionary<string, PeptideSpectralMatch> temp = new Dictionary<string, PeptideSpectralMatch>();

            MsDataScan fakeScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(fakeScan, 2, 0, "File", new CommonParameters());

            foreach (var peptide in peptides)
            {
                if (temp.TryGetValue(peptide.BaseSequence, out var psm))
                {
                    psm.AddOrReplace(peptide, 1, 0, true, new List<MatchedFragmentIon>());
                }
                else
                {
                    temp.Add(peptide.BaseSequence, new PeptideSpectralMatch(peptide, 0, 1, 0, scan, digestionParams, new List<MatchedFragmentIon>()));
                }
            }

            List<PeptideSpectralMatch> psms = temp.Values.ToList();

            foreach (var psm in psms)
            {
                psm.ResolveAllAmbiguities();
                psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            }

            // run parsimony
            ProteinParsimonyEngine parsimonyEngine = new ProteinParsimonyEngine(psms, false, new CommonParameters(), new List<string>());
            var parsimonyResults = (ProteinParsimonyResults)parsimonyEngine.Run();
            var proteinGroups = parsimonyResults.ProteinGroups;

            ProteinScoringAndFdrEngine proteinScoringAndFdrEngine = new ProteinScoringAndFdrEngine(proteinGroups, psms, true, false, true, new CommonParameters(), new List<string>());
            var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)proteinScoringAndFdrEngine.Run();
            proteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            // select the PSMs' proteins
            List<string> parsimonyProteinSequences = psms.SelectMany(p => p.BestMatchingPeptides.Select(v => v.Peptide.Protein)).Select(v => v.BaseSequence).Distinct().ToList();

            // check that correct proteins are in parsimony list
            Assert.Contains("AB--------", parsimonyProteinSequences);
            Assert.Contains("--C-------", parsimonyProteinSequences);
            Assert.Contains("-B-D---HHH--", parsimonyProteinSequences);
            Assert.Contains("----E----**", parsimonyProteinSequences);
            Assert.Contains("-B------I-", parsimonyProteinSequences);
            Assert.Contains("----EFG---", parsimonyProteinSequences);
            Assert.Contains("----EFG--J", parsimonyProteinSequences);
            Assert.AreEqual(8, parsimonyProteinSequences.Count);

            // sequence coverage test
            foreach (var proteinGroup in proteinGroups)
            {
                foreach (var coverage in proteinGroup.SequenceCoveragePercent)
                {
                    Assert.That(coverage <= 1.0);
                }
            }

            // test protein groups
            Assert.AreEqual(3, proteinGroups.Count);
            Assert.AreEqual(1, proteinGroups.First().Proteins.Count);
            Assert.AreEqual("AB--------", proteinGroups.First().Proteins.First().BaseSequence);
            Assert.AreEqual(2, proteinGroups.First().AllPsmsBelowOnePercentFDR.Count);
            Assert.AreEqual(2, proteinGroups.First().ProteinGroupScore);
        }

        [Test]
        public static void TestPTMOutput()
        {
            List<Modification> variableModifications = new List<Modification>();
            List<Modification> fixedModifications = new List<Modification>();

            ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
            variableModifications.Add(new Modification(_originalId: "resMod", _modificationType: "HaHa", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H")));

            var proteinList = new List<Protein> { new Protein("MNNNSKQQQ", "accession") };
            var protease = new Protease("CustomProtease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>
            {
                {variableModifications.Last(), 1 }
            };

            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1);

            var modPep = proteinList.First().Digest(digestionParams, fixedModifications, variableModifications).Last();
            HashSet<PeptideWithSetModifications> value = new HashSet<PeptideWithSetModifications> { modPep };
            PeptideWithSetModifications compactPeptide1 = value.First();
            Assert.AreEqual("QQQ", value.First().FullSequence);//this might be base

            var firstProtDigest = proteinList.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { firstProtDigest[0] };
            PeptideWithSetModifications compactPeptide2 = value2.First();
            Assert.AreEqual("MNNNSK", value2.First().FullSequence);//this might be base

            HashSet<PeptideWithSetModifications> value2mod = new HashSet<PeptideWithSetModifications> { firstProtDigest[1] };
            PeptideWithSetModifications compactPeptide2mod = value2mod.Last();
            Assert.AreEqual("MNNNS[HaHa:resMod on S]K", value2mod.Last().FullSequence);//this might be base

            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { firstProtDigest[2] };
            PeptideWithSetModifications compactPeptide3 = value3.First();
            Assert.AreEqual("NNNSK", value3.First().FullSequence);//this might be base
            HashSet<PeptideWithSetModifications> value3mod = new HashSet<PeptideWithSetModifications> { firstProtDigest[3] };

            PeptideWithSetModifications compactPeptide3mod = value3mod.Last();
            Assert.AreEqual("NNNS[HaHa:resMod on S]K", value3mod.Last().FullSequence);//this might be base

            var peptideList = new HashSet<PeptideWithSetModifications>();
            foreach (var protein in proteinList)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), variableModifications))
                {
                    peptideList.Add(peptide);
                }
            }

            MsDataScan jdfk = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass ms2scan = new Ms2ScanWithSpecificMass(jdfk, 2, 0, "File", new CommonParameters());
            
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            var match1 = new PeptideSpectralMatch(peptideList.ElementAt(0), 0, 10, 0, ms2scan, digestionParams, new List<MatchedFragmentIon>())
            {
            };
            match1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            var match2 = new PeptideSpectralMatch(peptideList.ElementAt(1), 0, 10, 0, ms2scan, digestionParams, new List<MatchedFragmentIon>())
            {
            };
            match2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            var match3 = new PeptideSpectralMatch(peptideList.ElementAt(1), 0, 10, 0, ms2scan, digestionParams, new List<MatchedFragmentIon>())
            {
            };
            match3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>
            {
                match1,
                match2,
                match3
            };

            psms.ForEach(p => p.ResolveAllAmbiguities());

            ProteinParsimonyEngine engine = new ProteinParsimonyEngine(psms, true, new CommonParameters(), new List<string> { "ff" });
            var cool = (ProteinParsimonyResults)engine.Run();
            var proteinGroups = cool.ProteinGroups;

            ProteinScoringAndFdrEngine f = new ProteinScoringAndFdrEngine(proteinGroups, psms, false, false, true, new CommonParameters(), new List<string>());
            f.Run();

            Assert.AreEqual("#aa5[resMod on S,info:occupancy=0.67(2/3)];", proteinGroups.First().ModsInfo[0]);
        }

        [Test]
        public static void TestProteinGroupsAccessionOutputOrder()
        {
            var p = new HashSet<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();

            // make protein B
            p.Add(new Protein("-----F----*", "B", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));

            // make protein A
            p.Add(new Protein("-----F----**", "A", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));

            // add protein B and A to the protein group
            ProteinGroup testGroup = new ProteinGroup(p, null, null);

            // test order is AB and not BA
            Assert.That(testGroup.ProteinGroupName.Equals("A|B"));
            Assert.That(testGroup.Proteins.First().Accession.Equals("B"));
        }
    }
}