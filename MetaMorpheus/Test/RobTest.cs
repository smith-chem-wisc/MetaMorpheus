﻿using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;
using Omics;
using Transcriptomics;

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
            List<DigestionMotif> digestionMotifs = new List<DigestionMotif>
            {
                new DigestionMotif("A", null, 1, null),
                new DigestionMotif("B", null, 1, null),
                new DigestionMotif("C", null, 1, null),
                new DigestionMotif("D", null, 1, null),
                new DigestionMotif("E", null, 1, null),
                new DigestionMotif("F", null, 1, null),
                new DigestionMotif("G", null, 1, null),
                new DigestionMotif("H", null, 1, null),
                new DigestionMotif("I", null, 1, null),
                new DigestionMotif("J", null, 1, null),
                new DigestionMotif("-", null, 1, null),
            };

            var protease = new Protease("test", CleavageSpecificity.Full, null, null, digestionMotifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1));

            // digest the proteins
            var peptides = new HashSet<PeptideWithSetModifications>();
            foreach (Protein protein in proteins)
            {
                foreach (PeptideWithSetModifications peptide in protein.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()))
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
            Dictionary<string, SpectralMatch> temp = new Dictionary<string, SpectralMatch>();

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
                    temp.Add(peptide.BaseSequence, new PeptideSpectralMatch(peptide, 0, 1, 0, scan, commonParameters, new List<MatchedFragmentIon>()));
                }
            }

            List<SpectralMatch> psms = temp.Values.ToList();

            foreach (var psm in psms)
            {
                psm.ResolveAllAmbiguities();
                psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            }

            // run parsimony
            ProteinParsimonyEngine parsimonyEngine = new ProteinParsimonyEngine(psms, false, new CommonParameters(), null, new List<string>());
            var parsimonyResults = (ProteinParsimonyResults)parsimonyEngine.Run();
            var proteinGroups = parsimonyResults.ProteinGroups;

            ProteinScoringAndFdrEngine proteinScoringAndFdrEngine = new ProteinScoringAndFdrEngine(proteinGroups, psms, true, false, true, new CommonParameters(), null, new List<string>());
            var proteinScoringAndFdrResults = (ProteinScoringAndFdrResults)proteinScoringAndFdrEngine.Run();
            proteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;

            // select the PSMs' proteins
            List<string> parsimonyProteinSequences = psms.SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer.Parent)).Select(v => v.BaseSequence).Distinct().ToList();

            // check that correct proteins are in parsimony list
            Assert.That(parsimonyProteinSequences.Contains("AB--------"));
            Assert.That(parsimonyProteinSequences.Contains("--C-------"));
            Assert.That(parsimonyProteinSequences.Contains("-B-D---HHH--"));
            Assert.That(parsimonyProteinSequences.Contains("-B------I-"));
            Assert.That(parsimonyProteinSequences.Contains("----EFG---"));
            Assert.That(parsimonyProteinSequences.Contains("----EFG--J"));
            Assert.That(parsimonyProteinSequences.Count, Is.EqualTo(6));

            // sequence coverage test
            foreach (var proteinGroup in proteinGroups)
            {
                foreach (var coverage in proteinGroup.SequenceCoverageFraction)
                {
                    Assert.That(coverage <= 1.0);
                }
            }

            // test protein groups
            Assert.That(proteinGroups.Count, Is.EqualTo(4));
            Assert.That(proteinGroups.First().Proteins.Count, Is.EqualTo(1));
            Assert.That(proteinGroups.First().Proteins.First().BaseSequence, Is.EqualTo("AB--------"));
            Assert.That(proteinGroups.First().AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(2));
            Assert.That(proteinGroups.First().ProteinGroupScore, Is.EqualTo(2));
        }

        [Test]
        public static void TestPTMOutput()
        {
            List<Modification> variableModifications = new List<Modification>();
            List<Modification> fixedModifications = new List<Modification>();

            ModificationMotif.TryGetMotif("S", out ModificationMotif motifS);
            ModificationMotif.TryGetMotif("I", out ModificationMotif motifI);
            variableModifications.Add(new Modification(_originalId: "resMod", _modificationType: "HaHa", _target: motifS, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H")));
            variableModifications.Add(new Modification(_originalId: "iModOne", _modificationType: "HaHa", _target: motifI, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H")));
            variableModifications.Add(new Modification(_originalId: "iModTwo", _modificationType: "HaHa", _target: motifI, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H")));

            var proteinList = new List<Protein> { new Protein("MNNNSKQQQI", "accession") };
            var protease = new Protease("CustomProtease", CleavageSpecificity.Full, null, null, new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) });
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>
            {
                {variableModifications.Last(), 1 }
            };

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1));
            var protDigest = proteinList.First().Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            int idx = 0;

            var pep1 = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep1.Single().FullSequence, Is.EqualTo("MNNNSK"));//this might be base

            var pep1mod = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep1mod.Single().FullSequence, Is.EqualTo("MNNNS[HaHa:resMod on S]K"));//this might be base

            var pep3 = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep3.Single().FullSequence, Is.EqualTo("NNNSK"));//this might be base

            var pep3mod = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep3mod.Single().FullSequence, Is.EqualTo("NNNS[HaHa:resMod on S]K"));//this might be base

            var pep4 = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep4.Single().FullSequence, Is.EqualTo("QQQI"));//this might be base

            var pep4mod1 = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep4mod1.Single().FullSequence, Is.EqualTo("QQQI[HaHa:iModOne on I]"));//this might be base

            var pep4mod2 = new HashSet<IBioPolymerWithSetMods> { protDigest[idx++] };
            Assert.That(pep4mod2.Single().FullSequence, Is.EqualTo("QQQI[HaHa:iModTwo on I]"));//this might be base

            var peptideList = new HashSet<IBioPolymerWithSetMods>();
            foreach (var peptide in proteinList.SelectMany(protein => protein.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications)))
            {               
                peptideList.Add(peptide);
            }

            MsDataScan jdfk = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass ms2scan = new Ms2ScanWithSpecificMass(jdfk, 2, 0, "File", new CommonParameters());
            
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            var match1 = new PeptideSpectralMatch(peptideList.ElementAt(0), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()){};
            match1.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match2 = new PeptideSpectralMatch(peptideList.ElementAt(1), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()){};
            match2.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match3 = new PeptideSpectralMatch(peptideList.ElementAt(1), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()){};
            match3.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match4 = new PeptideSpectralMatch(peptideList.ElementAt(4), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()) { };
            match4.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match5 = new PeptideSpectralMatch(peptideList.ElementAt(5), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()) { };
            match5.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match6 = new PeptideSpectralMatch(peptideList.ElementAt(6), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()) { };
            match6.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match44 = new PeptideSpectralMatch(peptideList.ElementAt(4), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()) { };
            match44.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match55 = new PeptideSpectralMatch(peptideList.ElementAt(5), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()) { };
            match55.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);
            var match66 = new PeptideSpectralMatch(peptideList.ElementAt(6), 0, 10, 0, ms2scan, commonParameters, new List<MatchedFragmentIon>()) { };
            match66.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);


            List<SpectralMatch> psms = new List<SpectralMatch>
            {
                match1,
                match2,
                match3,
                match4,
                match44,
                match5,
                match55,
                match6,
                match66
            };

            psms.ForEach(p => p.ResolveAllAmbiguities());

            ProteinParsimonyEngine engine = new ProteinParsimonyEngine(psms, true, new CommonParameters(), null, new List<string> { "ff" });
            var cool = (ProteinParsimonyResults)engine.Run();
            var proteinGroups = cool.ProteinGroups;

            ProteinScoringAndFdrEngine f = new ProteinScoringAndFdrEngine(proteinGroups, psms, false, false, true, new CommonParameters(), null, new List<string>());
            f.Run();

            Assert.That(proteinGroups.First().ModsInfo[0], Is.EqualTo("#aa5[resMod on S,info:occupancy=0.67(2/3)];#aa10[iModOne on I,info:occupancy=0.33(2/6)];#aa10[iModTwo on I,info:occupancy=0.33(2/6)]"));
        }

        [Test]
        public static void TestProteinGroupsAccessionOutputOrder()
        {
            var p = new HashSet<IBioPolymer>();
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

        [Test]
        public static void TestTranscriptGroupsAccessionOutputOrder()
        {
            var p = new HashSet<IBioPolymer>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();

            // make protein B
            p.Add(new RNA("AAAAACAAAAU", "B", isDecoy: true));

            // make protein A
            p.Add(new RNA("AAAAACAAAAUU", "A", isDecoy: true));

            // add protein B and A to the protein group
            ProteinGroup testGroup = new ProteinGroup(p, null, null);

            // test order is AB and not BA
            Assert.That(testGroup.ProteinGroupName.Equals("A|B"));
            Assert.That(testGroup.Proteins.First().Accession.Equals("B"));
        }
    }
}