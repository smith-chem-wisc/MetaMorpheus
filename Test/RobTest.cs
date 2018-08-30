using Chemistry;
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
            string[] proteinSequences = { "AB--------",   // 1: contains unique
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
            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch>();

            MsDataScan fakeScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);

            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(fakeScan, 2, 0, "File");

            foreach (var peptide in peptides)
            {
                switch (peptide.BaseSequence)
                {
                    case "A": psms.Add(new PeptideSpectralMatch(peptide, 0, 10, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "B": psms.Add(new PeptideSpectralMatch(peptide, 0, 9, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "C": psms.Add(new PeptideSpectralMatch(peptide, 0, 8, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "D": psms.Add(new PeptideSpectralMatch(peptide, 0, 7, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "E": psms.Add(new PeptideSpectralMatch(peptide, 0, 6, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "F": psms.Add(new PeptideSpectralMatch(peptide, 0, 5, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "G": psms.Add(new PeptideSpectralMatch(peptide, 0, 4, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "H": psms.Add(new PeptideSpectralMatch(peptide, 0, 3, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                    case "I": psms.Add(new PeptideSpectralMatch(peptide, 0, 2, 0, scan, digestionParams, new List<MatchedFragmentIon>())); break;
                }
            }
            
            foreach (var psm in psms)
            {
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
            List<string> parsimonyProteinSequences = psms.SelectMany(p => p.BestMatchingPeptideWithSetMods.Select(v => v.Pwsm.Protein)).Select(v => v.BaseSequence).Distinct().ToList();
            
            // check that correct proteins are in parsimony list
            Assert.Contains("AB--------", parsimonyProteinSequences);
            Assert.Contains("--C-------", parsimonyProteinSequences);
            Assert.Contains("-B-D---HHH--", parsimonyProteinSequences);
            Assert.Contains("-----F----*", parsimonyProteinSequences);
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
            Assert.AreEqual(4, proteinGroups.Count);
            Assert.AreEqual(1, proteinGroups.First().Proteins.Count);
            Assert.AreEqual("AB--------", proteinGroups.First().Proteins.First().BaseSequence);
            Assert.AreEqual(4, proteinGroups.First().AllPsmsBelowOnePercentFDR.Count);
            Assert.AreEqual(19, proteinGroups.First().ProteinGroupScore);
        }

        [Test]
        public static void TestFragments()
        {
            // creates some test proteins, digest, and fragment
            string[] sequences = { "GLSDGEWQQVLNVWGK" }; // just one peptide

            var peptides = new HashSet<PeptideWithSetModifications>();

            var p = new List<Protein>();
            for (int i = 0; i < sequences.Length; i++)
                p.Add(new Protein(sequences[i], (i + 1).ToString()));

            DigestionParams digestionParams = new DigestionParams();
            foreach (var protein in p)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                    peptides.Add(peptide);
            }

            Dictionary<DissociationType, List<double[]>> testDictionary = new Dictionary<DissociationType, List<double[]>>();

            foreach (var peptide in peptides)
            {
                CompactPeptide cp = new CompactPeptide(peptide, FragmentationTerminus.Both);
                foreach (DissociationType dt in Enum.GetValues(typeof(DissociationType)))
                {
                    if (testDictionary.ContainsKey(dt))
                    {
                        testDictionary[dt].Add(peptide.Fragment(dt, FragmentationTerminus.Both).Select(m => m.NeutralMass).ToArray());
                    }
                    else
                    {
                        testDictionary.Add(dt, new List<double[]> { peptide.Fragment(dt, FragmentationTerminus.Both).Select(m => m.NeutralMass).ToArray() });
                    }
                }
            }

            foreach (DissociationType dt in Enum.GetValues(typeof(DissociationType)))
            {
                switch (dt)
                {
                    case DissociationType.AnyActivationType:
                    case DissociationType.CID:
                    case DissociationType.Custom:
                    case DissociationType.ECD:
                    case DissociationType.ETD:
                    case DissociationType.EThCD:
                    case DissociationType.HCD:
                    case DissociationType.ISCID:
                    case DissociationType.MPD:
                    case DissociationType.PQD:
                    case DissociationType.Unknown:
                        Assert.IsTrue(true);
                        throw new Exception();
                    default:
                        break;
                }

            }

        }

        [Test]
        public static void TestNeutralLossFragments()
        {
            Protein p = new Protein("SR", "ac");

            ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);

            var cool = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var nice = cool.CompactPeptide(FragmentationTerminus.Both);
            Assert.AreEqual(2, nice.TerminalMasses.Where(t => t.Terminus == FragmentationTerminus.N).Count());
            Assert.AreEqual(1, nice.TerminalMasses.Where(t => t.Terminus == FragmentationTerminus.C).Count());
        }

        [Test]
        public static void TestPTMOutput()
        {
            List<Modification> variableModifications = new List<Modification>();
            List<Modification> fixedModifications = new List<Modification>();

            ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
            variableModifications.Add(new Modification(_id: "resMod", _modificationType: "HaHa", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H")));

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
            Assert.AreEqual("QQQ", value.First().Sequence);

            var firstProtDigest = proteinList.First().Digest(digestionParams, fixedModifications, variableModifications).ToList();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { firstProtDigest[0] };
            PeptideWithSetModifications compactPeptide2 = value2.First();
            Assert.AreEqual("MNNNSK", value2.First().Sequence);

            HashSet<PeptideWithSetModifications> value2mod = new HashSet<PeptideWithSetModifications> { firstProtDigest[1] };
            PeptideWithSetModifications compactPeptide2mod = value2mod.Last();
            Assert.AreEqual("MNNNS[HaHa:resMod]K", value2mod.Last().Sequence);

            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { firstProtDigest[2] };
            PeptideWithSetModifications compactPeptide3 = value3.First();
            Assert.AreEqual("NNNSK", value3.First().Sequence);
            HashSet<PeptideWithSetModifications> value3mod = new HashSet<PeptideWithSetModifications> { firstProtDigest[3] };

            PeptideWithSetModifications compactPeptide3mod = value3mod.Last();
            Assert.AreEqual("NNNS[HaHa:resMod]K", value3mod.Last().Sequence);

            var peptideList = new HashSet<PeptideWithSetModifications>();
            foreach (var protein in proteinList)
            {
                foreach (var peptide in protein.Digest(digestionParams, new List<Modification>(), variableModifications))
                {
                    peptideList.Add(peptide);
                }
            }
            
            MsDataScan jdfk = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass ms2scan = new Ms2ScanWithSpecificMass(jdfk, 2, 0, "File");

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
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

            ProteinParsimonyEngine engine = new ProteinParsimonyEngine(psms, true, new CommonParameters(), new List<string> { "ff" });
            var cool = (ProteinParsimonyResults)engine.Run();
            var proteinGroups = cool.ProteinGroups;

            ProteinScoringAndFdrEngine f = new ProteinScoringAndFdrEngine(proteinGroups, psms, false, false, true, new CommonParameters(), new List<string>());
            f.Run();

            Assert.AreEqual("#aa5[resMod,info:occupancy=0.67(2/3)];", proteinGroups.First().ModsInfo[0]);
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