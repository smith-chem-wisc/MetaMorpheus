using Chemistry;
using EngineLayer;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class GptmdEngineTest
    {
        [Test]
        [TestCase("NNNNN", "accession", @"not applied", 5)]
        [TestCase("NNNNN", "accession", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 4)]
        public static void TestGptmdEngine(string proteinSequence, string accession, string sequenceVariantDescription, int numModifiedResidues)
        {
            List<PeptideSpectralMatch> allResultingIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            allResultingIdentifications = new List<PeptideSpectralMatch>();
            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", "A", sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantProteins();

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(digestionParams, new List<Modification>(), variableModifications)).First();

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null), (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            allResultingIdentifications.Add(newPsm);

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(numModifiedResidues, res.Mods["accession"].Count);
        }

        [Test]
        [TestCase("NNNPPP", "accession", "A", @"not applied", 1, 6, 3, 3, 0)]
        [TestCase("NNNPPP", "accession", "A", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 1, 5, 2, 3, 0)]
        [TestCase("NNNPPP", "accession", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 2, 5, 2, 3, 1)]
        public static void TestCombos(string proteinSequence, string accession, string variantAA, string sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, int numModifiedResiduesNP)
        {
            List<PeptideSpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                                                                      new Modification(_originalId: "16",  _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", variantAA, sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantProteins();

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(digestionParams, new List<Modification>(), variableModifications)).First();

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943 + 15.994915).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            allIdentifications = new List<PeptideSpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(numModHashes, res.Mods.Count);
            Assert.AreEqual(numModifiedResidues, res.Mods["accession"].Count);
            Assert.AreEqual(numModifiedResiduesN, res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("21")).Count());
            Assert.AreEqual(numModifiedResiduesP, res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("16")).Count());
            res.Mods.TryGetValue("accession_N1P", out var hash);
            Assert.AreEqual(numModifiedResiduesNP, (hash ?? new HashSet<Tuple<int, Modification>>()).Count);
        }

        //[Test]
        //public static void LoadOriginalMismatchedModifications()
        //{
        //    var protein = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "oblm.xml"), true, DecoyType.Reverse, null, false, null, out var unknownModifications);
        //    Assert.AreEqual(0, protein[0].OneBasedPossibleLocalizedModifications.Count);
        //    protein[0].RestoreUnfilteredModifications();
        //    Assert.AreEqual(1, protein[0].OneBasedPossibleLocalizedModifications.Count);
        //}

        [Test]
        public static void TestSearchPtmVariantDatabase()
        {
            //Create Search Task
            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                },
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5))
            };

            //add task to task list
            var taskList = new List<(string, MetaMorpheusTask)> { ("task1", task1) };

            //create modification lists
            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where
                (b => task1.CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            //protein Creation (One with mod and one without)
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
            var variant = new SequenceVariation(3, "P", "K", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G|||||||||||||||||||\tGT:AD:DP\t1/1:30,30:30");
            Protein testProteinWithMod = new Protein("PEPTID", "accession1", sequenceVariations: new List<SequenceVariation> { variant });
            string variantAcc = VariantApplication.GetAccession(testProteinWithMod, new[] { variant });
            //First Write XML Database
            string xmlName = "oblm.xml";

            //Add Mod to list and write XML input database
            var modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(1, new Modification(_originalId: "acetyl on P", _modificationType: "type", _target: motifP, _monoisotopicMass: 42, _locationRestriction: "Anywhere.")),
            };
            var hashVar = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(3, new Modification(_originalId: "acetyl on K", _modificationType: "type", _target: motifK, _monoisotopicMass: 42, _locationRestriction: "Anywhere.")),
            };
            modList.Add(testProteinWithMod.Accession, hash);
            modList.Add(variantAcc, hashVar);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { testProteinWithMod }, xmlName);

            //now write MZML file
            var variantProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, null, false, null, out var unknownModifications);
            var variantProtein = variantProteins[0];
            var variantDecoy = variantProteins[1];
            Assert.AreEqual(0, unknownModifications.Count);

            Assert.AreEqual(2, variantProteins.Count); // target & decoy
            Assert.AreEqual(2, variantProteins[0].OneBasedPossibleLocalizedModifications.Count);
            List<int> foundResidueIndicies = variantProtein.OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            List<int> expectedResidueIndices = new List<int>() { 1, 3 };
            Assert.That(foundResidueIndicies, Is.EquivalentTo(expectedResidueIndices));
            Assert.AreEqual(2, variantDecoy.OneBasedPossibleLocalizedModifications.Count);
            foundResidueIndicies = variantDecoy.OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            expectedResidueIndices = new List<int>() { 4, 6 }; //originally modified residues are now at the end in the decoy
            Assert.That(foundResidueIndicies, Is.EquivalentTo(expectedResidueIndices));

            var thisOk = unknownModifications;//for debugging
            var commonParamsAtThisPoint = task1.CommonParameters.DigestionParams; //for debugging

            var digestedList = variantProteins[0].GetVariantProteins()[0].Digest(task1.CommonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(4, digestedList.Count);

            //Set Peptide with 1 mod at position 3
            PeptideWithSetModifications pepWithSetMods1 = digestedList[1];

            //Finally Write MZML file
            Assert.AreEqual("PEK[type:acetyl on K]TID", pepWithSetMods1.FullSequence);//this might be base sequence
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1 });
            string mzmlName = @"hello.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            //run!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName },
                new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
        }

        [Test]
        [TestCase("P", "PETID", "junk", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Unassigned.", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Anywhere.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "N-terminal.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "Peptide N-terminal.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "C-terminal.", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Peptide C-terminal.", 1, 5, 1, false)]
        [TestCase("E", "PETID", "Anywhere.", 2, 5, 2, true)]
        [TestCase("E", "PETID", "N-terminal.", 2, 5, 2, true)]
        [TestCase("E", "PETID", "Peptide N-terminal.", 2, 5, 2, false)]
        [TestCase("E", "PETID", "C-terminal.", 2, 5, 2, false)]
        [TestCase("E", "PETID", "Peptide C-terminal.", 2, 5, 2, false)]
        [TestCase("D", "PETID", "Anywhere.", 5, 5, 5, true)]
        [TestCase("D", "PETID", "N-terminal.", 5, 5, 5, false)]
        [TestCase("D", "PETID", "Peptide N-terminal.", 5, 5, 5, false)]
        [TestCase("D", "PETID", "C-terminal.", 5, 5, 5, true)]
        [TestCase("D", "PETID", "Peptide C-terminal.", 5, 5, 5, true)]
        public static void Test_GptmdEngineModFits(string targetAminoAcid, string proteinSequence, string locationRestriction, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex, bool result)
        {
            ModificationMotif.TryGetMotif(targetAminoAcid, out ModificationMotif motif);
            Modification attemptToLocalize = new Modification(null, null, null, null, _target: motif, _locationRestriction: locationRestriction, _chemicalFormula: null, _monoisotopicMass: 1, _databaseReference: null, _taxonomicRange: null, _keywords: null, _neutralLosses: null, _diagnosticIons: null, _fileOrigin: null);
            Dictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>();
            oneBasedModifications.Add(proteinOneBasedIndex, new List<Modification>() { attemptToLocalize });
            Protein protein = new Protein(proteinSequence, null, null, null, oneBasedModifications, null, null, null, false, false, null, null, null, null, null, null, "");

            Assert.AreEqual(result, GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
        }

        [Test]
        public static void TestUniProtGptmdConflict()
        {
            // this unit test checks to make sure GPTMD does not annotate mods at residues on 
            // proteins where the equivalent uniprot mod already exists

            Modification uniProtPhospho = GlobalVariables.AllModsKnown.First(p => p.ModificationType == "UniProt" && p.IdWithMotif.Contains("Phosphoserine"));
            Modification mmPhospho = GlobalVariables.AllModsKnown.First(p => p.ModificationType == "Common Biological" && p.IdWithMotif.Contains("Phosphorylation on S"));

            Protein protein = new Protein("PEPTIDESK", "test",
                oneBasedModifications: new Dictionary<int, List<Modification>>() { { 8, new List<Modification> { uniProtPhospho } } });

            PeptideWithSetModifications pep = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First(p => p.AllModsOneIsNterminus.Count == 0);

            // mod should not fit anywhere on the protein
            for (int i = 0; i < pep.Length; i++)
            {
                bool modFits = GptmdEngine.ModFits(mmPhospho, protein, i + 1, pep.Length, pep.OneBasedStartResidueInProtein + i);

                Assert.That(!modFits);
            }

            // the following code is just a control to make sure the phosphorylation actually does fit
            // at the given residue if the UniProt phosphorylation is not already present
            var someOtherSMod = GlobalVariables.AllModsKnown.Where(p => p.ModificationType == "Common Biological" && p.IdWithMotif.Contains("HexNAc on S")).First();

            protein = new Protein("PEPTIDESK", "test",
                oneBasedModifications: new Dictionary<int, List<Modification>>() { { 8, new List<Modification> { someOtherSMod } } });

            pep = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First(p => p.AllModsOneIsNterminus.Count == 0);

            // mod should fit at position 8
            for (int i = 0; i < pep.Length; i++)
            {
                bool modFits = GptmdEngine.ModFits(mmPhospho, protein, i + 1, pep.Length, pep.OneBasedStartResidueInProtein + i);

                if (i + 1 == 8)
                {
                    Assert.That(modFits);
                }
                else
                {
                    Assert.That(!modFits);
                }
            }
        }
    }
}