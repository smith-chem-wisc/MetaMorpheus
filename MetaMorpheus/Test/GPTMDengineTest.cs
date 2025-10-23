using Chemistry;
using EngineLayer;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;
using TaskLayer;
using UsefulProteomicsDatabases;
using EngineLayer.FdrAnalysis;
using EngineLayer.ClassicSearch;
using System.IO;
using System.Reflection;
using EngineLayer.DatabaseLoading;
using Omics.BioPolymer;
using Omics;
using Readers;

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
            List<SpectralMatch> allResultingIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            allResultingIdentifications = new List<SpectralMatch>();

            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", new CommonParameters()));

            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), fsp, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(0));

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", "A", sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantBioPolymers();
            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));

            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications)).First();

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null), (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<IBioPolymerWithSetMods> { modPep };
            SpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            newPsm.SetMs2Scan(scan.TheScan);
            allResultingIdentifications.Add(newPsm);

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null);
            res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(1));
            Assert.That(res.Mods["accession"].Count, Is.EqualTo(numModifiedResidues));
        }

        [Test]
        public static void TestGptmdEngineDissociationTypeAutodetect()
        {
            string origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            //var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters cp = new CommonParameters(maxThreadsToUsePerFile: 1, digestionParams: new DigestionParams());
            var commonParameters = cp.CloneWithNewDissociationType(DissociationType.Autodetect);
            SearchParameters SearchParameters = new SearchParameters();
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("SmallCalibratible_Yeast.mzML", commonParameters));
            Tolerance precursorMassTolerance = new PpmTolerance(20);
            var myMsDataFile = myFileManager.LoadFile(origDataFile, commonParameters);
            List<double> acceptableMassShifts = new List<double> { 0.984015583, 0.984015583 };
            MassDiffAcceptor searchModes = new DotMassDiffAcceptor("", acceptableMassShifts, precursorMassTolerance);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(myDatabase, true, DecoyType.Reverse, false, out var dbErrors, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, -1);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, @"TestData\SmallCalibratible_Yeast.mzML", commonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, fsp, null, new List<string>(), SearchParameters.WriteSpectralLibrary).Run();
            FdrAnalysisResults fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.Where(p => p != null).ToList(), 1,
                commonParameters, fsp, new List<string>()).Run());

            var nonNullPsms = allPsmsArray.Where(p => p != null).ToList();
            foreach(var psm in nonNullPsms)
            {
                psm.SetMs2Scan(listOfSortedms2Scans[psm.ScanIndex].TheScan);
            }
            GptmdParameters g = new GptmdParameters();
            List<Modification> gptmdModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => g.ListOfModsGptmd.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var reducedMods = new List<Modification>();
            foreach (var mod in gptmdModifications)
            {
                if (mod.IdWithMotif == "Deamidation on N" || mod.IdWithMotif == "Citrullination on R")
                {
                    reducedMods.Add(mod);
                }
            }
            
            var engine = new GptmdEngine(nonNullPsms, reducedMods, new List<Tuple<double, double>>(), new Dictionary<string, Tolerance> { { @"TestData\SmallCalibratible_Yeast.mzML", precursorMassTolerance } }, commonParameters, fsp, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.That(8, Is.EqualTo(res.Mods.Count));
        }

        /// <summary>
        /// Example of a sequence variation:
        /// @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30"
        /// Reference allele: A
        /// Alternate allele: G
        /// Snpeff annotation: ANN=G||||||||||||||||
        /// Allele Index: 1
        /// Format: GT:AD:DP
        /// Genotype: 1/1
        /// Allelic depths: 30,30
        /// Homozygous reference calls: 30
        /// Heterozygous calls: 30
        /// </summary>

        [Test]
        [TestCase("NNNPPP", "accession", "A", @"not applied", 1, 3, 0, 3, 0)]
        [TestCase("NNNPPP", "accession", "A", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 1, 3, 0, 3, 0)]
        [TestCase("NNNPPP", "accession", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 2, 3, 0, 3, 1)]
        public static void TestCombos(string proteinSequence, string accession, string variantAA, string sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, int numModifiedResiduesNP)
        {
            List<SpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                                                                      new Modification(_originalId: "16",  _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", variantAA, sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantBioPolymers();

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications)).First();

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943 + 15.994915).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<IBioPolymerWithSetMods> { modPep };
            SpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            match.SetMs2Scan(scan.TheScan);
            allIdentifications = new List<SpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(numModHashes));
            Assert.That(res.Mods["accession"].Count, Is.EqualTo(numModifiedResidues));
            Assert.That(res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("21")).Count(), Is.EqualTo(numModifiedResiduesN));
            Assert.That(res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("16")).Count(), Is.EqualTo(numModifiedResiduesP));
            res.Mods.TryGetValue("accession_N1P", out var hash);
            Assert.That((hash ?? new HashSet<Tuple<int, Modification>>()).Count, Is.EqualTo(numModifiedResiduesNP));
        }

        [Test]
        public static void TestPtmBeforeVariant()
        {
            
            List<SpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                                                                      new Modification(_originalId: "16",  _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            var parentProtein = new Protein("NNNPPP", "protein", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(6, 6, "P", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) });
            var variantProteins = parentProtein.GetVariantBioPolymers();

            CommonParameters commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(commonParameters.DigestionParams, new List<Modification>(), variableModifications)).First();

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943 + 15.994915).ToMz(1), 1, "filepath", new CommonParameters());

            var peptidesWithSetModifications = new List<IBioPolymerWithSetMods> { modPep };
            SpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            match.SetMs2Scan(scan.TheScan);
            allIdentifications = new List<SpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null);
            var res = (GptmdResults)engine.Run();
        }

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
            var variantProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, null, false, null, out var unknownModifications, minAlleleDepth: 0);
            var variantProtein = variantProteins[0];
            var variantDecoy = variantProteins[1];
            Assert.That(unknownModifications.Count, Is.EqualTo(0));

            Assert.That(variantProteins.Count, Is.EqualTo(2)); // target & decoy
            Assert.That(variantProteins[0].OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
            List<int> foundResidueIndicies = variantProtein.OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            List<int> expectedResidueIndices = new List<int>() { 1, 3 };
            Assert.That(foundResidueIndicies, Is.EquivalentTo(expectedResidueIndices));
            Assert.That(variantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
            foundResidueIndicies = variantDecoy.OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            expectedResidueIndices = new List<int>() { 4, 6 }; //originally modified residues are now at the end in the decoy
            Assert.That(foundResidueIndicies, Is.EquivalentTo(expectedResidueIndices));

            var thisOk = unknownModifications;//for debugging
            var commonParamsAtThisPoint = task1.CommonParameters.DigestionParams; //for debugging

            var digestedList = variantProteins[0].GetVariantBioPolymers()[0].Digest(task1.CommonParameters.DigestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.That(digestedList.Count, Is.EqualTo(4));

            //Set Peptide with 1 mod at position 3
            var pepWithSetMods1 = digestedList[1];

            //Finally Write MZML file
            Assert.That(pepWithSetMods1.FullSequence, Is.EqualTo("PEK[type:acetyl on K]TID"));//this might be base sequence
            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { pepWithSetMods1 });
            string mzmlName = @"hello.mzML";
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

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

            Assert.That(GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex), Is.EqualTo(result));
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
                bool modFits = GptmdEngine.ModFits(mmPhospho, protein, i + 1, pep.Length, pep.OneBasedStartResidue + i);

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
                bool modFits = GptmdEngine.ModFits(mmPhospho, protein, i + 1, pep.Length, pep.OneBasedStartResidue + i);

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

        [Test]
        public static void AddsModsToDecoys()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestDecoyWriting");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);

            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb_Reversed.fasta";
            var commonParameters = new CommonParameters();
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters = new() { (mzmlName, commonParameters) };
            var gptmdParams = new GptmdParameters
            {
                ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                    b.ModificationType.Equals("Common Artifact")
                    || b.ModificationType.Equals("Common Biological")
                    || b.ModificationType.Equals("Metal")
                    || b.ModificationType.Equals("Less Common")
                ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList()
            };

            // Load Information just like the task does
            var spectraFile = MsDataFileReader.GetDataFile(mzmlName).LoadAllStaticData();
            var ms2Scans = MetaMorpheusTask.GetMs2Scans(spectraFile, mzmlName, commonParameters)
                .OrderBy(b => b.PrecursorMass)
                .ToArray();
            SpectralMatch[] psmArray = new SpectralMatch[ms2Scans.Length];

            var modConversionMethod = typeof(MetaMorpheusTask).GetMethod("LoadModifications", BindingFlags.NonPublic | BindingFlags.Instance);
            object[] parameters = new object[] { "taskId", null, null, null };
            modConversionMethod!.Invoke(new SearchTask(), parameters);
            List<Modification> variableMods = (List<Modification>)parameters[1];
            List<Modification> fixedMods = (List<Modification>)parameters[2];
            List<string> localizableMods = (List<string>)parameters[3];
            List<Modification> gptmdModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b =>
                gptmdParams.ListOfModsGptmd.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            var combos = GptmdTask.LoadCombos(gptmdModifications, true).ToList();
            MassDiffAcceptor searchMode = new DotMassDiffAcceptor("", GptmdTask.GetAcceptableMassShifts(fixedMods, variableMods, gptmdModifications, combos), commonParameters.PrecursorMassTolerance);

            var dbLoadingEngine = new DatabaseLoadingEngine(commonParameters, fileSpecificParameters, ["t1"], [new(fastaName, false)], "t1", DecoyType.Reverse, true, localizableMods);
            var proteins = (dbLoadingEngine.Run() as DatabaseLoadingEngineResults)!.BioPolymers;
            // Run the search engine and FDR just like the task. 
            new ClassicSearchEngine(psmArray, ms2Scans, variableMods, fixedMods, null, null, null,
                proteins, searchMode, commonParameters, fileSpecificParameters, null, new List<string> { "t1", "Individual Spectra Files", mzmlName }, false).Run();

            new FdrAnalysisEngine(psmArray.Where(p => p != null).ToList(), searchMode.NumNotches, commonParameters, fileSpecificParameters, new List<string> { "t1" }, doPEP: false).Run();

            // Manually set all q value notches to be 0. This simulates a big dataset in which decoys would pass the filtering. 
            foreach (var psm in psmArray.Where(p => p != null))
            {
                psm.SetMs2Scan(ms2Scans[psm.ScanIndex].TheScan);
                var fdrInfo = psm.FdrInfo;
                psm.SetFdrValues(fdrInfo.CumulativeTarget, fdrInfo.CumulativeDecoy, fdrInfo.QValue, fdrInfo.CumulativeTargetNotch, fdrInfo.CumulativeDecoyNotch, 0, fdrInfo.PEP, fdrInfo.PEP_QValue);
            }

            Dictionary<string, HashSet<Tuple<int, Modification>>> allModDictionary = new();
            new GptmdEngine(psmArray.Where(p => p != null).ToList(), gptmdModifications, combos, new() { { mzmlName, commonParameters.PrecursorMassTolerance } }, commonParameters, fileSpecificParameters, new List<string> { "t1" }, allModDictionary, gptmdParams.GptmdFilters).Run();


            // Assert that mods were added to the decoys by the engine. 
            var decoyModDicts = allModDictionary.Where(kvp => kvp.Key.StartsWith("DECOY_")).ToList();
            Assert.That(decoyModDicts.Count, Is.GreaterThan(0));
            Assert.That(decoyModDicts.SelectMany(kvp => kvp.Value).Count, Is.GreaterThan(0));

            // Filter proteins and test filtering
            var toWrite = GptmdTask.GetBioPolymersToWrite(proteins, allModDictionary).ToList();
            Assert.That(toWrite.Count, Is.GreaterThan(0));
            Assert.That(toWrite.Count(b => b.Accession.StartsWith("DECOY_")), Is.GreaterThan(0));

            int decoyByAcCount = toWrite.Count(b => b.Accession.StartsWith("DECOY_"));
            int decoyByPropCount = toWrite.Count(b => b.IsDecoy);
            int decoysInModDictCount = allModDictionary.Count(kvp => kvp.Key.StartsWith("DECOY_"));
            Assert.That(decoyByAcCount, Is.EqualTo(decoyByPropCount));
            int targetByAcCount = toWrite.Count(b => !b.Accession.StartsWith("DECOY_"));
            int targetByPropCount = toWrite.Count(b => !b.IsDecoy);
            Assert.That(targetByAcCount, Is.EqualTo(targetByPropCount));

            Assert.That(proteins.Count(p => !p.IsDecoy), Is.EqualTo(targetByPropCount), "Not all targets are being retained :(");
            Assert.That(decoyByAcCount, Is.EqualTo(decoysInModDictCount), "Not all decoys with mods added are being retained :(");

            // write the xml 
            string outputFile = Path.Combine(outputFolder, "moddedDB.xml");
            ProteinDbWriter.WriteXmlDatabase(allModDictionary, toWrite, outputFile);

            // read the xml and do not generate decoys
            dbLoadingEngine = new DatabaseLoadingEngine(commonParameters, fileSpecificParameters, ["t1"], [new(outputFile, false)], "t1", DecoyType.None, true, localizableMods);
            var reloadedProteinsNoDecoyGen = (dbLoadingEngine.Run() as DatabaseLoadingEngineResults)!.BioPolymers;
            Assert.That(reloadedProteinsNoDecoyGen.Count, Is.EqualTo(targetByAcCount + decoysInModDictCount));
            var targetCount = reloadedProteinsNoDecoyGen.Count(p => !p.IsDecoy);
            Assert.That(targetCount, Is.EqualTo(targetByAcCount));

            // Only decoys with mods written should have been loaded back in 
            var decoyCount = reloadedProteinsNoDecoyGen.Count(p => p.IsDecoy);
            Assert.That(decoyCount, Is.EqualTo(decoysInModDictCount));
            foreach (var reloadedProt in reloadedProteinsNoDecoyGen.Where(p => p.IsDecoy))
            {
                Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.GreaterThan(0));
                Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.GreaterThan(0));
                Assert.That(allModDictionary.ContainsKey(reloadedProt.Accession));
            }

            // read the xml and generate decoys. 
            dbLoadingEngine = new DatabaseLoadingEngine(commonParameters, fileSpecificParameters, ["t1"], [new(outputFile, false)], "t1", DecoyType.Reverse, true, localizableMods);
            var reloadedProteinsWithDecoyGen = (dbLoadingEngine.Run() as DatabaseLoadingEngineResults)!.BioPolymers;

            // should be the same number of targets and decoys as before. 
            Assert.That(reloadedProteinsWithDecoyGen.Count(b => !b.IsDecoy), Is.EqualTo(reloadedProteinsWithDecoyGen.Count(b => b.IsDecoy)));
            Assert.That(reloadedProteinsWithDecoyGen.Count(b => b.Accession.StartsWith("DECOY_")), Is.EqualTo(reloadedProteinsWithDecoyGen.Count(b => !b.Accession.StartsWith("DECOY_"))));

            foreach (var reloadedProt in reloadedProteinsWithDecoyGen)
            {
                string targetAccession;
                string decoyAccession;
                if (reloadedProt.IsDecoy)
                {
                    decoyAccession = reloadedProt.Accession;
                    targetAccession = reloadedProt.Accession.Replace("DECOY_", "");

                }
                else
                {
                    targetAccession = reloadedProt.Accession;
                    decoyAccession = "DECOY_" + reloadedProt.Accession;
                }

                allModDictionary.TryGetValue(targetAccession, out var targetModSet);
                allModDictionary.TryGetValue(decoyAccession, out var decoyModSet);

                // Reloaded Target
                if (!reloadedProt.IsDecoy)
                {
                    // No mod added via GPTMD - No Modifications should be present
                    if (targetModSet == null)
                    {
                        Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(0));
                        Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.EqualTo(0));
                    }
                    // Mods added via GPTMD - Target Mods will match 1 to 1 with what gptmd added
                    else
                    {
                        int expectedMods = targetModSet.Count;
                        Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(expectedMods));
                        Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.EqualTo(expectedMods));
                    }
                }
                // Reloaded Decoy
                else
                {
                    // Target for decoy had no mods added via gptmd 
                    if (targetModSet == null)
                    {
                        // Decoy had no mods added in gptmd
                        if (decoyModSet == null)
                        {
                            int expectedMods = 0;
                            Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(expectedMods));
                            Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.EqualTo(expectedMods));
                        }
                        // Decoy had mods added in gptmd
                        else
                        {
                            int expectedMods = decoyModSet.Count;
                            Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(expectedMods));
                            Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.EqualTo(expectedMods));
                        }
                    }
                    // Target for decoy had mods added via gptmd 
                    else
                    {
                        // Decoy had no mods added in gptmd - mod count is equal to that of target
                        if (decoyModSet == null)
                        {
                            int expectedMods = targetModSet.Count;
                            Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(expectedMods));
                            Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.EqualTo(expectedMods));
                        }
                        // Decoy had mods added in gptmd - mod count is equal to target + decoy mods added from gptmd
                        else
                        {
                            int expectedMods = targetModSet.Count + decoyModSet.Count;
                            Assert.That(reloadedProt.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(expectedMods));
                            Assert.That(reloadedProt.OriginalNonVariantModifications.Count, Is.EqualTo(expectedMods));
                        }
                    }
                }
            }

            Directory.Delete(outputFolder, true);
        }
    }
}