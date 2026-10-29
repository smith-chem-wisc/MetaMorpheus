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
using Omics.Digestion;
using Readers;

namespace Test
{
    [TestFixture]
    public static class GptmdEngineTest
    {
        [Test]
        [TestCase("NNNNN", "accession", @"not applied", 5, "accession")]
        // the peptide spans the applied N1A, so its sites were seen only with the variant and go to the variant's own modifications
        [TestCase("NNNNN", "accession", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 4, "accession_N1A")]
        public static void TestGptmdEngine(string proteinSequence, string accession, string sequenceVariantDescription, int numModifiedResidues, string accessionOfSites)
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

            var recordingFilter = new RecordingGptmdFilter();
            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), null, new List<string>(), null, new List<IGptmdFilter> { recordingFilter });
            res = (GptmdResults)engine.Run();
            Assert.That(res.Mods.Count, Is.EqualTo(1));
            Assert.That(res.Mods[accessionOfSites].Count, Is.EqualTo(numModifiedResidues));
            Assert.That(recordingFilter.Sites, Is.Not.Empty);
            Assert.That(recordingFilter.Sites.All(site => site >= 1 && site <= modPep.Length), Is.True);
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
        [TestCase("NNNPPP", "accession", "A", @"not applied", 1, 3, 0, 3, "accession")]
        // the peptide spans the applied variant, so its sites (P1 on the variant, P4-P6 beside it) go to the variant's own modifications
        [TestCase("NNNPPP", "accession", "A", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 1, 3, 0, 3, "accession_N1A")]
        [TestCase("NNNPPP", "accession", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", 1, 4, 0, 4, "accession_N1P")]
        public static void TestCombos(string proteinSequence, string accession, string variantAA, string sequenceVariantDescription, int numModHashes, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP, string accessionOfSites)
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
            Assert.That(res.Mods[accessionOfSites].Count, Is.EqualTo(numModifiedResidues));
            Assert.That(res.Mods[accessionOfSites].Where(b => b.Item2.OriginalId.Equals("21")).Count(), Is.EqualTo(numModifiedResiduesN));
            Assert.That(res.Mods[accessionOfSites].Where(b => b.Item2.OriginalId.Equals("16")).Count(), Is.EqualTo(numModifiedResiduesP));
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
            List<int> foundResidueIndices = variantProtein.OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            List<int> expectedResidueIndices = new List<int>() { 1, 3 };
            Assert.That(foundResidueIndices, Is.EquivalentTo(expectedResidueIndices));
            Assert.That(variantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
            foundResidueIndices = variantDecoy.OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            expectedResidueIndices = new List<int>() { 4, 6 }; //originally modified residues are now at the end in the decoy
            Assert.That(foundResidueIndices, Is.EquivalentTo(expectedResidueIndices));

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

            var toWrite = GptmdTask.GetBioPolymersToWrite(proteins, allModDictionary, false).ToList();
            Assert.That(toWrite.Count, Is.GreaterThan(0));
            Assert.That(toWrite.Count(b => b.Accession.StartsWith("DECOY_")), Is.EqualTo(0));

            // Filter proteins and test filtering
            toWrite = GptmdTask.GetBioPolymersToWrite(proteins, allModDictionary).ToList();
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

        // The PEX14 case: QFPPSPSAPK carries UniProt's Phosphoserine at S5, and the spectra say S7.
        [Test]
        public static void TestGptmdAddsANeighbourThatExplainsTheSpectrumBetterThanTheDatabaseSite()
        {
            var res = RunGptmdOnPhosphoPsm(PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 7));

            Assert.That(res.Mods["accession"].Select(m => (m.Item1, m.Item2.IdWithMotif)),
                Is.EquivalentTo(new[] { (7, "Phosphoserine on S") }));
        }

        [Test]
        public static void TestGptmdAddsNothingWhenTheDatabaseSiteExplainsTheSpectrumBest()
        {
            var res = RunGptmdOnPhosphoPsm(PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 5));

            Assert.That(res.Mods.ContainsKey("accession"), Is.False);
        }

        // Without b5, b6, y4 and y5 nothing in the spectrum separates S5 from S7, so both go to the search.
        [Test]
        public static void TestGptmdAddsANeighbourTheSpectrumCannotTellFromTheDatabaseSite()
        {
            var res = RunGptmdOnPhosphoPsm(PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 7, keepFragment: NotSiteDetermining));

            Assert.That(res.Mods["accession"].Select(m => (m.Item1, m.Item2.IdWithMotif)),
                Is.EquivalentTo(new[] { (7, "Phosphoserine on S") }));
        }

        // UniProt's Phosphoserine carries no neutral loss; GPTMD's Phosphorylation on S carries H3PO4. The S7 that ties
        // UniProt's S5 is written as Phosphoserine, so in the search the two placements score alike on a spectrum whose
        // neutral-loss ions all hold both serines, and the site can be reported as ambiguous. Written under the GPTMD
        // name, S7 would collect those neutral losses alone and be called on its own.
        [Test]
        public static void TestGptmdWritesAMoveToTheSameAminoAcidAsTheDatabaseModSoTheSearchScoresBothPlacementsAlike()
        {
            var psms = PsmToDatabaseSites("QFPPSPSAPK", new() { { 5, UniProtMod("Phosphoserine on S") } }, new() { { 7, GptmdMod("Phosphorylation on S") } },
                keepFragment: NotSiteDetermining, neutralLosses: true);

            var res = RunGptmdOnPhosphoPsm(psms);
            Assert.That(res.Mods["accession"].Select(m => (m.Item1, m.Item2.IdWithMotif)), Is.EquivalentTo(new[] { (7, "Phosphoserine on S") }));

            var protein = (Protein)psms[0].BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.Parent;
            var reloaded = WriteAndReload(res.Mods, protein, nameof(TestGptmdWritesAMoveToTheSameAminoAcidAsTheDatabaseModSoTheSearchScoresBothPlacementsAlike)).Single();
            var placements = reloaded.Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .Where(p => p.AllModsOneIsNterminus.Count == 1).ToList();
            Assert.That(placements.Select(p => p.AllModsOneIsNterminus.Single().Key), Is.EquivalentTo(new[] { 6, 8 }));

            var scores = placements.Select(p => SearchScore(p, psms[0].Ms2Scan)).ToList();
            Assert.That(scores[0], Is.EqualTo(scores[1]).Within(1e-9));

            var underGptmdName = new Protein("QFPPSPSAPK", "gptmdName", oneBasedModifications: new Dictionary<int, List<Modification>> { { 7, new List<Modification> { GptmdMod("Phosphorylation on S") } } })
                .Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);
            Assert.That(SearchScore(underGptmdName, psms[0].Ms2Scan), Is.GreaterThan(scores[0] + 1), "control: the spectrum's neutral losses favour the GPTMD name");
        }

        [Test]
        public static void TestGptmdMovesADatabasePhosphoserineToAThreonine()
        {
            var res = RunGptmdOnPhosphoPsm(PsmToDatabasePhosphosite("QFPPSPTAPK", databaseSite: 5, trueSite: 7));

            Assert.That(res.Mods["accession"].Select(m => (m.Item1, m.Item2.IdWithMotif)),
                Is.EquivalentTo(new[] { (7, "Phosphorylation on T") }));
        }

        [Test]
        public static void TestGptmdLeavesDatabaseModsWhereTheyAreWhenAlternativeSitesAreOff()
        {
            var res = RunGptmdOnPhosphoPsm(PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 7), addAlternativeSites: false);

            Assert.That(res.Mods.ContainsKey("accession"), Is.False);
        }

        // An alternative site passes the same filters as any GPTMD candidate. ImprovedScoreFilter is strict, so it
        // drops a tie with the database site; the flanking filter needs ions on both sides of S7, which go with the
        // site-determining ions.
        [Test]
        [TestCase(typeof(ImprovedScoreFilter), true, true)]
        [TestCase(typeof(ImprovedScoreFilter), false, false)]
        [TestCase(typeof(FlankingIonCoverageFilter), true, true)]
        [TestCase(typeof(FlankingIonCoverageFilter), false, false)]
        [TestCase(typeof(DualDirectionalIonCoverageFilter), false, true)]
        [TestCase(typeof(UniDirectionalIonCoverageFilter), false, true)]
        public static void TestGptmdFiltersApplyToAlternativeSitesForDatabaseMods(Type filterType, bool siteDeterminingIonsInSpectrum, bool expectAdded)
        {
            var psms = PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 7,
                keepFragment: siteDeterminingIonsInSpectrum ? null : NotSiteDetermining);
            var res = RunGptmdOnPhosphoPsm(psms, new List<IGptmdFilter> { (IGptmdFilter)Activator.CreateInstance(filterType) });

            Assert.That(res.Mods.ContainsKey("accession") && res.Mods["accession"].Any(m => m.Item1 == 7), Is.EqualTo(expectAdded));
        }

        [Test]
        public static void TestGptmdAlternativeSitesForDatabaseModsAreOffByDefault()
        {
            Assert.That(new GptmdParameters().AddAlternativeSitesForDatabaseMods, Is.False);

            var psms = PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 7);
            var res = (GptmdResults)new GptmdEngine(psms, PhosphoMods(), new List<Tuple<double, double>>(),
                new Dictionary<string, Tolerance> { { "filepath", new PpmTolerance(10) } }, new CommonParameters(), null, new List<string>(), null).Run();

            Assert.That(res.Mods.ContainsKey("accession"), Is.False);
        }

        // A tie between two residues of the same amino acid is a site the spectrum cannot place, and both go to the
        // search. A tie between different amino acids is no evidence for the other one, so it is not added.
        [Test]
        public static void TestGptmdDoesNotMoveADatabaseModToAnotherAminoAcidOnATie()
        {
            var res = RunGptmdOnPhosphoPsm(PsmToDatabasePhosphosite("QFPPSPTAPK", databaseSite: 5, trueSite: 7, keepFragment: NotSiteDetermining));

            Assert.That(res.Mods.ContainsKey("accession"), Is.False);
        }

        // UniProt's 4-hydroxyproline and GPTMD's Hydroxylation on N share a formula (O). On a tie the hydroxyl stays on
        // proline: N7 is not added.
        [Test]
        public static void TestGptmdDoesNotMoveADatabaseHydroxyprolineToAnAsparagineOnATie()
        {
            Modification hydroxyproline = UniProtMod("4-hydroxyproline on P");
            var psms = PsmToDatabaseSites("QFPAPGNAPK", new() { { 5, hydroxyproline } }, new() { { 7, GptmdMod("Hydroxylation on N") } },
                keepFragment: p => NotDeterminingBetween(p, 5, 7, 10));

            var res = RunGptmd(psms, DefaultGptmdMods());

            Assert.That(res.Mods.TryGetValue("accession", out var added) ? added.Where(m => m.Item1 == 7) : [], Is.Empty);
        }

        // Each database modification is moved on its own; the other stays where the database put it, and a residue that
        // already carries one cannot take a second.
        [Test]
        public static void TestGptmdMovesEachOfTwoDatabaseModsOnOnePeptideOnItsOwn()
        {
            Modification phosphoserine = UniProtMod("Phosphoserine on S");
            Modification phosphoS = GptmdMod("Phosphorylation on S");
            var psms = PsmToDatabaseSites("QSFPSPSAPK", new() { { 2, phosphoserine }, { 5, phosphoserine } }, new() { { 2, phosphoS }, { 7, phosphoS } });

            var res = RunGptmd(psms, PhosphoMods());

            Assert.That(res.Mods["accession"].Select(m => (m.Item1, m.Item2.IdWithMotif)), Is.EquivalentTo(new[] { (7, "Phosphoserine on S") }));
        }

        // When a fixed modification targets the residue a database modification sits on, the search puts the fixed
        // modification back once the database one moves away, so the moved form is heavier than this precursor and
        // cannot match the spectrum it came from. Here the spectrum even prefers the moved form.
        [Test]
        [TestCase("QFPACGSAPK", "Phosphocysteine on C", 5, "Phosphorylation on S", 7, new[] { "Carbamidomethyl on C" })]
        [TestCase("QFPAGDLKR", "N6,N6-dimethyllysine on K", 8, "Dimethylation on R", 9, new[] { "TMT10 on K", "TMT10 on X" })]
        public static void TestGptmdDoesNotMoveADatabaseModOffAResidueThatTakesAFixedMod(string sequence, string databaseModId, int databaseSite,
            string movedModId, int movedSite, string[] fixedModIds)
        {
            var fixedMods = fixedModIds.Select(id => GlobalVariables.AllModsKnown.First(m => m.IdWithMotif == id
                && m.ModificationType is "Common Fixed" or "Multiplex Label")).ToList();
            var psms = PsmToDatabaseSites(sequence, new() { { databaseSite, UniProtMod(databaseModId) } }, new() { { movedSite, GptmdMod(movedModId) } },
                fixedMods: fixedMods);

            var withoutFixedMods = RunGptmd(PsmToDatabaseSites(sequence, new() { { databaseSite, UniProtMod(databaseModId) } },
                new() { { movedSite, GptmdMod(movedModId) } }), DefaultGptmdMods());
            Assert.That(withoutFixedMods.Mods["accession"].Any(m => m.Item1 == movedSite), Is.True, "control: the move is made when nothing fixed sits on the residue");

            var res = RunGptmd(psms, DefaultGptmdMods(), fixedMods: fixedMods);

            Assert.That(res.Mods.ContainsKey("accession"), Is.False);
        }

        // The ordinary pass's cleavage-blocking rule holds for moves too. UniProt's acetyl-K at K4 of AAAKAAAAK, and a
        // spectrum that cannot tell K4 from K9: with the search respecting cleavage-blocking modifications, the tie is not
        // moved onto K9, the residue trypsin cut after; with it off, K9 is added.
        [Test]
        [TestCase(true, new int[0])]
        [TestCase(false, new[] { 9 })]
        public static void TestGptmdDoesNotMoveACleavageBlockingDatabaseModOntoTheResidueTheProteaseCutAfter(bool respectCleavageBlockingModifications,
            int[] expectedSites)
        {
            var protein = new Protein("AAAKAAAAKAAAAR", "accession",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, new List<Modification> { UniProtMod("N6-acetyllysine on K") } } });
            var digestionParams = new DigestionParams(minPeptideLength: 5, respectCleavageBlockingModifications: respectCleavageBlockingModifications);
            var matched = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.BaseSequence == "AAAKAAAAK" && p.AllModsOneIsNterminus.ContainsKey(5));
            var source = new Protein("AAAKAAAAK", "source",
                    oneBasedModifications: new Dictionary<int, List<Modification>> { { 9, new List<Modification> { GptmdMod("Acetylation on K") } } })
                .Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);

            var res = RunGptmd(new List<SpectralMatch> { PsmWithSpectrum(matched, SpectrumOf(source, p => NotDeterminingBetween(p, 4, 9, 9))) },
                new List<Modification> { GptmdMod("Acetylation on K") });

            Assert.That(res.Mods.TryGetValue("accession", out var added) ? added.Select(m => m.Item1) : [], Is.EquivalentTo(expectedSites));
        }

        // A database site under a GPTMD name may come from an earlier round or from an annotation written with
        // MetaMorpheus's names, and the database cannot say which, so it is moved like any other when the spectrum
        // prefers another residue. Rounds still do not compound: see the next test.
        [Test]
        public static void TestGptmdMovesADatabaseSiteWrittenUnderAGptmdName()
        {
            var annotated = new Protein("QFPPSPSAPK", "accession",
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 5, new List<Modification> { GptmdMod("Phosphorylation on S") } } });
            var readBack = WriteAndReload(new(), annotated, nameof(TestGptmdMovesADatabaseSiteWrittenUnderAGptmdName)).Single();
            var matched = readBack.Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);
            var source = new Protein("QFPPSPSAPK", "source", oneBasedModifications: new Dictionary<int, List<Modification>> { { 7, new List<Modification> { GptmdMod("Phosphorylation on S") } } })
                .Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);

            var res = RunGptmdOnPhosphoPsm(new List<SpectralMatch> { PsmWithSpectrum(matched, SpectrumOf(source, null)) });

            Assert.That(res.Mods["accession"].Select(m => (m.Item1, m.Item2.IdWithMotif)), Is.EquivalentTo(new[] { (7, "Phosphorylation on S") }));
        }

        // The database offers S2 (an earlier round's site) and S5 (UniProt's); the PSM is to S5, and the spectrum says S2
        // but cannot tell S5 from S7. S7 ties the placement the PSM carries, not the best one the database offers, so it
        // is not added. Measuring moves against the PSM's own placement alone would add it, and each round would then
        // add sites a little further from what the spectra support.
        [Test]
        public static void TestGptmdComparesMovesWithEverySiteTheDatabaseAlreadyOffers()
        {
            var protein = new Protein("QSPPSPSAPK", "accession", oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 2, new List<Modification> { GptmdMod("Phosphorylation on S") } },
                { 5, new List<Modification> { UniProtMod("Phosphoserine on S") } },
            });
            var matched = protein.Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1 && p.AllModsOneIsNterminus.ContainsKey(6));
            var source = new Protein("QSPPSPSAPK", "source", oneBasedModifications: new Dictionary<int, List<Modification>> { { 2, new List<Modification> { GptmdMod("Phosphorylation on S") } } })
                .Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);

            var res = RunGptmdOnPhosphoPsm(new List<SpectralMatch> { PsmWithSpectrum(matched, SpectrumOf(source, NotSiteDetermining)) });

            Assert.That(res.Mods.ContainsKey("accession"), Is.False);
        }

        // Round one adds S7 beside UniProt's S5 (the spectrum cannot tell them apart). Written out and read back, the same
        // spectrum adds nothing in round two, whether it is matched to S5 or to S7.
        [Test]
        public static void TestGptmdAddsNothingInASecondRoundOverItsOwnDatabase()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestGptmdAddsNothingInASecondRoundOverItsOwnDatabase));
            Directory.CreateDirectory(outputFolder);
            string database = Path.Combine(outputFolder, "roundOneGPTMD.xml");

            var roundOnePsms = PsmToDatabasePhosphosite("QFPPSPSAPK", databaseSite: 5, trueSite: 7, keepFragment: NotSiteDetermining);
            var roundOne = RunGptmdOnPhosphoPsm(roundOnePsms);
            Assert.That(roundOne.Mods["accession"].Select(m => m.Item1), Is.EquivalentTo(new[] { 7 }));

            ProteinDbWriter.WriteXmlDatabase(roundOne.Mods, new List<Protein> { (Protein)roundOnePsms[0].BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer.Parent }, database);
            var reloaded = ProteinDbLoader.LoadProteinXML(database, true, DecoyType.None, GlobalVariables.AllModsKnown, false, null, out var unknownMods).Single();
            Assert.That(unknownMods, Is.Empty);
            Assert.That(reloaded.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(new[] { 5, 7 }));

            var spectrumSource = roundOnePsms[0];
            foreach (int matchedSite in new[] { 5, 7 })
            {
                var matched = reloaded.Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                    .First(p => p.AllModsOneIsNterminus.Count == 1 && p.AllModsOneIsNterminus.ContainsKey(matchedSite + 1));
                var roundTwo = RunGptmdOnPhosphoPsm(new List<SpectralMatch> { PsmWithSpectrum(matched, spectrumSource.Ms2Scan) });

                Assert.That(roundTwo.Mods.ContainsKey("accession"), Is.False, $"round two, matched to S{matchedSite}");
            }

            Directory.Delete(outputFolder, true);
        }

        // A new best keeps every earlier candidate still within tolerance of it, so the order candidates are tried in
        // does not change which are kept.
        [Test]
        public static void TestGptmdKeepsEveryCandidateWithinToleranceOfTheBestWhateverTheOrder()
        {
            var scored = new[] { ("a", 1.0), ("b", 1.08), ("c", 1.12) };
            int[][] orders = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

            foreach (var order in orders)
            {
                var inOrder = order.Select(i => scored[i]).ToList();
                Assert.That(GptmdEngine.WithinToleranceOfBest(inOrder, 0, 0.1), Is.EquivalentTo(new[] { "b", "c" }), string.Join(",", order));
            }

            // the baseline (the database's own placement) counts toward the best but is not returned
            Assert.That(GptmdEngine.WithinToleranceOfBest(scored, 1.15, 0.1), Is.EquivalentTo(new[] { "b", "c" }));
            Assert.That(GptmdEngine.WithinToleranceOfBest(scored, 1.3, 0.1), Is.Empty);
        }

        // A sequence variant carries its own modifications (SequenceVariation.OneBasedModifications), keyed in the
        // coordinates of the protein with that variant applied, and ProteinDbWriter writes the GPTMD dictionary entry
        // under the single-variant accession as that variant's subfeatures. A site on a peptide that spans a variant's
        // edited residues was seen only with the variant, so it goes there, even on a residue the variant does not edit.
        [Test]
        [TestCase("NNNPPP", 4, "P", "A", "accession_P4A", new[] { 1, 2, 3 }, "NNNAPP")]
        // after an insertion: N7 of the variant protein, which is N5 of the consensus
        [TestCase("PEPTNPEP", 2, "E", "EGG", "accession_E2EGG", new[] { 7 }, "PEGGPTNPEP")]
        // a deletion edits no residue it keeps: NNNINNNK's I-N junction exists only with ID4I, so all six sites go to it
        [TestCase("NNNIDNNNK", 4, "ID", "I", "accession_ID4I", new[] { 1, 2, 3, 5, 6, 7 }, "NNNINNNK")]
        public static void TestGptmdWritesASiteSeenOnlyOnAVariantPeptideToThatVariant(string sequence, int variantSite, string original, string variant,
            string expectedAccession, int[] expectedSites, string variantSequence)
        {
            var protein = new Protein(sequence, "accession", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(variantSite, original, variant, HomozygousVcf) });

            var res = RunGptmdOnUnexplainedMass(protein, GptmdMod("Hydroxylation on N"));

            Assert.That(res.Mods.Keys, Is.EquivalentTo(new[] { expectedAccession }));
            Assert.That(res.Mods[expectedAccession].Select(m => m.Item1), Is.EquivalentTo(expectedSites));

            var reloaded = WriteAndReloadVariant(res.Mods, protein, nameof(TestGptmdWritesASiteSeenOnlyOnAVariantPeptideToThatVariant));
            Assert.That(reloaded.BaseSequence, Is.EqualTo(variantSequence));
            Assert.That(reloaded.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(expectedSites));
            Assert.That(expectedSites.All(site => variantSequence[site - 1] == 'N'));
            Assert.That(reloaded.ConsensusVariant.OneBasedPossibleLocalizedModifications, Is.Empty);
        }

        // Two applied variants, the first an insertion. A site on the second variant's edited residue goes to that
        // variant alone, under its own single-variant accession (the writer drops any other accession) and at its
        // position with the insertion undone; a site the peptide spans both variants for goes to the consensus. Read
        // back, each lands on the residue it was found on.
        [Test]
        [TestCase("Hydroxylation on N", "accession_I12N", new[] { 12 }, new[] { 14 })]
        [TestCase("Phosphorylation on T", "accession", new[] { 4, 11 }, new[] { 6, 13 })]
        public static void TestGptmdWritesSitesOnAProteinWithTwoVariantsWhereTheyWereFound(string modId, string expectedAccession,
            int[] expectedSites, int[] sitesInVariantProtein)
        {
            var protein = new Protein("PEPTIDEPEPTIDE", "accession", sequenceVariations: new List<SequenceVariation>
            {
                new SequenceVariation(2, "E", "EGG", HomozygousVcf),
                new SequenceVariation(12, "I", "N", HomozygousVcf.Replace("50000000", "50000100")),
            });

            var res = RunGptmdOnUnexplainedMass(protein, GptmdMod(modId));

            Assert.That(res.Mods.Keys, Is.EquivalentTo(new[] { expectedAccession }));
            Assert.That(res.Mods[expectedAccession].Select(m => m.Item1), Is.EquivalentTo(expectedSites));

            var reloaded = WriteAndReloadVariant(res.Mods, protein, nameof(TestGptmdWritesSitesOnAProteinWithTwoVariantsWhereTheyWereFound));
            Assert.That(reloaded.BaseSequence, Is.EqualTo("PEGGPTIDEPEPTNDE"));
            Assert.That(reloaded.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(sitesInVariantProtein));
            Assert.That(reloaded.OneBasedPossibleLocalizedModifications.All(kv => kv.Value.All(m => m.Target.ToString() == reloaded.BaseSequence[kv.Key - 1].ToString())));
        }

        // A peptide that does not reach a variant's edited residues is the consensus protein's peptide too, so its sites
        // go to the consensus protein, at their consensus positions (the insertion upstream undone).
        [Test]
        public static void TestGptmdWritesASiteOnAPeptideAwayFromTheVariantToTheConsensus()
        {
            var protein = new Protein("PEPTIDEKNNNNNP", "accession", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(2, "E", "EGG", HomozygousVcf) });

            var res = RunGptmdOnUnexplainedMass(protein, GptmdMod("Hydroxylation on N"), peptide: "NNNNNP");

            Assert.That(res.Mods.Keys, Is.EquivalentTo(new[] { "accession" }));
            Assert.That(res.Mods["accession"].Select(m => m.Item1), Is.EquivalentTo(new[] { 9, 10, 11, 12, 13 }));
        }

        // A database modification moved onto a residue a variant creates goes to that variant; the database's own site,
        // upstream of the variant, stays on the consensus protein.
        [Test]
        public static void TestGptmdMovesADatabaseModOntoAResidueAVariantCreates()
        {
            var protein = new Protein("QFPPSPAAPK", "accession", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(7, "A", "S", HomozygousVcf) },
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 5, new List<Modification> { UniProtMod("Phosphoserine on S") } } });
            var matched = protein.GetVariantBioPolymers().Single()
                .Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);
            Assert.That(matched.BaseSequence, Is.EqualTo("QFPPSPSAPK"));
            var source = new Protein("QFPPSPSAPK", "source", oneBasedModifications: new Dictionary<int, List<Modification>> { { 7, new List<Modification> { GptmdMod("Phosphorylation on S") } } })
                .Digest(new CommonParameters().DigestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.AllModsOneIsNterminus.Count == 1);

            var res = RunGptmdOnPhosphoPsm(new List<SpectralMatch> { PsmWithSpectrum(matched, SpectrumOf(source, null)) });

            Assert.That(res.Mods.Keys, Is.EquivalentTo(new[] { "accession_A7S" }));
            Assert.That(res.Mods["accession_A7S"].Select(m => (m.Item1, m.Item2.IdWithMotif)), Is.EquivalentTo(new[] { (7, "Phosphoserine on S") }));

            var reloaded = WriteAndReloadVariant(res.Mods, protein, nameof(TestGptmdMovesADatabaseModOntoAResidueAVariantCreates));
            Assert.That(reloaded.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(new[] { 5, 7 }));
            Assert.That(reloaded.ConsensusVariant.OneBasedPossibleLocalizedModifications.Keys, Is.EquivalentTo(new[] { 5 }));
        }

        // mzLib tells applied variants apart by their VCF line, so two variants on one protein need two lines.
        private const string HomozygousVcf =@"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30";

        // GPTMD on one PSM to a peptide of protein's (variant-applied) form whose precursor is heavier by mod's mass. The
        // spectrum is a single peak, so every residue mod fits ties and all of them are added.
        private static GptmdResults RunGptmdOnUnexplainedMass(Protein protein, Modification mod, string peptide = null)
        {
            var commonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var matched = protein.GetVariantBioPolymers().Single()
                .Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>())
                .Where(p => peptide == null || p.BaseSequence == peptide)
                .MaxBy(p => p.Length);
            var scan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null,
                MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, (matched.MonoisotopicMass + mod.MonoisotopicMass.Value).ToMz(1), 1, "filepath", commonParameters);
            var psm = new PeptideSpectralMatch(matched, 0, 0, 0, scanWithMass, commonParameters, new List<MatchedFragmentIon>());
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            psm.SetMs2Scan(scan);
            return RunGptmd(new List<SpectralMatch> { psm }, new List<Modification> { mod });
        }

        // Writes protein with the GPTMD dictionary and reads it back as the search would (variants applied).
        private static List<Protein> WriteAndReload(Dictionary<string, HashSet<Tuple<int, Modification>>> mods, Protein protein, string name)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, name + ".xml");
            ProteinDbWriter.WriteXmlDatabase(mods, new List<Protein> { protein }, path);
            var reloaded = ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None, GlobalVariables.AllModsKnown, false, null, out var unknownMods, minAlleleDepth: 0);
            File.Delete(path);
            Assert.That(unknownMods, Is.Empty);
            return reloaded;
        }

        private static Protein WriteAndReloadVariant(Dictionary<string, HashSet<Tuple<int, Modification>>> mods, Protein protein, string name) =>
            WriteAndReload(mods, protein, name).Single(p => p.AppliedSequenceVariations.Count > 0);

        // b5, b6, y4 and y5 of a ten-residue peptide are the ions that tell a modification on residue 5 from one on 7.
        private static bool NotSiteDetermining(Product product) => NotDeterminingBetween(product, 5, 7, 10);

        // Whether a b or y ion is one that tells a modification on fromSite from one on toSite (one-based, fromSite < toSite).
        private static bool NotDeterminingBetween(Product product, int fromSite, int toSite, int peptideLength) =>
            !(product.ProductType == ProductType.b && product.FragmentNumber >= fromSite && product.FragmentNumber < toSite
              || product.ProductType == ProductType.y && product.FragmentNumber > peptideLength - toSite && product.FragmentNumber <= peptideLength - fromSite);

        // A one-peptide protein with UniProt's Phosphoserine at databaseSite, and a PSM to that phosphoform whose spectrum
        // holds the unit-intensity b and y ions of the same peptide phosphorylated at trueSite instead.
        private static List<SpectralMatch> PsmToDatabasePhosphosite(string sequence, int databaseSite, int trueSite, Func<Product, bool> keepFragment = null)
        {
            return PsmToDatabaseSites(sequence, new() { { databaseSite, UniProtMod("Phosphoserine on S") } },
                new() { { trueSite, GptmdMod($"Phosphorylation on {sequence[trueSite - 1]}") } }, keepFragment);
        }

        // A one-peptide protein carrying databaseMods (one-based residue -> modification), and a PSM to the peptidoform
        // carrying all of them (digested with the search's fixedMods), whose spectrum holds the unit-intensity b and y ions
        // of the same peptide carrying trueMods instead, with nothing fixed (and their neutral-loss ions, if asked).
        private static List<SpectralMatch> PsmToDatabaseSites(string sequence, Dictionary<int, Modification> databaseMods,
            Dictionary<int, Modification> trueMods, Func<Product, bool> keepFragment = null, List<Modification> fixedMods = null,
            bool neutralLosses = false)
        {
            var digestionParams = new CommonParameters().DigestionParams;
            var protein = new Protein(sequence, "accession",
                oneBasedModifications: databaseMods.ToDictionary(m => m.Key, m => new List<Modification> { m.Value }));
            var matched = protein.Digest(digestionParams, fixedMods ?? new List<Modification>(), new List<Modification>())
                .First(p => p.Length == sequence.Length
                    && databaseMods.All(m => p.AllModsOneIsNterminus.TryGetValue(m.Key + 1, out var mod) && mod == m.Value));

            var source = new Protein(sequence, "source",
                    oneBasedModifications: trueMods.ToDictionary(m => m.Key, m => new List<Modification> { m.Value }))
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.Length == sequence.Length && p.AllModsOneIsNterminus.Count == trueMods.Count);

            return new List<SpectralMatch> { PsmWithSpectrum(matched, SpectrumOf(source, keepFragment, neutralLosses)) };
        }

        // The unit-intensity singly charged b and y ions of a peptidoform, keeping those keepFragment allows.
        private static MsDataScan SpectrumOf(IBioPolymerWithSetMods source, Func<Product, bool> keepFragment, bool neutralLosses = false)
        {
            var products = new List<Product>();
            source.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            double[] mzs = products
                .Where(p => (neutralLosses || p.NeutralLoss == 0) && (p.ProductType == ProductType.b || p.ProductType == ProductType.y) && (keepFragment?.Invoke(p) ?? true))
                .Select(p => p.NeutralMass.ToMz(1)).OrderBy(mz => mz).ToArray();
            return new MsDataScan(new MzSpectrum(mzs, mzs.Select(_ => 1.0).ToArray(), false), 1, 2, true, Polarity.Positive, 1,
                new MzRange(50, 2000), null, MZAnalyzerType.Orbitrap, mzs.Length, null, null, "scan=1", double.NaN, 1, 1, double.NaN, 2,
                DissociationType.HCD, 0, double.NaN);
        }

        // A q-value-0 PSM of matched to spectrum, with the precursor at matched's own mass.
        private static SpectralMatch PsmWithSpectrum(IBioPolymerWithSetMods matched, MsDataScan spectrum)
        {
            var commonParameters = new CommonParameters();
            double precursorMz = matched.MonoisotopicMass.ToMz(1);
            var scan = new MsDataScan(spectrum.MassSpectrum, 1, 2, true, Polarity.Positive, 1,
                new MzRange(50, 2000), null, MZAnalyzerType.Orbitrap, spectrum.MassSpectrum.Size, null, null, "scan=1", precursorMz, 1, 1, precursorMz, 2,
                DissociationType.HCD, 0, precursorMz);
            var psm = new PeptideSpectralMatch(matched, 0, 0, 0, new Ms2ScanWithSpecificMass(scan, precursorMz, 1, "filepath", commonParameters),
                commonParameters, new List<MatchedFragmentIon>());
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            psm.SetMs2Scan(scan);
            return psm;
        }

        // The score the search gives a peptidoform on a spectrum: its HCD fragments, neutral losses included.
        private static double SearchScore(IBioPolymerWithSetMods peptidoform, MsDataScan spectrum)
        {
            var commonParameters = new CommonParameters();
            var products = new List<Product>();
            peptidoform.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            var scan = new Ms2ScanWithSpecificMass(spectrum, peptidoform.MonoisotopicMass.ToMz(1), 1, "filepath", commonParameters);
            return MetaMorpheusEngine.CalculatePeptideScore(spectrum, MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters));
        }

        private static Modification UniProtMod(string idWithMotif) =>
            GlobalVariables.AllModsKnown.First(m => m.ModificationType == "UniProt" && m.IdWithMotif == idWithMotif);

        private static Modification GptmdMod(string idWithMotif) =>
            GlobalVariables.AllModsKnown.First(m => m.ModificationType == "Common Biological" && m.IdWithMotif == idWithMotif);

        private static List<Modification> PhosphoMods() => GlobalVariables.AllModsKnown
            .Where(p => p.ModificationType == "Common Biological" && p.IdWithMotif.StartsWith("Phosphorylation on ")).ToList();

        private static List<Modification> DefaultGptmdMods()
        {
            var defaults = new GptmdParameters().ListOfModsGptmd;
            return GlobalVariables.AllModsKnown.Where(m => defaults.Contains((m.ModificationType, m.IdWithMotif))).ToList();
        }

        private static GptmdResults RunGptmdOnPhosphoPsm(List<SpectralMatch> psms, List<IGptmdFilter> filters = null, bool addAlternativeSites = true) =>
            RunGptmd(psms, PhosphoMods(), filters, addAlternativeSites);

        private static GptmdResults RunGptmd(List<SpectralMatch> psms, List<Modification> gptmdMods, List<IGptmdFilter> filters = null,
            bool addAlternativeSites = true, List<Modification> fixedMods = null)
        {
            return (GptmdResults)new GptmdEngine(psms, gptmdMods, new List<Tuple<double, double>>(),
                new Dictionary<string, Tolerance> { { "filepath", new PpmTolerance(10) } }, new CommonParameters(), null, new List<string>(), null,
                filters, addAlternativeSites, fixedMods).Run();
        }

        /// <summary>
        /// A modification that blocks the protease's cleavage cannot sit on a residue the protease cut after: trypsin
        /// does not cleave after an acetylated lysine, so a tryptic peptide cannot END in one. With the search's
        /// RespectCleavageBlockingModifications on, GPTMD must not localize such a modification there -- the same rule the
        /// search's digestion applies -- and must still try the peptide's other residues.
        ///
        /// The peptide is the protein's first nine residues: K4 is internal, K9 is the residue trypsin cut after, unless the
        /// protein ends there. Lys-N cuts BEFORE K, so its peptide ends at K9 only because K10 follows, and K9 is free.
        /// trypsin (K|) cuts before P too, so a K before P is refused like any other cut K; trypsin|P (K[P]|) never cuts
        /// there, so its K4 before P5 is no cleavage site. The spectrum is a stub, so every placement scores the same and
        /// GPTMD writes every site it allows.
        /// </summary>
        [Test]
        [TestCase("AAAKAAAAKAAAAR", "Acetylation", true, CleavageSpecificity.Full, "4", TestName = "BlockingModOnTheCutResidue_IsRefused")]
        [TestCase("AAAKAAAAKAAAAR", "Acetylation", false, CleavageSpecificity.Full, "4,9", TestName = "BlockingModOnTheCutResidue_IsKeptWithTheFlagOff")]
        [TestCase("AAAKAAAAK", "Acetylation", true, CleavageSpecificity.Full, "4,9", TestName = "BlockingModOnTheProteinCTerminus_IsKept")]
        [TestCase("AAAKAAAAKAAAAR", "Acetylation", true, CleavageSpecificity.Semi, "4,9", TestName = "BlockingModOnTheCutResidue_IsKeptInASemiSpecificSearch")]
        [TestCase("AAAKAAAAKAAAAR", "Methylation", true, CleavageSpecificity.Full, "4,9", TestName = "NonBlockingModOnTheCutResidue_IsKept")]
        [TestCase("AAAKAAAAKKAAAAR", "Acetylation", true, CleavageSpecificity.Full, "4,9", "Lys-N", TestName = "BlockingModOnTheLastResidue_IsKeptForAProteaseThatCutsBeforeIt")]
        [TestCase("AAAKAAAAKPAAAAR", "Acetylation", true, CleavageSpecificity.Full, "4", "trypsin", TestName = "BlockingModOnACutKBeforeP_IsRefused")]
        [TestCase("AAAKPAAAKAAAAR", "Acetylation", true, CleavageSpecificity.Full, "4", "trypsin|P", TestName = "BlockingModOnAKBeforePTheProteaseSkips_IsKept")]
        public static void GptmdRefusesACleavageBlockingModificationOnTheResidueTheProteaseCutAfter(string proteinSequence,
            string modificationName, bool respectCleavageBlockingModifications, CleavageSpecificity searchMode, string expectedSites,
            string protease = "trypsin")
        {
            ModificationMotif.TryGetMotif("K", out ModificationMotif motifK);
            double modificationMass = modificationName == "Acetylation" ? 42.010565 : 14.015650;
            var modification = new Modification(_originalId: modificationName, _modificationType: "Common Biological", _target: motifK,
                _locationRestriction: "Anywhere.", _monoisotopicMass: modificationMass);
            Assert.That(CleavageBlockingModifications.BlocksCleavageBy(modification, ProteaseDictionary.Dictionary["trypsin"]),
                Is.EqualTo(modificationName == "Acetylation"), "the test's premise: only acetyl blocks trypsin");

            var protein = new Protein(proteinSequence, "accession");
            var digestionParams = new DigestionParams(protease, minPeptideLength: 5, searchModeType: searchMode,
                respectCleavageBlockingModifications: respectCleavageBlockingModifications);
            var commonParameters = new CommonParameters(digestionParams: digestionParams);
            var peptide = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.BaseSequence == proteinSequence.Substring(0, 9) && p.OneBasedStartResidue == 1);

            MsDataScan stubScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null,
                MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            var scan = new Ms2ScanWithSpecificMass(stubScan,
                (new Proteomics.AminoAcidPolymer.Peptide(peptide.BaseSequence).MonoisotopicMass + modificationMass).ToMz(1), 1, "filepath", new CommonParameters());
            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            psm.SetMs2Scan(scan.TheScan);

            var engine = new GptmdEngine(new List<SpectralMatch> { psm }, new List<Modification> { modification }, new List<Tuple<double, double>>(),
                new Dictionary<string, Tolerance> { { "filepath", new PpmTolerance(10) } }, commonParameters, null, new List<string>(), null);
            var result = (GptmdResults)engine.Run();

            Assert.That(result.Mods["accession"].Select(site => site.Item1).OrderBy(site => site).ToArray(),
                Is.EqualTo(expectedSites.Split(',').Select(int.Parse).ToArray()));
        }

        /// <summary>
        /// What the stub spectrum above cannot show: an internal site written ONLY because the best-scoring residue was
        /// refused. Mods.txt's own Acetylation on K, and a spectrum of the b and y ions of AAAKAAAAK with the acetyl on
        /// K9, so K9 matches all 16 ions and K4 matches 6. Flag off, GPTMD writes K9 alone; flag on, K9 is refused and
        /// K4 is the best placement left.
        /// </summary>
        [Test]
        [TestCase(true, "4", TestName = "RealAcetylLysine_GoesToTheNextBestSiteWhenTheBestBlocksTheCut")]
        [TestCase(false, "9", TestName = "RealAcetylLysine_GoesToTheSpectrumsBestSiteWithTheFlagOff")]
        public static void GptmdWritesTheNextBestSiteWhenTheBestOneBlocksTheCleavage(bool respectCleavageBlockingModifications, string expectedSites)
        {
            Modification acetylK = GlobalVariables.AllModsKnown.Single(m => m.ModificationType == "Common Biological" && m.IdWithMotif == "Acetylation on K");
            Assert.That(CleavageBlockingModifications.BlocksCleavageBy(acetylK, ProteaseDictionary.Dictionary["trypsin"]), Is.True);

            var protein = new Protein("AAAKAAAAKAAAAR", "accession");
            var digestionParams = new DigestionParams(minPeptideLength: 5, respectCleavageBlockingModifications: respectCleavageBlockingModifications);
            var commonParameters = new CommonParameters(digestionParams: digestionParams);
            var peptide = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .First(p => p.BaseSequence == "AAAKAAAAK" && p.OneBasedStartResidue == 1);

            var acetylOnK9 = peptide.Localize(8, acetylK.MonoisotopicMass.Value);
            var products = new List<Product>();
            acetylOnK9.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            double[] mzs = products.Select(p => p.NeutralMass.ToMz(1)).OrderBy(mz => mz).ToArray();
            MsDataScan spectrum = new MsDataScan(new MzSpectrum(mzs, mzs.Select(_ => 1.0).ToArray(), false), 0, 1, true, Polarity.Positive, double.NaN, null, null,
                MZAnalyzerType.Orbitrap, mzs.Length, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(spectrum, acetylOnK9.MonoisotopicMass.ToMz(1), 1, "filepath", new CommonParameters());
            SpectralMatch psm = new PeptideSpectralMatch(peptide, 0, 0, 0, scan, commonParameters, new List<MatchedFragmentIon>());
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            psm.SetMs2Scan(scan.TheScan);

            var engine = new GptmdEngine(new List<SpectralMatch> { psm }, new List<Modification> { acetylK }, new List<Tuple<double, double>>(),
                new Dictionary<string, Tolerance> { { "filepath", new PpmTolerance(10) } }, commonParameters, null, new List<string>(), null);
            var result = (GptmdResults)engine.Run();

            Assert.That(result.Mods["accession"].Select(site => site.Item1).OrderBy(site => site).ToArray(),
                Is.EqualTo(expectedSites.Split(',').Select(int.Parse).ToArray()));
        }

        private sealed class RecordingGptmdFilter : IGptmdFilter
        {
            public List<int> Sites { get; } = new();

            public bool Passes(
                IBioPolymerWithSetMods candidatePeptide,
                SpectralMatch psm,
                double newScore,
                double originalScore,
                List<MatchedFragmentIon> matchedIons,
                int peptideOneBasedModSite,
                int peptideLength,
                Modification modAttemptingToAdd)
            {
                Sites.Add(peptideOneBasedModSite);
                return true;
            }
        }
    }
}
