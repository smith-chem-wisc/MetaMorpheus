using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace Test
{
    [TestFixture]
    public static class PsvTsvTest
    {
        [Test]
        public static void ReadOGlycoSinglePsms()
        {
            string psmFile = @"TestData\oglycoSinglePsms.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(2, parsedPsms.Count);
        }

        [Test]
        public static void ReadOGlycoPsms()
        {
            string psmFile = @"TestData\oglyco.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(9, parsedPsms.Count);
        }

        [Test]
        public static void ReadExcelEditedPsms()
        {
            string psmFile = @"TestData\ExcelEditedPeptide.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(1, parsedPsms.Count);
            IEnumerable<string> expectedIons = new string[] { "y3+1", "y4+1", "b4+1", "b5+1", "b6+1", "b8+1" };
            Assert.That(6 == parsedPsms[0].MatchedIons.Select(p => p.Annotation).Intersect(expectedIons).Count());
            Assert.That("TADDYTWEGDVGNDNAYQKFVK", Is.EqualTo(parsedPsms[0].FullSequence));
        }

        [Test]
        public static void MetaDrawLogicTestOglyco()
        {
            var metadrawLogic = new MetaDrawLogic();
            string psmFile = @"TestData\oglyco.psmtsv";
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(false, true);

            Assert.That(!errors.Any());
        }

        [Test]
        public static void CrosslinkPsmFromTsvTest()
        {
            string psmFile = @"XlTestData\XL_Intralinks_MIons.tsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(1, parsedPsms.Count);
            Assert.That(parsedPsms[0].UniqueSequence, Is.EqualTo("EKVLTSSAR(2)SLGKVGTR(4)"));
        }

        [Test]
        public static void CrosslinkPsmFromTsvToLibrarySpectrumTest()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\XL_Intralinks_MIons.tsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = PsmTsvReader.ReadTsv(psmTsvPath, out warnings).ToList();
            Assert.That(warnings.Count == 0);

            CrosslinkLibrarySpectrum librarySpectrum = psms[0].ToLibrarySpectrum() as CrosslinkLibrarySpectrum;
            Assert.IsNotNull(librarySpectrum);
            Assert.AreEqual("Name: EKVLTSSAR(2)SLGKVGTR(4)/4", librarySpectrum.ToString().Split('\n')[0].Trim());

            // This test would be better if MatchedIon.equals method worked, but it breaks because the mz comparison is implemented incorrectly.
            CollectionAssert.AreEquivalent(librarySpectrum.MatchedFragmentIons.Select(ion => ion.Annotation), psms[0].MatchedIons.Select(ion => ion.Annotation));
            CollectionAssert.AreEquivalent(librarySpectrum.BetaPeptideSpectrum.MatchedFragmentIons.Select(ion => ion.Annotation), psms[0].BetaPeptideMatchedIons.Select(ion => ion.Annotation));
        }

        [Test]
        public static void TestPsmFromTsvDisambiguatingConstructor()
        {
            // initialize values
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = PsmTsvReader.ReadTsv(psmTsvPath, out warnings);
            PsmFromTsv psm = psms.First();

            // non ambiguous construction should not be successful
            string fullSeq = psm.FullSequence;
            fullSeq = fullSeq.Substring(0, fullSeq.Length - 1);
            PsmFromTsv modifiedPsm = new(psm, fullSeq);
            Assert.That(modifiedPsm.FullSequence == fullSeq);

            // disambiguation construction
            var ambiguousPsms = psms.Where(p => p.FullSequence.Contains('|'));
            PsmFromTsv ambiguousPsm = ambiguousPsms.First();
            var fullSeqStrings = ambiguousPsm.FullSequence.Split('|');

            PsmFromTsv modifiedAmbiguousPsm = new(ambiguousPsm, fullSeqStrings[0]);
            List<string[]> test = new();
            foreach (var ambPsm in ambiguousPsms)
            {
                PsmFromTsv disambiguatedPSM = new(ambPsm, ambPsm.FullSequence.Split("|")[0]);
                Assert.That(disambiguatedPSM.StartAndEndResiduesInProtein == ambPsm.StartAndEndResiduesInProtein.Split("|")[0]);
                Assert.That(disambiguatedPSM.BaseSeq == ambPsm.BaseSeq.Split("|")[0]);
                Assert.That(disambiguatedPSM.EssentialSeq == ambPsm.EssentialSeq.Split("|")[0]);
                Assert.That(disambiguatedPSM.ProteinAccession == ambPsm.ProteinAccession.Split("|")[0]);
                Assert.That(disambiguatedPSM.PeptideMonoMass == ambPsm.PeptideMonoMass.Split("|")[0]);
                Assert.That(disambiguatedPSM.MassDiffDa == ambPsm.MassDiffDa.Split("|")[0]);
                Assert.That(disambiguatedPSM.MassDiffPpm == ambPsm.MassDiffPpm.Split("|")[0]);
                Assert.That(disambiguatedPSM.ProteinName == ambPsm.ProteinName.Split("|")[0]);
                Assert.That(disambiguatedPSM.GeneName == ambPsm.GeneName.Split("|")[0]);

                for (int i = 0; i < ambPsm.MatchedIons.Count; i++)
                {
                    Assert.That(disambiguatedPSM.MatchedIons[i] == ambPsm.MatchedIons[i]);
                }

                if (ambPsm.StartAndEndResiduesInProtein.Split("|").Count() > 1)
                {
                    for (int i = 1; i < ambPsm.StartAndEndResiduesInProtein.Split("|").Count(); i++)
                    {
                        disambiguatedPSM = new(ambPsm, ambPsm.FullSequence.Split("|")[i], i);
                        Assert.That(disambiguatedPSM.StartAndEndResiduesInProtein == ambPsm.StartAndEndResiduesInProtein.Split("|")[i]);
                        Assert.That(disambiguatedPSM.BaseSeq == ambPsm.BaseSeq.Split("|")[i]);
                        Assert.That(disambiguatedPSM.EssentialSeq == ambPsm.EssentialSeq.Split("|")[i]);
                        Assert.That(disambiguatedPSM.ProteinAccession == ambPsm.ProteinAccession.Split("|")[i]);
                        Assert.That(disambiguatedPSM.ProteinName == ambPsm.ProteinName.Split("|")[i]);
                        Assert.That(disambiguatedPSM.GeneName == ambPsm.GeneName.Split("|")[i]);

                        if (ambPsm.PeptideMonoMass.Split("|").Count() == 1)
                        {
                            Assert.That(disambiguatedPSM.PeptideMonoMass == ambPsm.PeptideMonoMass.Split("|")[0]);
                            Assert.That(disambiguatedPSM.MassDiffDa == ambPsm.MassDiffDa.Split("|")[0]);
                            Assert.That(disambiguatedPSM.MassDiffPpm == ambPsm.MassDiffPpm.Split("|")[0]);
                        }
                        else
                        {
                            Assert.That(disambiguatedPSM.PeptideMonoMass == ambPsm.PeptideMonoMass.Split("|")[i]);
                            Assert.That(disambiguatedPSM.MassDiffDa == ambPsm.MassDiffDa.Split("|")[i]);
                            Assert.That(disambiguatedPSM.MassDiffPpm == ambPsm.MassDiffPpm.Split("|")[i]);
                        }
                    }
                }
            }
        }

        [Test]
        public static void TestParseModification()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = PsmTsvReader.ReadTsv(psmTsvPath, out warnings).Take(20).ToList();
            Assert.That(warnings.Count == 0);

            // psm with single modificaiton
            PsmFromTsv singleMod = psms[0];
            var modDict = PsmFromTsv.ParseModifications(singleMod.FullSequence);
            Assert.That(modDict.Count == 1);
            Assert.That(modDict.ContainsKey(37));
            Assert.That(modDict.Values.First().Contains("Common Fixed:Carbamidomethyl on C"));

            // psm with two modifications
            PsmFromTsv twoMods = psms[15];
            modDict = PsmFromTsv.ParseModifications(twoMods.FullSequence);
            Assert.That(modDict.Count == 2);
            Assert.That(modDict.ContainsKey(0) && modDict.ContainsKey(104));
            Assert.That(modDict[0].Count == 1);
            Assert.That(modDict[0].Contains("UniProt:N-acetylserine on S"));
            Assert.That(modDict[104].Count == 1);
            Assert.That(modDict[104].Contains("UniProt:N5-methylglutamine on Q"));


            // psm with two mods on the same amino acid
            string fullSeq = "[Common Fixed:Carbamidomethyl on C]|[UniProt:N-acetylserine on S]KPRKIEEIKDFLLTARRKDAKSVKIKKNKDNVKFK";
            modDict = PsmFromTsv.ParseModifications(fullSeq);
            Assert.That(modDict.Count == 1);
            Assert.That(modDict.ContainsKey(0));
            Assert.That(modDict[0].Count == 2);
            Assert.That(modDict[0].Contains("Common Fixed:Carbamidomethyl on C"));
            Assert.That(modDict[0].Contains("UniProt:N-acetylserine on S"));
        }

        [Test]
        public static void TestRemoveSpecialCharacters()
        {
            // successful removal of the default character
            string toRemove = "ANDVHAO|CNVASDF|ABVCUAE";
            int length = toRemove.Length;
            PsmFromTsv.RemoveSpecialCharacters(ref toRemove);
            Assert.That(toRemove.Length == length - 2);
            Assert.That(toRemove.Equals("ANDVHAOCNVASDFABVCUAE"));

            // does not remove default character when prompted otherwise
            toRemove = "ANDVHAO|CNVASDF|ABVCUAE";
            PsmFromTsv.RemoveSpecialCharacters(ref toRemove, specialCharacter: @"\[");
            Assert.That(toRemove.Length == length);
            Assert.That(toRemove.Equals("ANDVHAO|CNVASDF|ABVCUAE"));

            // replaces default symbol when prompted
            PsmFromTsv.RemoveSpecialCharacters(ref toRemove, replacement: @"%");
            Assert.That(toRemove.Length == length);
            Assert.That(toRemove.Equals("ANDVHAO%CNVASDF%ABVCUAE"));

            // replaces inputted symbol with non-default symbol
            PsmFromTsv.RemoveSpecialCharacters(ref toRemove, replacement: @"=", specialCharacter: @"%");
            Assert.That(toRemove.Length == length);
            Assert.That(toRemove.Equals("ANDVHAO=CNVASDF=ABVCUAE"));
        }

        [Test]
        public static void TestToString()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = PsmTsvReader.ReadTsv(psmTsvPath, out warnings).Take(3).ToList();
            Assert.That(warnings.Count == 0);

            Assert.That(psms[0].FullSequence.Equals(psms[0].ToString()));
            Assert.That(psms[1].FullSequence.Equals(psms[1].ToString()));
            Assert.That(psms[2].FullSequence.Equals(psms[2].ToString()));
        }

        [Test]
        public static void TestSimpleToLibrarySpectrum() 
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchResults.psmtsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = PsmTsvReader.ReadTsv(psmTsvPath, out warnings).Take(3).ToList();
            Assert.That(warnings.Count == 0);

            string librarySpectrum = psms[0].ToLibrarySpectrum().ToString();

            string expectedLibrarySpectrum = File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\simple.msp"));

            //not a great way to test equality but we are experiencing a great deal of 10th digit rounding differences
            Assert.AreEqual(Regex.Matches(expectedLibrarySpectrum, "ppm").Count, Regex.Matches(librarySpectrum, "ppm").Count);

            
            //the code below tests the addition and correct output for neutral loss fragments
            Product p = new Product(ProductType.bWaterLoss, FragmentationTerminus.N, 1, 1, 1, 18);
            MatchedFragmentIon matchedIon = new(p, 1, 1, 1);
            psms[0].MatchedIons.Add(matchedIon);
            string librarySpectrumWithNeutralLoss = psms[0].ToLibrarySpectrum().ToString();

            Assert.That(librarySpectrumWithNeutralLoss.Contains("WaterLoss"));
        }
        [Test]
        public static void TestPsmSortFunction()
        {
            string[] sequences = {
                "ABCKPEPR",
                "BRPEPR",
                "ARPEPR",
                "PEPPER",
                "PEPTIDE",
                "PRTIEN"
            };

            var p = new List<Protein>();
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();
            for (int i = 0; i < sequences.Length; i++)
            {
                p.Add(new Protein(sequences[i], (i + 1).ToString(), null, gn, new Dictionary<int, List<Modification>>()));
            }

            CommonParameters commonParams = new CommonParameters(digestionParams: new DigestionParams(protease: "trypsin", minPeptideLength: 1));

            PeptideWithSetModifications pepOne = new PeptideWithSetModifications(protein: p.ElementAt(0), digestionParams: commonParams.DigestionParams, oneBasedStartResidueInProtein: 5, oneBasedEndResidueInProtein: 8, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ABCKPEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepTwo = new PeptideWithSetModifications(protein: p.ElementAt(1), digestionParams: commonParams.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "BRPEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepThree = new PeptideWithSetModifications(protein: p.ElementAt(2), digestionParams: commonParams.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "ARPEPR", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepFour = new PeptideWithSetModifications(protein: p.ElementAt(3), digestionParams: commonParams.DigestionParams, oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 3, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPPER", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepFive = new PeptideWithSetModifications(protein: p.ElementAt(4), digestionParams: commonParams.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PEPTIDE", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);
            PeptideWithSetModifications pepSix = new PeptideWithSetModifications(protein: p.ElementAt(5), digestionParams: commonParams.DigestionParams, oneBasedStartResidueInProtein: 3, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "PRTIEN", missedCleavages: 0, allModsOneIsNterminus: new Dictionary<int, Modification>(), numFixedMods: 0);

            MsDataScan scanNumberOne = new MsDataScan(new MzSpectrum(new double[] { 10 }, new double[] { 1 }, false), 1, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", 10, 2, 100, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            MsDataScan scanNumberTwo = new MsDataScan(new MzSpectrum(new double[] { 20 }, new double[] { 1 }, false), 2, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=2", 20, 2, 100, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            MsDataScan scanNumberThree = new MsDataScan(new MzSpectrum(new double[] { 20 }, new double[] { 1 }, false), 3, 2, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=3", 20, 2, 100, double.NaN, null, DissociationType.AnyActivationType, 0, null);

            Ms2ScanWithSpecificMass ms2ScanOneMzTen = new Ms2ScanWithSpecificMass(scanNumberOne, 10, 2, "File", new CommonParameters());
            Ms2ScanWithSpecificMass ms2ScanOneMzTwenty = new Ms2ScanWithSpecificMass(scanNumberTwo, 20, 2, "File", new CommonParameters());
            Ms2ScanWithSpecificMass ms2ScanTwoMzTwenty = new Ms2ScanWithSpecificMass(scanNumberTwo, 20, 2, "File", new CommonParameters());
            Ms2ScanWithSpecificMass ms2ScanThreeMzTwenty = new Ms2ScanWithSpecificMass(scanNumberThree, 20, 2, "File", new CommonParameters());

            //highest score
            PeptideSpectralMatch psmOne = new(pepOne, 0, 10, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());
            
            //second highest score 9. delta score 1
            PeptideSpectralMatch psmTwo = new(pepOne, 0, 8, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());
            psmTwo.AddOrReplace(pepTwo, 9, 0, true, new List<MatchedFragmentIon>(), 0);

            //second highest score 9. delta score 0.1
            PeptideSpectralMatch psmThree = new(pepOne, 0, 8.90000000000000000000000000000, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());
            psmThree.AddOrReplace(pepTwo, 9, 0, true, new List<MatchedFragmentIon>(), 0);

            //third highest score delta score 1. ppm error -888657.54 low
            PeptideSpectralMatch psmFour = new(pepOne, 0, 7, 0, ms2ScanOneMzTwenty, commonParams,
                new List<MatchedFragmentIon>());
            psmFour.AddOrReplace(pepFour, 8, 0, true, new List<MatchedFragmentIon>(), 0);

            //third highest score 8. delta score 1. ppm error -947281.29 high
            PeptideSpectralMatch psmFive = new(pepOne, 0, 7, 0, ms2ScanOneMzTen, commonParams,
                new List<MatchedFragmentIon>());
            psmFive.AddOrReplace(pepFour, 8, 0, true, new List<MatchedFragmentIon>(), 0);

            //fourth highest score 7. delta score 1. same ppm error
            PeptideSpectralMatch psmSix = new(pepOne, 0, 6, 0, ms2ScanTwoMzTwenty, commonParams,
                new List<MatchedFragmentIon>());
            psmSix.AddOrReplace(pepSix, 7, 0, true, new List<MatchedFragmentIon>(), 0);

            //fourth highest score 7. delta score 1. same ppm error
            PeptideSpectralMatch psmSeven = new(pepOne, 0, 6, 0, ms2ScanThreeMzTwenty, commonParams,
                new List<MatchedFragmentIon>());
            psmSeven.AddOrReplace(pepSix, 7, 0, true, new List<MatchedFragmentIon>(), 0);

            List<PeptideSpectralMatch> psms = new List<PeptideSpectralMatch> { psmFour, psmOne, psmThree, psmSeven, psmTwo, psmFive, psmSix };
            psms.ForEach(j => j.ResolveAllAmbiguities());
            psms.ForEach(j => j.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0));

            List<PeptideSpectralMatch> orderedPsms = psms.OrderByDescending(p=>p).ToList();

            Assert.AreEqual(10, orderedPsms[0].Score);
            Assert.AreEqual(9, orderedPsms[1].Score);
            Assert.AreEqual(9, orderedPsms[2].Score);
            Assert.AreEqual(8, orderedPsms[3].Score);
            Assert.AreEqual(8, orderedPsms[4].Score);
            Assert.AreEqual(7, orderedPsms[5].Score);
            Assert.AreEqual(7, orderedPsms[6].Score);

            Assert.AreEqual(5, orderedPsms[0].RunnerUpScore);
            Assert.AreEqual(8.9, orderedPsms[1].RunnerUpScore);
            Assert.AreEqual(8, orderedPsms[2].RunnerUpScore);
            Assert.AreEqual(7, orderedPsms[3].RunnerUpScore);
            Assert.AreEqual(7, orderedPsms[4].RunnerUpScore);
            Assert.AreEqual(6, orderedPsms[5].RunnerUpScore);
            Assert.AreEqual(6, orderedPsms[6].RunnerUpScore);

            Assert.IsTrue(Math.Abs(orderedPsms[3].PrecursorMassErrorPpm.First()) < Math.Abs(orderedPsms[4].PrecursorMassErrorPpm.First()));
            Assert.IsTrue(orderedPsms[5].ScanNumber < orderedPsms[6].ScanNumber);
        }
    }
}