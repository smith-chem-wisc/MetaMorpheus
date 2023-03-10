using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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
    }
}