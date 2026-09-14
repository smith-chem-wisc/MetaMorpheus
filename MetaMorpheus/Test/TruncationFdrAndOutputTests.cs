using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.Truncation;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    /// <summary>
    /// Phase 3 tests: pooled FDR/PEP over the inherited-intact + truncation PSMs (#15, #18).
    /// Output-file tests are added alongside as that step lands.
    /// </summary>
    [TestFixture]
    public class TruncationFdrAndOutputTests
    {
        private CommonParameters _cp;

        [OneTimeSetUp]
        public void Setup()
        {
            _cp = new CommonParameters(
                dissociationType: DissociationType.HCD,
                precursorMassTolerance: new PpmTolerance(10),
                productMassTolerance: new PpmTolerance(20));
        }

        [Test]
        public void PooledFdr_AssignsQValues()
        {
            // Deliberate target/decoy ordering by score (desc): T T T D T D.
            var pool = new List<SpectralMatch>
            {
                MakeTargetPsm("PEPTIDEAAAR", 10.0, 1),
                MakeTargetPsm("PEPTIDEBBBR", 9.0, 2),
                MakeTargetPsm("PEPTIDECCCR", 8.0, 3),
                MakeDecoyPsm("PEPTIDEDDDR", 7.0, 4),
                MakeTargetPsm("PEPTIDEEEER", 6.0, 5),
                MakeDecoyPsm("PEPTIDEFFFR", 5.0, 6)
            };

            List<SpectralMatch> result = TruncationFdr.RunPooledFdr(pool, _cp, 1,
                new List<(string, CommonParameters)> { ("synthetic", _cp) }, "TruncationTask", doPep: false);

            var byScoreDesc = result.OrderByDescending(p => p.Score).ToList();

            // q-values are populated, within [0,1], and monotonically non-decreasing as score decreases.
            Assert.That(byScoreDesc.All(p => p.FdrInfo.QValue >= 0 && p.FdrInfo.QValue <= 1), Is.True);
            for (int i = 1; i < byScoreDesc.Count; i++)
            {
                Assert.That(byScoreDesc[i].FdrInfo.QValue, Is.GreaterThanOrEqualTo(byScoreDesc[i - 1].FdrInfo.QValue - 1e-9));
            }
            // The top, decoy-free target has q-value 0.
            Assert.That(byScoreDesc.First().FdrInfo.QValue, Is.EqualTo(0).Within(1e-9));
            Assert.That(byScoreDesc.First().IsDecoy, Is.False);
        }

        [Test]
        public void PepRuns_OrSkipsCleanly_OnSmallPool()
        {
            var pool = new List<SpectralMatch>
            {
                MakeTargetPsm("PEPTIDEAAAR", 10.0, 1),
                MakeDecoyPsm("PEPTIDEBBBR", 5.0, 2)
            };

            // doPep:true on a tiny pool must not throw (PEP needs >1000 PSMs/protease to train) and q-values still compute.
            List<SpectralMatch> result = null;
            Assert.DoesNotThrow(() => result = TruncationFdr.RunPooledFdr(pool, _cp, 1,
                new List<(string, CommonParameters)> { ("synthetic", _cp) }, "TruncationTask", doPep: true));

            Assert.That(result.All(p => p.FdrInfo.QValue >= 0 && p.FdrInfo.QValue <= 1), Is.True);

            // Pin the branch, not just absence-of-throw: a 2-PSM pool can't train PEP, so PEP must be
            // SKIPPED -- every row keeps the "PEP not computed" sentinel (PEP_QValue == 2), rather than
            // silently emitting an untrained PEP q-value (#18).
            Assert.That(result.All(p => p.FdrInfo.PEP_QValue == 2), Is.True);
        }

        [Test]
        public void OutputFiles_HaveExpectedColumns()
        {
            var psms = new List<SpectralMatch>
            {
                MakePsm(MakeProteoform("PEPTIDEKAAAA", "ACC1", "C-terminal truncation(1-8)"), 10.0, 1),
                MakePsm(MakeProteoform("PEPTIDERBBBB", "ACC2", "N-terminal truncation(3-12)"), 9.0, 2)
            };
            var processed = TruncationFdr.RunPooledFdr(psms, _cp, 1,
                new List<(string, CommonParameters)> { ("synthetic", _cp) }, "TruncationTask", doPep: false);

            string dir = TestContext.CurrentContext.TestDirectory;
            string psmsPath = Path.Combine(dir, "AllTruncatedPSMs_test.psmtsv");
            string proteoformsPath = Path.Combine(dir, "AllTruncatedProteoforms_test.psmtsv");
            try
            {
                TruncationOutput.WritePsms(processed, psmsPath);
                TruncationOutput.WriteProteoforms(processed, proteoformsPath);

                string[] psmLines = File.ReadAllLines(psmsPath);
                string[] header = psmLines[0].Split('\t');
                Assert.That(header, Does.Contain("Full Sequence"));
                Assert.That(header, Does.Contain("Description"));
                Assert.That(header, Does.Contain("QValue"));
                Assert.That(header, Does.Contain("Decoy/Contaminant/Target"));
                Assert.That(psmLines.Length - 1, Is.EqualTo(2)); // one row per PSM

                int descIdx = Array.IndexOf(header, "Description");
                Assert.That(psmLines.Skip(1).Any(l => l.Split('\t')[descIdx].Contains("truncation(")), Is.True);

                string[] proteoformLines = File.ReadAllLines(proteoformsPath);
                Assert.That(proteoformLines.Length - 1, Is.EqualTo(2)); // two distinct truncated full sequences
            }
            finally
            {
                File.Delete(psmsPath);
                File.Delete(proteoformsPath);
            }
        }

        [Test]
        public void IntactInheritedRow_HasFullLengthLabel()
        {
            var psms = new List<SpectralMatch> { MakePsm(MakeProteoform("PEPTIDEKAAAA", "ACC1", "full-length"), 10.0, 1) };
            var processed = TruncationFdr.RunPooledFdr(psms, _cp, 1,
                new List<(string, CommonParameters)> { ("synthetic", _cp) }, "TruncationTask", doPep: false);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "AllTruncatedPSMs_fulllen_test.psmtsv");
            try
            {
                TruncationOutput.WritePsms(processed, path);
                string[] lines = File.ReadAllLines(path);
                int descIdx = Array.IndexOf(lines[0].Split('\t'), "Description");
                Assert.That(lines[1].Split('\t')[descIdx], Is.EqualTo("full-length"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void ProteinAccessionAmbiguity_PreservedInOutput()
        {
            Ms2ScanWithSpecificMass scan = BuildMinimalScan(1, MakeProteoform("PEPTIDE", "x").MonoisotopicMass, _cp);
            var pA1 = MakeProteoform("PEPTIDE", "A1", "C-terminal truncation(1-7)");
            var pA2 = MakeProteoform("PEPTIDE", "A2", "C-terminal truncation(1-7)");
            var psmA1 = new PeptideSpectralMatch(pA1, 0, 10.0, 0, scan, _cp, new List<MatchedFragmentIon>());
            var psmA2 = new PeptideSpectralMatch(pA2, 0, 10.0, 0, scan, _cp, new List<MatchedFragmentIon>());

            var collapsed = TruncationPass3.CollapseDuplicateTruncations(new[]
            {
                new TruncationPsm { SpectralMatch = psmA1, TruncatedForm = pA1, ScanIndex = 0, ProteinAccessions = "A1", TruncationProductType = TruncationPass3.CTerminalTruncation },
                new TruncationPsm { SpectralMatch = psmA2, TruncatedForm = pA2, ScanIndex = 0, ProteinAccessions = "A2", TruncationProductType = TruncationPass3.CTerminalTruncation }
            });
            Assert.That(collapsed.Count, Is.EqualTo(1));

            var processed = TruncationFdr.RunPooledFdr(collapsed, _cp, 1,
                new List<(string, CommonParameters)> { ("synthetic", _cp) }, "TruncationTask", doPep: false);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "AllTruncatedProteoforms_amb_test.psmtsv");
            try
            {
                TruncationOutput.WriteProteoforms(processed, path);
                string[] lines = File.ReadAllLines(path);
                int accIdx = Array.IndexOf(lines[0].Split('\t'), "Accession");
                Assert.That(lines.Length - 1, Is.EqualTo(1)); // one proteoform row
                string accCell = lines[1].Split('\t')[accIdx];
                Assert.That(accCell, Does.Contain("A1"));
                Assert.That(accCell, Does.Contain("A2"));
                Assert.That(accCell, Does.Contain("|"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        // ---------- helpers ----------

        private SpectralMatch MakeTargetPsm(string sequence, double score, int scanNumber)
        {
            return MakePsm(MakeProteoform(sequence, "ACC_" + sequence), score, scanNumber);
        }

        private SpectralMatch MakeDecoyPsm(string sequence, double score, int scanNumber)
        {
            var target = MakeProteoform(sequence, "ACC_" + sequence);
            var decoy = target.GetReverseDecoyFromTarget(new int[target.BaseSequence.Length]);
            return MakePsm(decoy, score, scanNumber);
        }

        private SpectralMatch MakePsm(PeptideWithSetModifications peptide, double score, int scanNumber)
        {
            Ms2ScanWithSpecificMass scan = BuildMinimalScan(scanNumber, peptide.MonoisotopicMass, _cp);
            return new PeptideSpectralMatch(peptide, 0, score, scanNumber - 1, scan, _cp, new List<MatchedFragmentIon>());
        }

        private static PeptideWithSetModifications MakeProteoform(string sequence, string accession, string peptideDescription = "top-down")
        {
            var protein = new Protein(sequence, accession);
            var digestionParams = new DigestionParams(protease: "top-down", minPeptideLength: 1, maxPeptideLength: 100000);
            return new PeptideWithSetModifications(protein, digestionParams, 1, sequence.Length,
                CleavageSpecificity.Full, peptideDescription, 0, new Dictionary<int, Modification>(), 0);
        }

        private static Ms2ScanWithSpecificMass BuildMinimalScan(int scanNumber, double precursorMass, CommonParameters cp)
        {
            var spectrum = new MzSpectrum(new[] { 500.0 }, new[] { 1000.0 }, false);
            var msDataScan = new MsDataScan(
                massSpectrum: spectrum, oneBasedScanNumber: scanNumber, msnOrder: 2, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: scanNumber, scanWindowRange: new MzRange(0, 1_000_000),
                scanFilter: "f", mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 1000.0,
                injectionTime: 1.0, noiseData: null, nativeId: $"scan={scanNumber}");

            var envelope = new IsotopicEnvelope(new List<(double mz, double intensity)> { (500.0, 1000.0) }, 499.0, 1, 1000.0, 0);
            return new Ms2ScanWithSpecificMass(msDataScan, precursorMass.ToMz(1), 1, "synthetic.mzML", cp,
                neutralExperimentalFragments: new[] { envelope });
        }
    }
}
