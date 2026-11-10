using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectrumMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// A spectral library MetaMorpheus writes stores each precursor at the theoretical m/z of the form it identified
    /// (that form's monoisotopic mass at the scan's charge), not at the m/z the scan measured. The scans here take the
    /// 13C peak for the monoisotopic peak, so the measured m/z is 1.003 Da / z high and a one-notch acceptor still
    /// identifies them.
    /// </summary>
    [TestFixture]
    public static class SpectralLibraryPrecursorMzTests
    {
        private const string OxidizedSequence = "PEPTM[Common Variable:Oxidation on M]IDEK";

        private static double ThirteenCPeakMz(IHasMass form, int charge) => (form.MonoisotopicMass + Constants.C13MinusC12).ToMz(charge);

        private static PeptideWithSetModifications Peptide(string fullSequence) =>
            new(fullSequence, GlobalVariables.AllModsKnownDictionary,
                p: new Protein(new PeptideWithSetModifications(fullSequence, GlobalVariables.AllModsKnownDictionary).BaseSequence, "accession"));

        /// <summary>
        /// A match to <paramref name="peptide"/> on a scan whose precursor is the peptide's 13C peak (notch 1).
        /// </summary>
        private static PeptideSpectralMatch ThirteenCPickedMatch(PeptideWithSetModifications peptide, int charge)
        {
            var commonParameters = new CommonParameters();
            MsDataScan ms2 = new TestDataFile(peptide).GetOneBasedScan(2);
            var scan = new Ms2ScanWithSpecificMass(ms2, ThirteenCPeakMz(peptide, charge), charge, "file", commonParameters);

            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters);

            var match = new PeptideSpectralMatch(peptide, 1, MetaMorpheusEngine.CalculatePeptideScore(ms2, matchedIons), 0, scan, commonParameters, matchedIons);
            match.ResolveAllAmbiguities();
            return match;
        }

        [Test]
        [TestCase("PEPTIDEK", 2)]
        [TestCase(OxidizedSequence, 3)]
        public static void LibrarySpectrumStoresTheTheoreticalPrecursorMzNotTheMeasuredOne(string fullSequence, int charge)
        {
            PeptideWithSetModifications peptide = Peptide(fullSequence);
            PeptideSpectralMatch match = ThirteenCPickedMatch(peptide, charge);

            LibrarySpectrum spectrum = match.ToLibrarySpectrum();

            Assert.That(spectrum.PrecursorMz, Is.EqualTo(peptide.MonoisotopicMass.ToMz(charge)).Within(1e-9));
            Assert.That(match.ScanPrecursorMonoisotopicPeakMz - spectrum.PrecursorMz, Is.EqualTo(Constants.C13MinusC12 / charge).Within(1e-9));
            Assert.That(spectrum.Sequence, Is.EqualTo(fullSequence));
            Assert.That(spectrum.ChargeState, Is.EqualTo(charge));
            Assert.That(spectrum.MatchedFragmentIons, Is.EqualTo(match.MatchedFragmentIons));
            Assert.That(spectrum.RetentionTime, Is.EqualTo(match.ScanRetentionTime));
            Assert.That(spectrum.IsDecoy, Is.False);
        }

        [Test]
        public static void LibrarySpectrumIsNullWhenTheCandidatesDifferInMass()
        {
            PeptideSpectralMatch match = ThirteenCPickedMatch(Peptide("PEPTMIDEK"), 2);
            match.AddOrReplace(Peptide(OxidizedSequence), match.Score, 0, true, match.MatchedFragmentIons);
            match.ResolveAllAmbiguities();

            Assert.That(match.BioPolymerWithSetModsMonoisotopicMass, Is.Null);
            Assert.That(match.ToLibrarySpectrum(), Is.Null);
        }

        /// <summary>
        /// Writes a library, then updates another, through the search task, from scans whose precursor is the 13C peak.
        /// Every spectrum either writes must carry the theoretical precursor m/z.
        /// </summary>
        [Test]
        public static void SearchTaskWritesAndUpdatesLibrariesAtTheTheoreticalPrecursorMz()
        {
            const int charge = 2;
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SpectralLibraryPrecursorMz");
            if (Directory.Exists(folder))
                Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);

            var proteins = new List<Protein> { new("PEPTIDEK", "P1"), new("SEQENTIALK", "P2") };
            var peptides = proteins.Select(p => new PeptideWithSetModifications(p.BaseSequence, new Dictionary<string, Modification>(), p: p)).ToList();
            string database = Path.Combine(folder, "proteins.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins, database);
            string spectra = Path.Combine(folder, "thirteenCPicks.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(ThirteenCPickedFile(peptides, charge), spectra, false);

            // write a library
            RunSearch(folder, "Write", spectra, new DbForTask(database, false), parameters => parameters.WriteSpectralLibrary = true);

            var psms = SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(folder, "Write", "AllPSMs.psmtsv"), out _);
            Assert.That(psms.Count, Is.EqualTo(2));
            foreach (var psm in psms)
            {
                var peptide = peptides.Single(p => p.FullSequence == psm.FullSequence);
                Assert.That(double.Parse(psm.Notch, CultureInfo.InvariantCulture), Is.EqualTo(1.0029).Within(1e-9)); // the one-notch acceptor's notch
                Assert.That(psm.PrecursorMz, Is.EqualTo(ThirteenCPeakMz(peptide, charge)).Within(1e-4));
            }
            AssertTheoreticalPrecursorMz(Directory.GetFiles(Path.Combine(folder, "Write"), "SpectralLibrary*.msp").Single(), peptides, charge);

            // update a library written at the measured m/z with fewer fragments: PEPTIDEK is replaced, SEQENTIALK is added
            var oldMatch = ThirteenCPickedMatch(peptides[0], charge);
            string oldLibrary = Path.Combine(folder, "measuredMz.msp");
            File.WriteAllText(oldLibrary, new LibrarySpectrum(oldMatch.FullSequence, oldMatch.ScanPrecursorMonoisotopicPeakMz, charge,
                oldMatch.MatchedFragmentIons.Take(3).ToList(), oldMatch.ScanRetentionTime).ToString() + Environment.NewLine);

            RunSearch(folder, "Update", spectra, new DbForTask(database, false), parameters => parameters.UpdateSpectralLibrary = true,
                new DbForTask(oldLibrary, false));

            string updatedLibrary = Directory.GetFiles(Path.Combine(folder, "Update"), "updateSpectralLibrary*.msp").Single();
            var library = new SpectralLibrary(new List<string> { updatedLibrary });
            Assert.That(library.TryGetSpectrum("PEPTIDEK", charge, out var replaced));
            Assert.That(replaced.MatchedFragmentIons.Count, Is.GreaterThan(3));
            library.CloseConnections();
            AssertTheoreticalPrecursorMz(updatedLibrary, peptides, charge);

            // The search task leaves the library it updated open; nothing references it once the task is done.
            GC.Collect();
            GC.WaitForPendingFinalizers();
            Directory.Delete(folder, true);
        }

        private static void RunSearch(string folder, string taskId, string spectra, DbForTask database, Action<SearchParameters> library,
            DbForTask spectralLibrary = null)
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(doPrecursorDeconvolution: false, useProvidedPrecursorInfo: true),
                SearchParameters = new SearchParameters
                {
                    MassDiffAcceptorType = MassDiffAcceptorType.OneMM,
                    DoLabelFreeQuantification = false,
                    WriteMzId = false,
                }
            };
            library(task.SearchParameters);
            var databases = new List<DbForTask> { database };
            if (spectralLibrary != null)
                databases.Add(spectralLibrary);
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { (taskId, task) }, new List<string> { spectra }, databases, folder).Run();
        }

        private static void AssertTheoreticalPrecursorMz(string libraryPath, List<PeptideWithSetModifications> peptides, int charge)
        {
            var library = new SpectralLibrary(new List<string> { libraryPath });
            foreach (var peptide in peptides)
            {
                Assert.That(library.TryGetSpectrum(peptide.FullSequence, charge, out var spectrum), peptide.FullSequence);
                Assert.That(spectrum.PrecursorMz, Is.EqualTo(peptide.MonoisotopicMass.ToMz(charge)).Within(1e-6), peptide.FullSequence);
            }
            Assert.That(library.GetAllLibrarySpectra().Count(), Is.EqualTo(peptides.Count));
            library.CloseConnections();
        }

        /// <summary>
        /// One MS1 and one MS2 scan per peptide. The MS1 holds only the 13C and 13C2 peaks and the MS2 names the 13C peak
        /// as its monoisotopic precursor, as an instrument does when it misses the monoisotopic peak.
        /// </summary>
        private static MsDataFile ThirteenCPickedFile(List<PeptideWithSetModifications> peptides, int charge)
        {
            var scans = new List<MsDataScan>();
            foreach (var peptide in peptides)
            {
                int ms1Number = scans.Count + 1;
                double retentionTime = ms1Number;
                double pickedMz = ThirteenCPeakMz(peptide, charge);

                var ms1 = new MzSpectrum(new[] { pickedMz, (peptide.MonoisotopicMass + 2 * Constants.C13MinusC12).ToMz(charge) }, new[] { 1.0, 0.6 }, false);
                scans.Add(new MsDataScan(ms1, ms1Number, 1, true, Polarity.Positive, retentionTime, new MzRange(100, 2000), "FTMS + p NSI Full ms",
                    MZAnalyzerType.Orbitrap, ms1.SumOfAllY, 1, null, "scan=" + ms1Number));

                var products = new List<Product>();
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
                double[] fragmentMz = products.Select(p => p.NeutralMass.ToMz(1)).OrderBy(mz => mz).ToArray();
                var ms2 = new MzSpectrum(fragmentMz, fragmentMz.Select(_ => 1.0).ToArray(), false);
                scans.Add(new MsDataScan(ms2, ms1Number + 1, 2, true, Polarity.Positive, retentionTime + 0.01, new MzRange(100, 2000), "FTMS + p NSI d Full ms2",
                    MZAnalyzerType.Orbitrap, ms2.SumOfAllY, 1, null, "scan=" + (ms1Number + 1), pickedMz, charge, 1, pickedMz, 2, DissociationType.HCD,
                    ms1Number, pickedMz));
            }
            return new GenericMsDataFile(scans.ToArray(), new SourceFile("no nativeID format", "mzML format", null, null, null));
        }

        /// <summary>
        /// Name, charge and precursor m/z of each spectrum in a MetaMorpheus .msp, in file order.
        /// </summary>
        internal static List<(string Name, int Charge, double PrecursorMz)> ReadMspPrecursors(string mspPath)
        {
            var precursors = new List<(string, int, double)>();
            string name = null;
            foreach (string line in File.ReadLines(mspPath))
            {
                if (line.StartsWith("Name: "))
                    name = line.Substring("Name: ".Length);
                else if (line.StartsWith("MW: "))
                {
                    int slash = name.LastIndexOf('/');
                    precursors.Add((name.Substring(0, slash), int.Parse(name.Substring(slash + 1), CultureInfo.InvariantCulture),
                        double.Parse(line.Substring("MW: ".Length), CultureInfo.InvariantCulture)));
                }
            }
            return precursors;
        }
    }
}
