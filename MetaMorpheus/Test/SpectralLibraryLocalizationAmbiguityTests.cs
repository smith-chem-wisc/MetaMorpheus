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
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// A match ambiguous only in where its modification sits (PEPAIDET[Phospho]SEQVENAK or PEPAIDETS[Phospho]EQVENAK)
    /// keeps its mass but has no full sequence. A spectral library entry is named by its full sequence, so such a match
    /// has no entry to write. Before, every such match at one charge shared the group key (null, charge), and the best of
    /// them was written under an empty name ("Name: /2").
    /// </summary>
    [TestFixture]
    public static class SpectralLibraryLocalizationAmbiguityTests
    {
        private const int Charge = 2;
        private const string PhosphoOnT = "PEPAIDET[Common Biological:Phosphorylation on T]SEQVENAK";
        private const string PhosphoOnS = "PEPAIDETS[Common Biological:Phosphorylation on S]EQVENAK";

        [Test]
        public static void SearchTaskWritesAndUpdatesLibrariesWithoutLocalizationAmbiguousMatches()
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, "SpectralLibraryLocalizationAmbiguity");
            if (Directory.Exists(folder))
                Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);

            var ambiguous = new Protein("PEPAIDETSEQVENAK", "P1");
            var unambiguous = new Protein("PEPTIDEK", "P2");
            string database = Path.Combine(folder, "proteins.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ambiguous, unambiguous }, database);

            // The phosphopeptide's scan holds only the fragments its two sites share: b1-b7 (PEPAIDE) and y1-y7 (EQVENAK).
            var phosphoOnT = new PeptideWithSetModifications(PhosphoOnT, GlobalVariables.AllModsKnownDictionary, p: ambiguous);
            var peptideTidek = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Modification>(), p: unambiguous);
            string spectra = Path.Combine(folder, "localizationAmbiguous.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(DataFile(new List<(PeptideWithSetModifications, Func<Product, bool>)>
            {
                (phosphoOnT, p => p.NeutralLoss == 0 && p.FragmentNumber <= 7),
                (peptideTidek, p => true),
            }), spectra, false);

            // write a library
            RunSearch(folder, "Write", spectra, database, parameters => parameters.WriteSpectralLibrary = true);

            // the match is ambiguous only in its site, and the psmtsv still reports it
            var psms = SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(folder, "Write", "AllPSMs.psmtsv"), out _);
            Assert.That(psms.Count, Is.EqualTo(2));
            Assert.That(psms.Single(p => p.BaseSeq == "PEPAIDETSEQVENAK").FullSequence.Split('|'), Is.EquivalentTo(new[] { PhosphoOnT, PhosphoOnS }));

            // update a library that holds PEPTIDEK
            var oldEntry = new LibrarySpectrum("PEPTIDEK", peptideTidek.MonoisotopicMass.ToMz(Charge), Charge,
                MatchedIons(peptideTidek).Take(3).ToList(), 1);
            string oldLibrary = Path.Combine(folder, "old.msp");
            File.WriteAllText(oldLibrary, oldEntry + Environment.NewLine);

            RunSearch(folder, "Update", spectra, database, parameters => parameters.UpdateSpectralLibrary = true, oldLibrary);

            // neither library holds the ambiguous match, under its sites' names or under an empty one
            Assert.Multiple(() =>
            {
                Assert.That(LibraryNames(Directory.GetFiles(Path.Combine(folder, "Write"), "SpectralLibrary*.msp").Single()),
                    Is.EqualTo(new[] { "PEPTIDEK/2" }), "written library");
                Assert.That(LibraryNames(Directory.GetFiles(Path.Combine(folder, "Update"), "updateSpectralLibrary*.msp").Single()),
                    Is.EqualTo(new[] { "PEPTIDEK/2" }), "updated library");
            });

            // The search task leaves the library it updated open; nothing references it once the task is done.
            GC.Collect();
            GC.WaitForPendingFinalizers();
            Directory.Delete(folder, true);
        }

        private static void RunSearch(string folder, string taskId, string spectra, string database, Action<SearchParameters> library,
            string spectralLibrary = null)
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(doPrecursorDeconvolution: false, useProvidedPrecursorInfo: true,
                    listOfModsVariable: new List<(string, string)>
                    {
                        ("Common Biological", "Phosphorylation on S"),
                        ("Common Biological", "Phosphorylation on T"),
                    }),
                SearchParameters = new SearchParameters
                {
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                    DoLabelFreeQuantification = false,
                    WriteMzId = false,
                }
            };
            library(task.SearchParameters);
            var databases = new List<DbForTask> { new(database, false) };
            if (spectralLibrary != null)
                databases.Add(new DbForTask(spectralLibrary, false));
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { (taskId, task) }, new List<string> { spectra }, databases, folder).Run();
        }

        /// <summary>
        /// The name of each spectrum in an .msp, in file order, read straight from its "Name: " lines.
        /// </summary>
        private static List<string> LibraryNames(string mspPath) =>
            File.ReadLines(mspPath).Where(l => l.StartsWith("Name: ")).Select(l => l.Substring("Name: ".Length)).ToList();

        private static List<Product> Fragments(PeptideWithSetModifications peptide)
        {
            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            return products;
        }

        private static List<MatchedFragmentIon> MatchedIons(PeptideWithSetModifications peptide) =>
            Fragments(peptide).Select(p => new MatchedFragmentIon(p, p.NeutralMass.ToMz(1), 1, 1)).ToList();

        /// <summary>
        /// One MS1 and one MS2 scan per peptide, the MS2 holding the peptide's fragments that pass the filter, at charge 1.
        /// </summary>
        private static MsDataFile DataFile(List<(PeptideWithSetModifications Peptide, Func<Product, bool> Keep)> peptides)
        {
            var scans = new List<MsDataScan>();
            foreach (var (peptide, keep) in peptides)
            {
                int ms1Number = scans.Count + 1;
                double retentionTime = ms1Number;
                double precursorMz = peptide.MonoisotopicMass.ToMz(Charge);

                var ms1 = new MzSpectrum(new[] { precursorMz, (peptide.MonoisotopicMass + Constants.C13MinusC12).ToMz(Charge) }, new[] { 1.0, 0.6 }, false);
                scans.Add(new MsDataScan(ms1, ms1Number, 1, true, Polarity.Positive, retentionTime, new MzRange(100, 2000), "FTMS + p NSI Full ms",
                    MZAnalyzerType.Orbitrap, ms1.SumOfAllY, 1, null, "scan=" + ms1Number));

                double[] fragmentMz = Fragments(peptide).Where(keep).Select(p => p.NeutralMass.ToMz(1)).OrderBy(mz => mz).ToArray();
                var ms2 = new MzSpectrum(fragmentMz, fragmentMz.Select(_ => 1.0).ToArray(), false);
                scans.Add(new MsDataScan(ms2, ms1Number + 1, 2, true, Polarity.Positive, retentionTime + 0.01, new MzRange(100, 2000), "FTMS + p NSI d Full ms2",
                    MZAnalyzerType.Orbitrap, ms2.SumOfAllY, 1, null, "scan=" + (ms1Number + 1), precursorMz, Charge, 1, precursorMz, 2, DissociationType.HCD,
                    ms1Number, precursorMz));
            }
            return new GenericMsDataFile(scans.ToArray(), new SourceFile("no nativeID format", "mzML format", null, null, null));
        }
    }
}
