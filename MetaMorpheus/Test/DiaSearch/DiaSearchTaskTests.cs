// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.DiaSearch;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.DiaSearch
{
    [TestFixture]
    public class DiaSearchTaskTests
    {
        private string _outputFolder;

        [SetUp]
        public void SetUp()
        {
            _outputFolder = Path.Combine(Path.GetTempPath(), "DiaSearchTaskTests_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_outputFolder);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_outputFolder))
            {
                try { Directory.Delete(_outputFolder, true); }
                catch { /* best effort cleanup */ }
            }
        }

        #region Helper Methods

        /// <summary>
        /// Creates synthetic DIA MS2 scans with known isolation windows and fragment peaks.
        /// Windows: 400-425, 425-450, 450-475, 475-500 m/z (4 windows).
        /// 10 scans per window, RT 1.0–10.0 min.
        /// Each scan contains 50 peaks spread across the window's fragment range.
        /// </summary>
        private static MsDataScan[] CreateSyntheticDiaScans()
        {
            var scans = new List<MsDataScan>();
            double[] windowCenters = { 412.5, 437.5, 462.5, 487.5 };
            double windowWidth = 25.0;
            int scanNumber = 1;

            foreach (double center in windowCenters)
            {
                for (int s = 0; s < 10; s++)
                {
                    double rt = 1.0 + s; // RT from 1.0 to 10.0
                    int peakCount = 50;
                    double[] mzArray = new double[peakCount];
                    double[] intensityArray = new double[peakCount];

                    for (int p = 0; p < peakCount; p++)
                    {
                        mzArray[p] = 200.0 + p * 10.0; // m/z 200–690
                        intensityArray[p] = 1000.0 + p * 100.0;
                    }

                    Array.Sort(mzArray, intensityArray);

                    var scan = new MsDataScan(
                        new MzSpectrum(mzArray, intensityArray, false),
                        scanNumber++, 2, true, Polarity.Positive, rt,
                        new MzRange(100, 2000), "FTMS", MZAnalyzerType.Orbitrap,
                        intensityArray.Sum(), 20.0, null, $"scan={scanNumber - 1}",
                        isolationMZ: center, isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD);

                    scans.Add(scan);
                }
            }

            return scans.ToArray();
        }

        /// <summary>
        /// Creates a synthetic spectral library with precursors that map to the synthetic DIA windows.
        /// </summary>
        private static List<LibrarySpectrum> CreateSyntheticLibrary(int count = 5, bool includeDecoys = true)
        {
            var library = new List<LibrarySpectrum>();
            double[] precursorMzValues = { 410.0, 415.0, 435.0, 460.0, 485.0 };

            for (int i = 0; i < Math.Min(count, precursorMzValues.Length); i++)
            {
                // Create fragment ions with m/z values that exist in our synthetic scans
                var fragments = new List<MatchedFragmentIon>();
                for (int f = 0; f < 6; f++)
                {
                    double fragMz = 200.0 + f * 10.0; // matches synthetic scan peaks
                    var product = new Product(ProductType.b, FragmentationTerminus.N,
                        fragMz - 1.00727647, f + 1, f + 1, 0);
                    fragments.Add(new MatchedFragmentIon(product, fragMz, 1000.0 + f * 200.0, 1));
                }

                library.Add(new LibrarySpectrum(
                    sequence: $"PEPTIDE{(char)('A' + i)}K",
                    precursorMz: precursorMzValues[i],
                    chargeState: 2,
                    peaks: fragments,
                    rt: 5.0, // center of our RT range
                    isDecoy: false));

                if (includeDecoys)
                {
                    library.Add(new LibrarySpectrum(
                        sequence: $"DECOY_PEPTIDE{(char)('A' + i)}K",
                        precursorMz: precursorMzValues[i] + 0.5,
                        chargeState: 2,
                        peaks: fragments,
                        rt: 5.0,
                        isDecoy: true));
                }
            }

            return library;
        }

        #endregion

        #region MetaMorpheusDiaSearchParameters Tests

        [Test]
        public void MetaMorpheusDiaSearchParameters_DefaultValues()
        {
            var p = new MetaMorpheusDiaSearchParameters();

            Assert.That(p.PpmTolerance, Is.EqualTo(20f));
            Assert.That(p.RtToleranceMinutes, Is.EqualTo(5.0f));
            Assert.That(p.MinFragmentsRequired, Is.EqualTo(3));
            Assert.That(p.MinScoreThreshold, Is.EqualTo(0.0f));
            Assert.That(p.PreferGpu, Is.False);
            Assert.That(p.MaxThreads, Is.EqualTo(-1));
            Assert.That(p.SpectralLibraryPath, Is.EqualTo(""));
            Assert.That(p.WriteDecoyResults, Is.False);
            Assert.That(p.WriteDiagnostics, Is.True);
        }

        [Test]
        public void MetaMorpheusDiaSearchParameters_ToMzLibParameters_MapsCorrectly()
        {
            var mmParams = new MetaMorpheusDiaSearchParameters
            {
                PpmTolerance = 15f,
                RtToleranceMinutes = 3.0f,
                MinFragmentsRequired = 5,
                MinScoreThreshold = 0.5f,
                MaxThreads = 8,
                PreferGpu = true
            };

            var mzLibParams = mmParams.ToMzLibParameters();

            Assert.That(mzLibParams.PpmTolerance, Is.EqualTo(15f));
            Assert.That(mzLibParams.RtToleranceMinutes, Is.EqualTo(3.0f));
            Assert.That(mzLibParams.MinFragmentsRequired, Is.EqualTo(5));
            Assert.That(mzLibParams.MinScoreThreshold, Is.EqualTo(0.5f));
            Assert.That(mzLibParams.MaxThreads, Is.EqualTo(8));
            Assert.That(mzLibParams.PreferGpu, Is.True);
        }

        [Test]
        public void MetaMorpheusDiaSearchParameters_TomlRoundTrip()
        {
            var original = new MetaMorpheusDiaSearchParameters
            {
                PpmTolerance = 15f,
                RtToleranceMinutes = 3.5f,
                MinFragmentsRequired = 4,
                MinScoreThreshold = 0.3f,
                PreferGpu = true,
                MaxThreads = 4,
                SpectralLibraryPath = @"C:\data\library.msp",
                WriteDecoyResults = true,
                WriteDiagnostics = false
            };

            // Write to TOML
            string tomlPath = Path.Combine(_outputFolder, "test_params.toml");
            Toml.WriteFile(original, tomlPath, MetaMorpheusTask.tomlConfig);

            // Read back
            var loaded = Toml.ReadFile<MetaMorpheusDiaSearchParameters>(tomlPath, MetaMorpheusTask.tomlConfig);

            Assert.That(loaded.PpmTolerance, Is.EqualTo(original.PpmTolerance));
            Assert.That(loaded.RtToleranceMinutes, Is.EqualTo(original.RtToleranceMinutes));
            Assert.That(loaded.MinFragmentsRequired, Is.EqualTo(original.MinFragmentsRequired));
            Assert.That(loaded.MinScoreThreshold, Is.EqualTo(original.MinScoreThreshold));
            Assert.That(loaded.PreferGpu, Is.EqualTo(original.PreferGpu));
            Assert.That(loaded.MaxThreads, Is.EqualTo(original.MaxThreads));
            Assert.That(loaded.SpectralLibraryPath, Is.EqualTo(original.SpectralLibraryPath));
            Assert.That(loaded.WriteDecoyResults, Is.EqualTo(original.WriteDecoyResults));
            Assert.That(loaded.WriteDiagnostics, Is.EqualTo(original.WriteDiagnostics));
        }

        #endregion

        #region DiaSearchTask TOML Tests

        [Test]
        public void DiaSearchTask_TomlRoundTrip()
        {
            var task = new DiaSearchTask();
            task.DiaSearchParameters.PpmTolerance = 15f;
            task.DiaSearchParameters.RtToleranceMinutes = 3.0f;
            task.DiaSearchParameters.MinFragmentsRequired = 5;
            task.DiaSearchParameters.SpectralLibraryPath = @"C:\data\library.msp";

            string tomlPath = Path.Combine(_outputFolder, "DiaSearchTask.toml");
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);

            // Verify the TOML file contains expected content
            string tomlContent = File.ReadAllText(tomlPath);
            Assert.That(tomlContent, Does.Contain("PpmTolerance"));
            Assert.That(tomlContent, Does.Contain("RtToleranceMinutes"));

            // Read back as DiaSearchTask
            var loaded = Toml.ReadFile<DiaSearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);

            Assert.That(loaded.DiaSearchParameters.PpmTolerance, Is.EqualTo(15f));
            Assert.That(loaded.DiaSearchParameters.RtToleranceMinutes, Is.EqualTo(3.0f));
            Assert.That(loaded.DiaSearchParameters.MinFragmentsRequired, Is.EqualTo(5));
            Assert.That(loaded.DiaSearchParameters.SpectralLibraryPath, Is.EqualTo(@"C:\data\library.msp"));
        }

        #endregion

        #region DiaEngine Direct Tests (Synthetic Data)

        [Test]
        public void DiaEngine_RunsWithSyntheticData_ProducesResults()
        {
            var scans = CreateSyntheticDiaScans();
            var library = CreateSyntheticLibrary(count: 3, includeDecoys: false);
            var mzLibParams = new MassSpectrometry.Dia.DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 1 // low threshold for synthetic data
            };

            using var engine = new DiaEngine(
                scans, library, mzLibParams,
                new CommonParameters(),
                new List<(string, CommonParameters)>(),
                new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results, Is.Not.Null);
            Assert.That(results.DiaResults, Is.Not.Null);
            // With synthetic data where fragment peaks match, we should get results
            Assert.That(results.DiaResults.Count, Is.GreaterThan(0));
        }

        [Test]
        public void DiaEngine_EmptyLibrary_ProducesEmptyResults()
        {
            var scans = CreateSyntheticDiaScans();
            var emptyLibrary = new List<LibrarySpectrum>();
            var mzLibParams = new MassSpectrometry.Dia.DiaSearchParameters();

            using var engine = new DiaEngine(
                scans, emptyLibrary, mzLibParams,
                new CommonParameters(),
                new List<(string, CommonParameters)>(),
                new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.DiaResults, Is.Not.Null);
            Assert.That(results.DiaResults.Count, Is.EqualTo(0));
        }

        [Test]
        public void DiaEngine_PrecursorOutsideAllWindows_NoResults()
        {
            var scans = CreateSyntheticDiaScans(); // windows: 400-425, 425-450, 450-475, 475-500
            var fragments = new List<MatchedFragmentIon>();
            for (int f = 0; f < 6; f++)
            {
                var product = new Product(ProductType.b, FragmentationTerminus.N,
                    200.0 + f * 10.0, f + 1, f + 1, 0);
                fragments.Add(new MatchedFragmentIon(product, 200.0 + f * 10.0, 1000.0, 1));
            }

            var library = new List<LibrarySpectrum>
            {
                new LibrarySpectrum("OUTSIDEPEPTIDE", 900.0, 2, fragments, 5.0) // way outside all windows
            };

            var mzLibParams = new MassSpectrometry.Dia.DiaSearchParameters { MinFragmentsRequired = 1 };

            using var engine = new DiaEngine(
                scans, library, mzLibParams,
                new CommonParameters(),
                new List<(string, CommonParameters)>(),
                new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.DiaResults.Count, Is.EqualTo(0));
        }

        [Test]
        public void DiaEngine_NullScans_Throws()
        {
            Assert.Throws<ArgumentNullException>(() =>
            {
                new DiaEngine(
                    null,
                    new List<LibrarySpectrum>(),
                    new MassSpectrometry.Dia.DiaSearchParameters(),
                    new CommonParameters(),
                    new List<(string, CommonParameters)>(),
                    new List<string>());
            });
        }

        [Test]
        public void DiaEngine_NullLibrary_Throws()
        {
            Assert.Throws<ArgumentNullException>(() =>
            {
                new DiaEngine(
                    CreateSyntheticDiaScans(),
                    null,
                    new MassSpectrometry.Dia.DiaSearchParameters(),
                    new CommonParameters(),
                    new List<(string, CommonParameters)>(),
                    new List<string>());
            });
        }

        [Test]
        public void DiaEngine_NullParams_Throws()
        {
            Assert.Throws<ArgumentNullException>(() =>
            {
                new DiaEngine(
                    CreateSyntheticDiaScans(),
                    new List<LibrarySpectrum>(),
                    null,
                    new CommonParameters(),
                    new List<(string, CommonParameters)>(),
                    new List<string>());
            });
        }

        [Test]
        public void DiaEngine_TargetDecoyCount()
        {
            var scans = CreateSyntheticDiaScans();
            var library = CreateSyntheticLibrary(count: 3, includeDecoys: true);
            var mzLibParams = new MassSpectrometry.Dia.DiaSearchParameters
            {
                MinFragmentsRequired = 1,
                MinScoreThreshold = 0.0f
            };

            using var engine = new DiaEngine(
                scans, library, mzLibParams,
                new CommonParameters(),
                new List<(string, CommonParameters)>(),
                new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.TargetCount, Is.GreaterThanOrEqualTo(0));
            Assert.That(results.DecoyCount, Is.GreaterThanOrEqualTo(0));
            Assert.That(results.TargetCount + results.DecoyCount, Is.EqualTo(results.DiaResults.Count));
            // Verify the aliases used by Phase 8 tests also work
            Assert.That(results.TargetLibraryCount, Is.EqualTo(results.TargetCount));
            Assert.That(results.DecoyLibraryCount, Is.EqualTo(results.DecoyCount));
        }

        #endregion

        #region ConvertLibrarySpectra Tests

        [Test]
        public void ConvertLibrarySpectra_CorrectCount()
        {
            var library = CreateSyntheticLibrary(count: 3, includeDecoys: false);
            var inputs = DiaEngine.ConvertLibrarySpectra(library);
            Assert.That(inputs.Length, Is.EqualTo(library.Count));
        }

        [Test]
        public void ConvertLibrarySpectra_PreservesSequenceAndCharge()
        {
            var library = CreateSyntheticLibrary(count: 1, includeDecoys: false);
            var inputs = DiaEngine.ConvertLibrarySpectra(library);

            Assert.That(inputs[0].Sequence, Is.EqualTo(library[0].Sequence));
            Assert.That(inputs[0].ChargeState, Is.EqualTo(library[0].ChargeState));
        }

        [Test]
        public void ConvertLibrarySpectra_PrecursorMzAsDouble()
        {
            var library = CreateSyntheticLibrary(count: 1, includeDecoys: false);
            var inputs = DiaEngine.ConvertLibrarySpectra(library);

            Assert.That(inputs[0].PrecursorMz, Is.EqualTo(library[0].PrecursorMz));
        }

        [Test]
        public void ConvertLibrarySpectra_FragmentsAsFloat()
        {
            var library = CreateSyntheticLibrary(count: 1, includeDecoys: false);
            var inputs = DiaEngine.ConvertLibrarySpectra(library);

            Assert.That(inputs[0].FragmentMzs.Length, Is.EqualTo(library[0].MatchedFragmentIons.Count));
            Assert.That(inputs[0].FragmentIntensities.Length, Is.EqualTo(library[0].MatchedFragmentIons.Count));

            for (int f = 0; f < inputs[0].FragmentMzs.Length; f++)
            {
                Assert.That(inputs[0].FragmentMzs[f],
                    Is.EqualTo((float)library[0].MatchedFragmentIons[f].Mz).Within(0.001f));
                Assert.That(inputs[0].FragmentIntensities[f],
                    Is.EqualTo((float)library[0].MatchedFragmentIons[f].Intensity).Within(0.1f));
            }
        }

        [Test]
        public void ConvertLibrarySpectra_NullRetentionTime()
        {
            var fragments = new List<MatchedFragmentIon>();
            var product = new Product(ProductType.b, FragmentationTerminus.N, 199.27, 1, 1, 0);
            fragments.Add(new MatchedFragmentIon(product, 200.0, 1000.0, 1));

            var library = new List<LibrarySpectrum>
            {
                new LibrarySpectrum("PEPTIDEK", 400.0, 2, fragments, rt: null)
            };

            var inputs = DiaEngine.ConvertLibrarySpectra(library);

            Assert.That(inputs[0].RetentionTime, Is.Null);
        }

        [Test]
        public void ConvertLibrarySpectra_DecoyFlag()
        {
            var library = CreateSyntheticLibrary(count: 1, includeDecoys: true);
            var inputs = DiaEngine.ConvertLibrarySpectra(library);

            var target = inputs.First(x => !x.IsDecoy);
            var decoy = inputs.First(x => x.IsDecoy);

            Assert.That(target.IsDecoy, Is.False);
            Assert.That(decoy.IsDecoy, Is.True);
        }

        [Test]
        public void ConvertLibrarySpectra_EmptyList()
        {
            var inputs = DiaEngine.ConvertLibrarySpectra(new List<LibrarySpectrum>());
            Assert.That(inputs.Length, Is.EqualTo(0));
        }

        [Test]
        public void ConvertLibrarySpectra_NullList()
        {
            var inputs = DiaEngine.ConvertLibrarySpectra(null);
            Assert.That(inputs.Length, Is.EqualTo(0));
        }

        #endregion

        #region TSV Output Tests

        [Test]
        public void TsvOutput_HasCorrectColumns()
        {
            var scans = CreateSyntheticDiaScans();
            var library = CreateSyntheticLibrary(count: 2, includeDecoys: false);
            var mzLibParams = new MassSpectrometry.Dia.DiaSearchParameters { MinFragmentsRequired = 1 };

            using var engine = new DiaEngine(
                scans, library, mzLibParams,
                new CommonParameters(),
                new List<(string, CommonParameters)>(),
                new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            // Write TSV manually to test the format
            string tsvPath = Path.Combine(_outputFolder, "test_results.tsv");
            using (var writer = new StreamWriter(tsvPath))
            {
                // Write header
                writer.WriteLine(string.Join("\t",
                    "File Name", "Sequence", "Precursor m/z", "Charge", "Window ID",
                    "Is Decoy", "Dot Product Score", "Spectral Angle Score",
                    "Fragments Detected", "Fragments Queried", "Fragment Detection Rate",
                    "Library RT", "RT Window Start", "RT Window End",
                    "XIC Point Counts", "Extracted Intensities"));

                foreach (var r in results.DiaResults)
                {
                    writer.WriteLine(string.Join("\t",
                        "test_file",
                        r.Sequence,
                        r.PrecursorMz.ToString("F4", CultureInfo.InvariantCulture),
                        r.ChargeState,
                        r.WindowId,
                        r.IsDecoy ? "TRUE" : "FALSE",
                        float.IsNaN(r.DotProductScore) ? "NaN" : r.DotProductScore.ToString("F4", CultureInfo.InvariantCulture),
                        float.IsNaN(r.SpectralAngleScore) ? "NaN" : r.SpectralAngleScore.ToString("F4", CultureInfo.InvariantCulture),
                        r.FragmentsDetected,
                        r.FragmentsQueried,
                        r.FragmentDetectionRate.ToString("F4", CultureInfo.InvariantCulture),
                        r.LibraryRetentionTime?.ToString("F2", CultureInfo.InvariantCulture) ?? "",
                        r.RtWindowStart.ToString("F2", CultureInfo.InvariantCulture),
                        r.RtWindowEnd.ToString("F2", CultureInfo.InvariantCulture),
                        string.Join(";", r.XicPointCounts),
                        string.Join(";", r.ExtractedIntensities.Select(x => x.ToString("G6", CultureInfo.InvariantCulture)))));
                }
            }

            // Validate
            string[] lines = File.ReadAllLines(tsvPath);
            Assert.That(lines.Length, Is.GreaterThan(1)); // header + at least one result

            string[] header = lines[0].Split('\t');
            Assert.That(header.Length, Is.EqualTo(16));
            Assert.That(header[0], Is.EqualTo("File Name"));
            Assert.That(header[1], Is.EqualTo("Sequence"));
            Assert.That(header[5], Is.EqualTo("Is Decoy"));
            Assert.That(header[15], Is.EqualTo("Extracted Intensities"));

            // Validate data row column count
            if (lines.Length > 1)
            {
                string[] dataRow = lines[1].Split('\t');
                Assert.That(dataRow.Length, Is.EqualTo(16));
            }
        }

        #endregion
    }
}
