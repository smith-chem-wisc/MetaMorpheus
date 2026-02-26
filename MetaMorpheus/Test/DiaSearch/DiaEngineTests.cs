// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaSearch;
using MassSpectrometry;
using MassSpectrometry.Dia;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using MzLibUtil;

namespace Test.DiaSearch
{
    [TestFixture]
    public class DiaEngineTests
    {
        #region Helpers

        private static MsDataScan CreateMs2Scan(
            double[] mzs, double[] intensities,
            int oneBasedScanNumber, double rt,
            double isolationMz, double isolationWidth)
        {
            var spectrum = new MzSpectrum(mzs, intensities, shouldCopy: false);
            return new MsDataScan(
                massSpectrum: spectrum,
                oneBasedScanNumber: oneBasedScanNumber,
                msnOrder: 2,
                isCentroid: true,
                polarity: Polarity.Positive,
                retentionTime: rt,
                scanWindowRange: new MzRange(mzs.First(), mzs.Last()),
                scanFilter: null,
                mzAnalyzer: MZAnalyzerType.Orbitrap,
                totalIonCurrent: intensities.Sum(),
                injectionTime: null,
                noiseData: null,
                nativeId: $"scan={oneBasedScanNumber}",
                isolationMZ: isolationMz,
                isolationWidth: isolationWidth,
                dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: 1,
                selectedIonMonoisotopicGuessMz: isolationMz);
        }

        private static LibrarySpectrum CreateLibrarySpectrum(
            string sequence, double precursorMz, int charge,
            double[] fragMzs, double[] fragIntensities,
            double? rt = null, bool isDecoy = false)
        {
            var ions = new List<MatchedFragmentIon>();
            for (int i = 0; i < fragMzs.Length; i++)
            {
                var product = new Product(ProductType.b, FragmentationTerminus.N, fragMzs[i], i + 1, 0, 0);
                ions.Add(new MatchedFragmentIon(product, fragMzs[i], fragIntensities[i], 1));
            }
            return new LibrarySpectrum(sequence, precursorMz, charge, ions, rt, isDecoy);
        }

        private static MsDataScan[] CreateSimpleDiaScans()
        {
            var scans = new MsDataScan[10];
            var rng = new Random(42);
            for (int i = 0; i < 10; i++)
            {
                double rt = 10.0 + i * 0.5;
                var mzs = new double[] { 300.0, 400.0, 500.0, 600.0, 700.0 };
                var intensities = new double[5];
                for (int j = 0; j < 5; j++)
                    intensities[j] = 1000.0 + rng.NextDouble() * 500.0;
                scans[i] = CreateMs2Scan(mzs, intensities, i + 1, rt,
                    isolationMz: 500.0, isolationWidth: 25.0);
            }
            return scans;
        }

        #endregion

        #region ConvertLibrarySpectra Tests

        [Test]
        public void ConvertLibrarySpectra_CorrectCount()
        {
            var spectra = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("PEPTIDE", 500.0, 2, new[] { 300.0, 400.0, 500.0 }, new[] { 1000.0, 800.0, 500.0 }, rt: 12.0),
                CreateLibrarySpectrum("ANOTHER", 510.0, 2, new[] { 350.0, 450.0 }, new[] { 900.0, 700.0 }, rt: 13.0)
            };
            var inputs = DiaEngine.ConvertLibrarySpectra(spectra);
            Assert.That(inputs.Length, Is.EqualTo(2));
        }

        [Test]
        public void ConvertLibrarySpectra_PreservesSequenceAndCharge()
        {
            var spectra = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("PEPTIDE", 500.0, 3, new[] { 300.0 }, new[] { 1000.0 }, rt: 12.0)
            };
            var inputs = DiaEngine.ConvertLibrarySpectra(spectra);
            Assert.That(inputs[0].Sequence, Is.EqualTo("PEPTIDE"));
            Assert.That(inputs[0].ChargeState, Is.EqualTo(3));
        }

        [Test]
        public void ConvertLibrarySpectra_PreservesPrecursorMzAsDouble()
        {
            var spectra = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("TEST", 500.123456789, 2, new[] { 300.0 }, new[] { 1000.0 })
            };
            var inputs = DiaEngine.ConvertLibrarySpectra(spectra);
            Assert.That(inputs[0].PrecursorMz, Is.EqualTo(500.123456789));
        }

        [Test]
        public void ConvertLibrarySpectra_ConvertsFragmentsToFloat()
        {
            var spectra = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("TEST", 500.0, 2, new[] { 300.123, 400.456 }, new[] { 1000.789, 800.012 })
            };
            var inputs = DiaEngine.ConvertLibrarySpectra(spectra);
            Assert.That(inputs[0].FragmentMzs.Length, Is.EqualTo(2));
            Assert.That(inputs[0].FragmentMzs[0], Is.EqualTo(300.123f).Within(0.001f));
            Assert.That(inputs[0].FragmentIntensities[0], Is.EqualTo(1000.789f).Within(0.01f));
        }

        [Test]
        public void ConvertLibrarySpectra_HandlesNullRetentionTime()
        {
            var spectra = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("TEST", 500.0, 2, new[] { 300.0 }, new[] { 1000.0 }, rt: null)
            };
            var inputs = DiaEngine.ConvertLibrarySpectra(spectra);
            Assert.That(inputs[0].RetentionTime.HasValue, Is.False);
        }

        [Test]
        public void ConvertLibrarySpectra_PreservesDecoyFlag()
        {
            var spectra = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("TARGET", 500.0, 2, new[] { 300.0 }, new[] { 1000.0 }, isDecoy: false),
                CreateLibrarySpectrum("DECOY", 510.0, 2, new[] { 350.0 }, new[] { 900.0 }, isDecoy: true)
            };
            var inputs = DiaEngine.ConvertLibrarySpectra(spectra);
            Assert.That(inputs[0].IsDecoy, Is.False);
            Assert.That(inputs[1].IsDecoy, Is.True);
        }

        [Test]
        public void ConvertLibrarySpectra_EmptyList_ReturnsEmptyArray()
        {
            var inputs = DiaEngine.ConvertLibrarySpectra(new List<LibrarySpectrum>());
            Assert.That(inputs.Length, Is.EqualTo(0));
        }

        #endregion

        #region DiaEngine Pipeline Tests

        [Test]
        public void DiaEngine_RunsWithoutException()
        {
            var scans = CreateSimpleDiaScans();
            var library = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("PEPTIDE", 500.0, 2, new[] { 300.0, 400.0, 500.0 }, new[] { 1000.0, 800.0, 500.0 }, rt: 12.0)
            };
            var diaParams = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 5.0f, MinFragmentsRequired = 1, MaxThreads = 1 };
            var engine = new DiaEngine(scans, library, diaParams, new CommonParameters(),
                new List<(string, CommonParameters)>(), new List<string> { "TestTask", "DiaEngine" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results, Is.Not.Null);
            Assert.That(results.DiaResults, Is.Not.Null);
        }

        [Test]
        public void DiaEngine_ProducesResults_WhenFragmentsMatch()
        {
            var scans = CreateSimpleDiaScans();
            var library = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("PEPTIDE", 500.0, 2, new[] { 300.0, 400.0, 500.0 }, new[] { 1000.0, 800.0, 500.0 }, rt: 12.0)
            };
            var diaParams = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 5.0f, MinFragmentsRequired = 1, MaxThreads = 1 };
            var engine = new DiaEngine(scans, library, diaParams, new CommonParameters(),
                new List<(string, CommonParameters)>(), new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.DiaResults.Count, Is.GreaterThan(0),
                "Should detect at least one precursor when fragment m/z values match scan peaks");
        }

        [Test]
        public void DiaEngine_EmptyLibrary_ReturnsEmptyResults()
        {
            var scans = CreateSimpleDiaScans();
            var diaParams = new DiaSearchParameters { MaxThreads = 1 };
            var engine = new DiaEngine(scans, new List<LibrarySpectrum>(), diaParams, new CommonParameters(),
                new List<(string, CommonParameters)>(), new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.DiaResults.Count, Is.EqualTo(0));
        }

        [Test]
        public void DiaEngine_PrecursorOutsideWindow_Skipped()
        {
            var scans = CreateSimpleDiaScans(); // window: 487.5–512.5
            var library = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("OUTSIDE", 800.0, 2, new[] { 300.0, 400.0, 500.0 }, new[] { 1000.0, 800.0, 500.0 }, rt: 12.0)
            };
            var diaParams = new DiaSearchParameters { MaxThreads = 1, MinFragmentsRequired = 1 };
            var engine = new DiaEngine(scans, library, diaParams, new CommonParameters(),
                new List<(string, CommonParameters)>(), new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.DiaResults.Count, Is.EqualTo(0),
                "Precursor outside all DIA windows should produce no results");
        }

        [Test]
        public void DiaEngine_NullScans_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() =>
                new DiaEngine(null, new List<LibrarySpectrum>(), new DiaSearchParameters(),
                    new CommonParameters(), new List<(string, CommonParameters)>(), new List<string>()));
        }

        [Test]
        public void DiaEngine_NullLibrary_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() =>
                new DiaEngine(new MsDataScan[0], null, new DiaSearchParameters(),
                    new CommonParameters(), new List<(string, CommonParameters)>(), new List<string>()));
        }

        #endregion

        #region DiaEngineResults Tests

        [Test]
        public void DiaEngineResults_CountsTargetsAndDecoys()
        {
            var scans = CreateSimpleDiaScans();
            var library = new List<LibrarySpectrum>
            {
                CreateLibrarySpectrum("TARGET", 500.0, 2, new[] { 300.0, 400.0, 500.0 }, new[] { 1000.0, 800.0, 500.0 }, rt: 12.0, isDecoy: false),
                CreateLibrarySpectrum("DECOY", 500.0, 2, new[] { 300.0, 400.0, 500.0 }, new[] { 1000.0, 800.0, 500.0 }, rt: 12.0, isDecoy: true)
            };
            var diaParams = new DiaSearchParameters { MinFragmentsRequired = 1, MaxThreads = 1 };
            var engine = new DiaEngine(scans, library, diaParams, new CommonParameters(),
                new List<(string, CommonParameters)>(), new List<string> { "TestTask" });

            var results = (DiaEngineResults)engine.Run();

            Assert.That(results.TargetLibraryCount, Is.EqualTo(1));
            Assert.That(results.DecoyLibraryCount, Is.EqualTo(1));
        }

        #endregion
    }
}
