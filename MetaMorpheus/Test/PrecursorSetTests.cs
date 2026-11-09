using Easy.Common.Extensions;
using EngineLayer.Util;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using EngineLayer;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    [TestFixture]
    public class PrecursorSetTests
    {
        private PrecursorSet _precursorSet;
        private Tolerance _tolerance;

        [SetUp]
        public void Setup()
        {
            _tolerance = new PpmTolerance(10);
            _precursorSet = new PrecursorSet(_tolerance);
        }

        [Test]
        public void Constructor_InitializesProperties()
        {
            // Assert
            Assert.That(_precursorSet.DeconPeakTolerance, Is.EqualTo(_tolerance));
            Assert.That(_precursorSet.Count, Is.EqualTo(0));
            Assert.That(_precursorSet.PrecursorDictionary, Is.Not.Null);
            Assert.That(_precursorSet.PrecursorDictionary.Count, Is.EqualTo(0));
        }

        [Test]
        public void Add_NullPrecursor_ReturnsFalse()
        {
            // Act
            bool result = _precursorSet.Add(null);

            // Assert
            Assert.That(result, Is.False);
            Assert.That(_precursorSet.Count, Is.EqualTo(0));
        }

        [Test]
        public void Add_ValidPrecursor_ReturnsTrue()
        {
            // Arrange
            var precursor = new Precursor(500.5, 2, 999.0, 1000.0, 3);

            // Act
            bool result = _precursorSet.Add(precursor);

            // Assert
            Assert.That(result, Is.True);
            Assert.That(_precursorSet.Count, Is.EqualTo(1));
            Assert.That(_precursorSet.PrecursorDictionary.ContainsKey((int)Math.Round(precursor.MonoisotopicPeakMz * 100.0)));
        }

        [Test]
        public void Add_DuplicatePrecursor_AddsAndSanitizesSuccessfully()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5, 2, 999.0, 2000.0, 4); // Same m/z and charge, different intensity

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.True);
            Assert.That(_precursorSet.Count, Is.EqualTo(2));

            // After sanitization, only one should remain
            Assert.That(_precursorSet.IsDirty);
            _precursorSet.Sanitize();
            Assert.That(!_precursorSet.IsDirty);
            Assert.That(_precursorSet.Count, Is.EqualTo(1));
        }

        [Test]
        public void Add_EquivalentPrecursorWithinTolerance_AddsAndSanitizesSuccessfully()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5 + (_tolerance.Value * 500.5 / 1e6 * 0.5), 2, 999.0, 2000.0, 4); // Within tolerance

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.True);
            Assert.That(_precursorSet.Count, Is.EqualTo(2));

            // After sanitization, only one should remain
            Assert.That(_precursorSet.IsDirty);
            _precursorSet.Sanitize();
            Assert.That(!_precursorSet.IsDirty);
            Assert.That(_precursorSet.Count, Is.EqualTo(1));
        }

        [Test]
        public void Add_DifferentChargeSameMz_ReturnsTrue()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5, 3, 1500.0, 2000.0, 4); // Same m/z, different charge

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.True);
            Assert.That(_precursorSet.Count, Is.EqualTo(2));
        }

        [Test]
        public void Add_DifferentMzOutsideTolerance_ReturnsTrue()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5 + (_tolerance.Value * 500.5 / 1e6 * 2), 2, 999.0, 2000.0, 4); // Outside tolerance

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.True);
            Assert.That(_precursorSet.Count, Is.EqualTo(2));
        }

        [Test]
        public void ContainsEquivalent_NullPrecursor_ReturnsFalse()
        {
            // Act
            bool result = _precursorSet.ContainsEquivalent(null, out int integerKey);

            // Assert
            Assert.That(result, Is.False);
        }

        [Test]
        public void ContainsEquivalent_EquivalentPrecursor_ReturnsTrue()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5 + (_tolerance.Value * 500.5 / 1e6 * 0.5), 2, 999.0, 2000.0, 4); // Within tolerance
            _precursorSet.Add(precursor1);

            // Act
            bool result = _precursorSet.ContainsEquivalent(precursor2, out int integerKey);

            // Assert
            Assert.That(result, Is.True);
            Assert.That(integerKey, Is.EqualTo((int)Math.Round(precursor2.MonoisotopicPeakMz * 100.0)));
        }

        [Test]
        public void ContainsEquivalent_DifferentMzOutsideTolerance_ReturnsFalse()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5 + (_tolerance.Value * 500.5 / 1e6 * 2), 2, 999.0, 2000.0, 4); // Outside tolerance
            _precursorSet.Add(precursor1);

            // Act
            bool result = _precursorSet.ContainsEquivalent(precursor2, out int integerKey);

            // Assert
            Assert.That(result, Is.False);
            Assert.That(integerKey, Is.EqualTo((int)Math.Round(precursor2.MonoisotopicPeakMz * 100.0)));
        }

        [Test]
        public void ContainsEquivalent_DifferentCharge_ReturnsFalse()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5, 3, 1500.0, 2000.0, 4); // Same m/z, different charge
            _precursorSet.Add(precursor1);

            // Act
            bool result = _precursorSet.ContainsEquivalent(precursor2, out int integerKey);

            // Assert
            Assert.That(result, Is.False);
            Assert.That(integerKey, Is.EqualTo((int)Math.Round(precursor2.MonoisotopicPeakMz * 100.0)));
        }

        [Test]
        public void GetEnumerator_ReturnsAllPrecursors()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(600.5, 3, 1500.0, 2000.0, 4);
            var precursor3 = new Precursor(700.5, 4, 2000.0, 3000.0, 5);
            _precursorSet.Add(precursor1);
            _precursorSet.Add(precursor2);
            _precursorSet.Add(precursor3);

            // Act
            var allPrecursors = _precursorSet.ToList();

            // Assert
            Assert.That(allPrecursors.Count, Is.EqualTo(3));
            Assert.That(allPrecursors, Contains.Item(precursor1));
            Assert.That(allPrecursors, Contains.Item(precursor2));
            Assert.That(allPrecursors, Contains.Item(precursor3));
        }

        [Test]
        public void Add_MultiplePrecursorsInSameBucket_StoresCorrectly()
        {
            // Arrange
            // These precursors will have the same integer key but are outside tolerance of each other
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.505, 2, 999.0, 2000.0, 4);
            _tolerance = new PpmTolerance(1); // Very tight tolerance
            _precursorSet = new PrecursorSet(_tolerance);

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.True);
            Assert.That(_precursorSet.Count, Is.EqualTo(2));
            
            // Both precursors should be in the same bucket
            int intKey = (int)Math.Round(precursor1.MonoisotopicPeakMz * 100.0);
            Assert.That(_precursorSet.PrecursorDictionary.ContainsKey(intKey));
            Assert.That(_precursorSet.PrecursorDictionary[intKey].Count, Is.EqualTo(2));
        }

        [Test]
        public void ContainsEquivalent_ChecksAdjacentBuckets()
        {
            // Arrange
            // Create a precursor at a bucket boundary
            double mz = 500.995; // This will round to 50100
            var precursor1 = new Precursor(mz, 2, 999.0, 1000.0, 3);
            _precursorSet.Add(precursor1);

            // Create a precursor in the previous bucket but within tolerance
            double deltaMz = (_tolerance.Value * mz / 1e6 * 0.5); // Within tolerance
            var precursor2 = new Precursor(mz - deltaMz, 2, 999.0, 2000.0, 4);

            // Act
            bool result = _precursorSet.ContainsEquivalent(precursor2, out int integerKey);

            // Assert
            Assert.That(result, Is.True);
            // Verify that the integer keys are different
            Assert.That(integerKey, Is.EqualTo((int)Math.Round(precursor2.MonoisotopicPeakMz * 100.0)));
            Assert.That((int)Math.Round(precursor1.MonoisotopicPeakMz * 100.0), Is.Not.EqualTo(integerKey));
        }

        [Test]
        public void Count_ReflectsAddedItems()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(600.5, 3, 1500.0, 2000.0, 4);
            var precursor3 = new Precursor(700.5, 4, 2000.0, 3000.0, 5);
            
            // Act & Assert
            Assert.That(_precursorSet.Count, Is.EqualTo(0));
            
            _precursorSet.Add(precursor1);
            Assert.That(_precursorSet.Count, Is.EqualTo(1));
            
            _precursorSet.Add(precursor2);
            Assert.That(_precursorSet.Count, Is.EqualTo(2));
            
            _precursorSet.Add(precursor3);
            Assert.That(_precursorSet.Count, Is.EqualTo(3));
            
            // Try adding a duplicate
            _precursorSet.Add(precursor1);
            Assert.That(_precursorSet.Count, Is.EqualTo(4), "Count should change when adding duplicate");

            _precursorSet.Sanitize();
            Assert.That(_precursorSet.Count, Is.EqualTo(3), "Count should reflect unique precursors after sanitization");
        }

        [Test]
        public void Enumerate_SanitizesIfNeeded()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(600.5, 3, 1500.0, 2000.0, 4);
            var precursor3 = new Precursor(700.5, 4, 2000.0, 3000.0, 5);

            // Act & Assert
            _precursorSet.Add(precursor1);
            _precursorSet.Add(precursor1);
            _precursorSet.Add(precursor2);
            _precursorSet.Add(precursor3);

            Assert.That(_precursorSet.Count, Is.EqualTo(4));
            Assert.That(_precursorSet.IsDirty, Is.True);

            int count = 0;
            foreach (var precursor in _precursorSet)
            {
                count++;
            }
            Assert.That(count, Is.EqualTo(3), "Enumeration should reflect unique precursors after sanitization");
            Assert.That(_precursorSet.IsDirty, Is.False, "Set should not be dirty after enumeration");
        }

        [Test]
        public void Precursor_AllConstructorsProduceIdenticalPrecursors()
        {
            // Test values
            double mz = 500.0;
            int charge = 2;
            double mass = mz.ToMass(charge);
            double intensity = 1000.0;
            int envelopePeakCount = 4;
            double? fractionalIntensity = 0.5;

            // Create an IsotopicEnvelope for the envelope-based constructor
            var envelope = new IsotopicEnvelope(
                new List<(double mz, double intensity)>
                {
                    (mz, intensity / 2.0),
                    (mz + 1, intensity / 4.0),
                    (mz + 2, intensity / 8.0),
                    (mz + 3, intensity / 8.0)
                },
                mass, charge, intensity, envelopePeakCount);

            // Construct using all available constructors
            var p1 = new Precursor(mz, charge, mass, intensity, envelopePeakCount, fractionalIntensity);
            var p2 = new Precursor(mz, charge, intensity, envelopePeakCount, fractionalIntensity);
            var p3 = new Precursor(charge, mass, intensity, envelopePeakCount, fractionalIntensity);
            var p4 = new Precursor(envelope, intensity, fractionalIntensity);
            var p5 = new Precursor(envelope, null, fractionalIntensity);

            // Compare all pairs for equality
            Assert.That(p1.Equals(p2, new PpmTolerance(1)), Is.True, "p1 != p2");
            Assert.That(p1.Equals(p3, new PpmTolerance(1)), Is.True, "p1 != p3");
            Assert.That(p1.Equals(p4, new PpmTolerance(1)), Is.True, "p1 != p4");
            Assert.That(p1.Equals(p5, new PpmTolerance(1)), Is.True, "p1 != p5");
            Assert.That(p2.Equals(p3, new PpmTolerance(1)), Is.True, "p2 != p3");
            Assert.That(p2.Equals(p4, new PpmTolerance(1)), Is.True, "p2 != p4");
            Assert.That(p2.Equals(p5, new PpmTolerance(1)), Is.True, "p2 != p5");
            Assert.That(p3.Equals(p4, new PpmTolerance(1)), Is.True, "p3 != p4");
            Assert.That(p3.Equals(p5, new PpmTolerance(1)), Is.True, "p3 != p5");
            Assert.That(p4.Equals(p5, new PpmTolerance(1)), Is.True, "p4 != p5");

            // Compare all pairs for hash code equivalence
            Assert.That(p1.GetHashCode(), Is.EqualTo(p2.GetHashCode()), "p1 hash != p2 hash");
            Assert.That(p1.GetHashCode(), Is.EqualTo(p3.GetHashCode()), "p1 hash != p3 hash");
            Assert.That(p1.GetHashCode(), Is.EqualTo(p4.GetHashCode()), "p1 hash != p4 hash");
            Assert.That(p1.GetHashCode(), Is.EqualTo(p5.GetHashCode()), "p1 hash != p5 hash");
            Assert.That(p2.GetHashCode(), Is.EqualTo(p3.GetHashCode()), "p2 hash != p3 hash");
            Assert.That(p2.GetHashCode(), Is.EqualTo(p4.GetHashCode()), "p2 hash != p4 hash");
            Assert.That(p2.GetHashCode(), Is.EqualTo(p5.GetHashCode()), "p2 hash != p5 hash");
            Assert.That(p3.GetHashCode(), Is.EqualTo(p4.GetHashCode()), "p3 hash != p4 hash");
            Assert.That(p3.GetHashCode(), Is.EqualTo(p5.GetHashCode()), "p3 hash != p5 hash");
            Assert.That(p4.GetHashCode(), Is.EqualTo(p5.GetHashCode()), "p4 hash != p5 hash");

            // Check all properties
            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p2.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p2.Charge));
            Assert.That(p1.MonoisotopicMass, Is.EqualTo(p2.MonoisotopicMass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p2.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p2.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p2.FractionalIntensity));

            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p3.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p3.Charge));
            Assert.That(p1.MonoisotopicMass, Is.EqualTo(p3.MonoisotopicMass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p3.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p3.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p3.FractionalIntensity));

            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p4.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p4.Charge));
            Assert.That(p1.MonoisotopicMass, Is.EqualTo(p4.MonoisotopicMass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p4.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p4.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p4.FractionalIntensity));

            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p5.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p5.Charge));
            Assert.That(p1.MonoisotopicMass, Is.EqualTo(p5.MonoisotopicMass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p5.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p5.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p5.FractionalIntensity));
        }

        #region Precursor Merging and Filtering

        private static IsotopicEnvelope GetEnvelope(string sequence, int charge, int intensityMultiplier = 100, double percentIntensityCutoff = 0.01)
        {
            var withSetMods = new PeptideWithSetModifications(sequence, GlobalVariables.AllRnaModsKnownDictionary);
            var distribution = IsotopicDistribution.GetDistribution(withSetMods.ThisChemicalFormula);

            // Create a list of (mz, intensity) where intensity > 0.01
            var peaks = new List<(double mz, double intensity)>();
            for (int i = 0; i < distribution.Masses.Length; i++)
            {
                if (distribution.Intensities[i] > percentIntensityCutoff)
                {
                    peaks.Add((distribution.Masses[i].ToMz(charge), distribution.Intensities[i] * intensityMultiplier));
                }
            }

            return new IsotopicEnvelope(peaks, withSetMods.MonoisotopicMass, charge, peaks.Sum(p => p.intensity), 0.5);
        }

        // Creates a low harmonic envelope of the new charge. 
        private static IsotopicEnvelope GetHarmonicEnvelope(IsotopicEnvelope envelope, int newCharge)
        {
            if (envelope.Charge % newCharge != 0 || newCharge >= envelope.Charge)
                throw new ArgumentException("New charge must be a lower harmonic of the original charge.");

            var chargeRatio = envelope.Charge / newCharge;
            var newPeaks = envelope.Peaks
                .Where((_, index) => index % chargeRatio == 0) // Take every nth peak based on charge ratio
                .ToList();
            var newMass = envelope.MonoisotopicMass.ToMz(envelope.Charge).ToMass(newCharge);

            return new IsotopicEnvelope(newPeaks, newMass, newCharge, newPeaks.Sum(p => p.intensity), 0.5);
        }

        [Test]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        public static void AreSequentialAtSpacing_ReturnTrue(string sequence, int charge)
        {
            var envelope = GetEnvelope(sequence, charge);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).ToList();
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var precursor1 = new Precursor(envelope1);
            var precursor2 = new Precursor(envelope2);
            var expectedSpacing = Constants.C13MinusC12 / charge;

            Assert.That(PrecursorSet.AreSequentialAtSpacing(precursor1, precursor2, expectedSpacing, new PpmTolerance(5) ), Is.True);
        }

        [Test]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        public static void AreSequentialAtSpacing_ReturnFalse(string sequence, int charge)
        {
            var envelope = GetEnvelope(sequence, charge);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).Skip(2).ToList(); // Skip peaks to ensure not sequential
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var precursor1 = new Precursor(envelope1);
            var precursor2 = new Precursor(envelope2);
            var expectedSpacing = Constants.C13MinusC12 / charge;

            Assert.That(PrecursorSet.AreSequentialAtSpacing(precursor1, precursor2, expectedSpacing, new PpmTolerance(5)), Is.False);
        }

        [Test]
        public static void AreSequentialAtSpacing_ReturnFalseOnFailureConditions()
        {
            string sequence = "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE";
            int charge = 2;
            var envelope = GetEnvelope(sequence, charge);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).ToList(); 
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var precursor1 = new Precursor(envelope1);
            var precursor2 = new Precursor(envelope2);
            var expectedSpacing = Constants.C13MinusC12 / charge;

            // Empty peaks in envelopes
            envelope1.Peaks.Clear();
            precursor1 = new Precursor(envelope1);
            Assert.That(PrecursorSet.AreSequentialAtSpacing(precursor1, precursor2, expectedSpacing, new PpmTolerance(5)), Is.False);
            envelope1.Peaks.AddRange(envelope.Peaks.Where(p => p.mz <= mostIntense));

            envelope2.Peaks.Clear();
            precursor2 = new Precursor(envelope2);
            Assert.That(PrecursorSet.AreSequentialAtSpacing(precursor1, precursor2, expectedSpacing, new PpmTolerance(5)), Is.False);
            envelope2.Peaks.AddRange(envelope.Peaks.Where(p => p.mz > mostIntense));

            // null envelope
            var nullEnvelopePrecursor = new Precursor(1, 2, 100, 2, 0, 0);
            Assert.That(PrecursorSet.AreSequentialAtSpacing(nullEnvelopePrecursor, precursor2, expectedSpacing, new PpmTolerance(5)), Is.False);
            Assert.That(PrecursorSet.AreSequentialAtSpacing(precursor1, nullEnvelopePrecursor, expectedSpacing, new PpmTolerance(5)), Is.False);

            // different charges
            var diffChargeEnvelope = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge + 1, peaks2.Sum(p => p.intensity), 0.1);
            var precursorDiffCharge = new Precursor(diffChargeEnvelope);
            Assert.That(PrecursorSet.AreSequentialAtSpacing(precursor1, precursorDiffCharge, expectedSpacing, new PpmTolerance(5)), Is.False);
        }

        [Test]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        public static void MergePrecursors_BaseCase(string sequence, int charge)
        {
            var tol = new PpmTolerance(5);
            var envelope = GetEnvelope(sequence, charge, 100, 0.05);
            var originalPrecursor = new Precursor(envelope);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).ToList();
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var leftPrecursor = new Precursor(envelope1);
            var rightPrecursor = new Precursor(envelope2);

            var mergedPrecursor = PrecursorSet.MergePrecursors(leftPrecursor, rightPrecursor, charge, tol);

            Assert.That(mergedPrecursor.MonoisotopicPeakMz, Is.EqualTo(originalPrecursor.MonoisotopicPeakMz).Within(1e-6));
            Assert.That(mergedPrecursor.Charge, Is.EqualTo(originalPrecursor.Charge));
            Assert.That(mergedPrecursor.MonoisotopicMass, Is.EqualTo(originalPrecursor.MonoisotopicMass).Within(1e-6));
            Assert.That(mergedPrecursor.Intensity, Is.EqualTo(originalPrecursor.Intensity).Within(1e-6));
            Assert.That(mergedPrecursor.EnvelopePeakCount, Is.EqualTo(originalPrecursor.EnvelopePeakCount));
            Assert.That(mergedPrecursor.FractionalIntensity, Is.EqualTo(originalPrecursor.FractionalIntensity));
        }

        [Test]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        public static void MergeSplitEnvelopes_MergesCorrectly(string sequence, int charge)
        {
            var envelope = GetEnvelope(sequence, charge, 100, 0.05);
            var originalPrecursor = new Precursor(envelope);
            var precursorToNotMerge = new Precursor(1000.0, 2, 100, 4);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).ToList();
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var precursor1 = new Precursor(envelope1);
            var precursor2 = new Precursor(envelope2);

            List<Precursor> toMerge = new List<Precursor> { precursor1, precursor2, precursorToNotMerge };
            PrecursorSet.MergeSplitEnvelopes(in toMerge, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(toMerge.Count, Is.EqualTo(2), "Should have merged two precursors into one");
            var mergedPrecursor = toMerge.First(p => p != precursorToNotMerge);
            Assert.That(mergedPrecursor.MonoisotopicPeakMz, Is.EqualTo(originalPrecursor.MonoisotopicPeakMz).Within(1e-6));
            Assert.That(mergedPrecursor.Charge, Is.EqualTo(originalPrecursor.Charge));
            Assert.That(mergedPrecursor.MonoisotopicMass, Is.EqualTo(originalPrecursor.MonoisotopicMass).Within(1e-6));
            Assert.That(mergedPrecursor.Intensity, Is.EqualTo(originalPrecursor.Intensity).Within(1e-6));
            Assert.That(mergedPrecursor.EnvelopePeakCount, Is.EqualTo(originalPrecursor.EnvelopePeakCount));
            Assert.That(mergedPrecursor.FractionalIntensity, Is.EqualTo(originalPrecursor.FractionalIntensity));

            // ensure it works bidirectionally 
            toMerge = new List<Precursor> { precursorToNotMerge, precursor2, precursor1,  };
            PrecursorSet.MergeSplitEnvelopes(in toMerge, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(toMerge.Count, Is.EqualTo(2), "Should have merged two precursors into one");
            mergedPrecursor = toMerge.First(p => p != precursorToNotMerge);
            Assert.That(mergedPrecursor.MonoisotopicPeakMz, Is.EqualTo(originalPrecursor.MonoisotopicPeakMz).Within(1e-6));
            Assert.That(mergedPrecursor.Charge, Is.EqualTo(originalPrecursor.Charge));
            Assert.That(mergedPrecursor.MonoisotopicMass, Is.EqualTo(originalPrecursor.MonoisotopicMass).Within(1e-6));
            Assert.That(mergedPrecursor.Intensity, Is.EqualTo(originalPrecursor.Intensity).Within(1e-6));
            Assert.That(mergedPrecursor.EnvelopePeakCount, Is.EqualTo(originalPrecursor.EnvelopePeakCount));
            Assert.That(mergedPrecursor.FractionalIntensity, Is.EqualTo(originalPrecursor.FractionalIntensity));
        }

        [Test]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        public static void Sanitize_MergeSplitEnvelopes_MergesCorrectly(string sequence, int charge)
        {
            var envelope = GetEnvelope(sequence, charge, 100, 0.05);
            var originalPrecursor = new Precursor(envelope);
            var precursorToNotMerge = new Precursor(1000.0, 2, 100, 4);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).ToList();
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var precursor1 = new Precursor(envelope1);
            var precursor2 = new Precursor(envelope2);

            PrecursorSet set = new PrecursorSet(new PpmTolerance(5), applyEnvelopeFiltering: true)
            {
                precursor1,
                precursor2,
                precursorToNotMerge
            };
            set.Sanitize();

            var toMerge = set.ToList();

            Assert.That(toMerge.Count, Is.EqualTo(2), "Should have merged two precursors into one");
            var mergedPrecursor = toMerge.First(p => p != precursorToNotMerge);
            Assert.That(mergedPrecursor.MonoisotopicPeakMz, Is.EqualTo(originalPrecursor.MonoisotopicPeakMz).Within(1e-6));
            Assert.That(mergedPrecursor.Charge, Is.EqualTo(originalPrecursor.Charge));
            Assert.That(mergedPrecursor.MonoisotopicMass, Is.EqualTo(originalPrecursor.MonoisotopicMass).Within(1e-6));
            Assert.That(mergedPrecursor.Intensity, Is.EqualTo(originalPrecursor.Intensity).Within(1e-6));
            Assert.That(mergedPrecursor.EnvelopePeakCount, Is.EqualTo(originalPrecursor.EnvelopePeakCount));
            Assert.That(mergedPrecursor.FractionalIntensity, Is.EqualTo(originalPrecursor.FractionalIntensity));

            // ensure it works bidirectionally 
            set.Clear();
            set.Add(precursorToNotMerge);
            set.Add(precursor2);
            set.Add(precursor1);
            set.Sanitize();
            toMerge = set.ToList();

            Assert.That(toMerge.Count, Is.EqualTo(2), "Should have merged two precursors into one");
            mergedPrecursor = toMerge.First(p => p != precursorToNotMerge);
            Assert.That(mergedPrecursor.MonoisotopicPeakMz, Is.EqualTo(originalPrecursor.MonoisotopicPeakMz).Within(1e-6));
            Assert.That(mergedPrecursor.Charge, Is.EqualTo(originalPrecursor.Charge));
            Assert.That(mergedPrecursor.MonoisotopicMass, Is.EqualTo(originalPrecursor.MonoisotopicMass).Within(1e-6));
            Assert.That(mergedPrecursor.Intensity, Is.EqualTo(originalPrecursor.Intensity).Within(1e-6));
            Assert.That(mergedPrecursor.EnvelopePeakCount, Is.EqualTo(originalPrecursor.EnvelopePeakCount));
            Assert.That(mergedPrecursor.FractionalIntensity, Is.EqualTo(originalPrecursor.FractionalIntensity));
        }

        [Test]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 2)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5)]
        [TestCase("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 10)]
        public static void MergeSplitEnvelopes_DoesNotMergesOnNonSequentialIsotopicPeaks(string sequence, int charge)
        {
            var envelope = GetEnvelope(sequence, charge, 100, 0.04);

            // split into two envelopes with at most abundant peak
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).Skip(1).ToList(); // Skip peaks to ensure not sequential
            var envelope1 = new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1);

            var precursor1 = new Precursor(envelope1);
            var precursor2 = new Precursor(envelope2);

            List<Precursor> toMerge = new List<Precursor> { precursor1, precursor2 };
            PrecursorSet.MergeSplitEnvelopes(in toMerge, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(toMerge.Count, Is.EqualTo(2), "Should not have merged two precursors into one");
        }

        [Test]
        [TestCase(4, 2)]
        [TestCase(15, 3)]
        [TestCase(15, 5)]
        [TestCase(12, 3)]
        [TestCase(12, 4)]
        [TestCase(12, 6)]
        public static void FilterLowHarmonics(int baseCharge, int harmonicCharge)
        {
            var tol = new PpmTolerance(5);
            string sequence = "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE";
            var baseEnvelope = GetEnvelope(sequence, baseCharge, 100, 0.05);
            var basePrecursor = new Precursor(baseEnvelope);

            // should never get filtered out
            var diffChargeDiffMass = new Precursor(GetEnvelope("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTID", baseCharge / 2, 100, 0.05));
            var sameChargeDiffMass = new Precursor(GetEnvelope("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTID", baseCharge, 100, 0.05));

            var lowHarmonicEnvelope = GetHarmonicEnvelope(baseEnvelope, harmonicCharge);
            var lowHarmonicPrecursor = new Precursor(lowHarmonicEnvelope);

            List<Precursor> precursors = new List<Precursor> { basePrecursor, lowHarmonicPrecursor, diffChargeDiffMass, sameChargeDiffMass };
            PrecursorSet.RemoveLowHarmonics(in precursors, tol);

            Assert.That(precursors.Contains(basePrecursor), Is.True, "Base precursor should not be filtered out");
            Assert.That(precursors.Contains(lowHarmonicPrecursor), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(precursors.Contains(diffChargeDiffMass), Is.True, "Different charge and mass precursor should not be filtered out");
            Assert.That(precursors.Contains(sameChargeDiffMass), Is.True, "Same charge different mass precursor should not be filtered out");

            // Assure it works in revers
            precursors = new List<Precursor> { lowHarmonicPrecursor, basePrecursor, sameChargeDiffMass, diffChargeDiffMass };
            PrecursorSet.RemoveLowHarmonics(in precursors, tol);
            Assert.That(precursors.Contains(basePrecursor), Is.True, "Base precursor should not be filtered out");
            Assert.That(precursors.Contains(lowHarmonicPrecursor), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(precursors.Contains(diffChargeDiffMass), Is.True, "Different charge and mass precursor should not be filtered out");
            Assert.That(precursors.Contains(sameChargeDiffMass), Is.True, "Same charge different mass precursor should not be filtered out");
        }

        [Test]
        public static void Sanitize_CleansHarmonicsAndSplit()
        {
            var tol = new PpmTolerance(5);
            string sequence = "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE";
            var baseEnvelope = GetEnvelope(sequence, 8, 100, 0.05);
            var basePrecursor = new Precursor(baseEnvelope);

            // split into two envelopes with at most abundant peak
            var mostIntense = baseEnvelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = baseEnvelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = baseEnvelope.Peaks.Where(p => p.mz > mostIntense).ToList();
            var envelope1 = new IsotopicEnvelope(peaks1, baseEnvelope.MonoisotopicMass, baseEnvelope.Charge, peaks1.Sum(p => p.intensity), 0.1);
            var envelope2 = new IsotopicEnvelope(peaks2, baseEnvelope.MonoisotopicMass, baseEnvelope.Charge, peaks2.Sum(p => p.intensity), 0.1);
            var splitLeft = new Precursor(envelope1);
            var splitRight = new Precursor(envelope2);

            var lowHarmonicEnvelope = GetHarmonicEnvelope(baseEnvelope, 4);
            var lowHarmonicPrecursor = new Precursor(lowHarmonicEnvelope);
            var lowHarmonicEnvelope2 = GetHarmonicEnvelope(baseEnvelope, 2);
            var lowHarmonicPrecursor2 = new Precursor(lowHarmonicEnvelope2);

            var diffChargeDiffMass = new Precursor(GetEnvelope("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTID", 4, 100, 0.05));
            var sameChargeDiffMass = new Precursor(GetEnvelope("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTID", 8, 100, 0.05));

            PrecursorSet set = new PrecursorSet(tol, applyEnvelopeFiltering: true)
            {
                splitLeft,
                splitRight,
                lowHarmonicPrecursor,
                lowHarmonicPrecursor2,
                diffChargeDiffMass,
                sameChargeDiffMass
            };

            set.Sanitize();

            var survivors = set.ToList();
            Assert.That(survivors.Count, Is.EqualTo(3), "Should have merged two precursors into one and removed harmonics");
            Assert.That(survivors.Contains(basePrecursor), Is.True, "Base precursor should not be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor2), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(diffChargeDiffMass), Is.True, "Different charge and mass precursor should not be filtered out");
            Assert.That(survivors.Contains(sameChargeDiffMass), Is.True, "Same charge different mass precursor should not be filtered out");

            // Reverse the order
            set.Clear();
            set.Add(sameChargeDiffMass);
            set.Add(diffChargeDiffMass);
            set.Add(lowHarmonicPrecursor2);
            set.Add(lowHarmonicPrecursor);
            set.Add(splitRight);
            set.Add(splitLeft);
            set.Sanitize();

            survivors = set.ToList();
            Assert.That(survivors.Count, Is.EqualTo(3), "Should have merged two precursors into one and removed harmonics");
            Assert.That(survivors.Contains(basePrecursor), Is.True, "Base precursor should not be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor2), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(diffChargeDiffMass), Is.True, "Different charge and mass precursor should not be filtered out");
            Assert.That(survivors.Contains(sameChargeDiffMass), Is.True, "Same charge different mass precursor should not be filtered out");


            // Repeat with base added in to test successful merge AND duplicate removal 
            set.Clear();
            set.Add(basePrecursor);
            set.Add(sameChargeDiffMass);
            set.Add(diffChargeDiffMass);
            set.Add(lowHarmonicPrecursor2);
            set.Add(lowHarmonicPrecursor);
            set.Add(splitRight);
            set.Add(splitLeft);
            set.Sanitize();

            survivors = set.ToList();
            Assert.That(survivors.Count, Is.EqualTo(3), "Should have merged two precursors into one and removed harmonics and duplicate");
            Assert.That(survivors.Contains(basePrecursor), Is.True, "Base precursor should not be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor2), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(diffChargeDiffMass), Is.True, "Different charge and mass precursor should not be filtered out");
            Assert.That(survivors.Contains(sameChargeDiffMass), Is.True, "Same charge different mass precursor should not be filtered out");

            // reverse with base added in
            set.Clear();
            set.Add(splitLeft);
            set.Add(splitRight);
            set.Add(lowHarmonicPrecursor);
            set.Add(lowHarmonicPrecursor2);
            set.Add(diffChargeDiffMass);
            set.Add(sameChargeDiffMass);
            set.Add(basePrecursor);
            set.Sanitize();

            survivors = set.ToList();
            Assert.That(survivors.Count, Is.EqualTo(3), "Should have merged two precursors into one and removed harmonics and duplicate");
            Assert.That(survivors.Contains(basePrecursor), Is.True, "Base precursor should not be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(lowHarmonicPrecursor2), Is.False, "Low harmonic precursor should be filtered out");
            Assert.That(survivors.Contains(diffChargeDiffMass), Is.True, "Different charge and mass precursor should not be filtered out");
            Assert.That(survivors.Contains(sameChargeDiffMass), Is.True, "Same charge different mass precursor should not be filtered out");
        }

        // Builds an envelope of evenly spaced peaks, listed in m/z order, with an explicit monoisotopic m/z
        private static Precursor GetLadderPrecursor(double monoMz, double firstPeakMz, int peakCount, int charge, double intensity, double score = 0.5)
        {
            double spacing = Constants.C13MinusC12 / charge;
            var peaks = Enumerable.Range(0, peakCount)
                .Select(k => (firstPeakMz + k * spacing, intensity))
                .ToList();
            var envelope = new IsotopicEnvelope(peaks, monoMz.ToMass(charge), charge, peaks.Sum(p => p.Item2), score);
            return new Precursor(envelope) { DeconvolutionScore = score };
        }

        [Test]
        public static void Sanitize_DoesNotFilterByDefault()
        {
            var envelope = GetEnvelope("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", 5, 100, 0.05);
            var mostIntense = envelope.Peaks.MaxBy(p => p.intensity).mz;
            var peaks1 = envelope.Peaks.Where(p => p.mz <= mostIntense).ToList();
            var peaks2 = envelope.Peaks.Where(p => p.mz > mostIntense).ToList();
            var splitLeft = new Precursor(new IsotopicEnvelope(peaks1, envelope.MonoisotopicMass, envelope.Charge, peaks1.Sum(p => p.intensity), 0.1));
            var splitRight = new Precursor(new IsotopicEnvelope(peaks2, envelope.MonoisotopicMass + 1.0, envelope.Charge, peaks2.Sum(p => p.intensity), 0.1));

            var set = new PrecursorSet(new PpmTolerance(5)) { splitLeft, splitRight };
            set.Sanitize();

            Assert.That(set.ApplyEnvelopeFiltering, Is.False);
            Assert.That(set.Count, Is.EqualTo(2), "Split envelopes should only be merged when filtering is switched on");
            Assert.That(new CommonParameters().ApplyPrecursorEnvelopeFiltering, Is.False);
        }

        [Test]
        public static void MergeSplitEnvelopes_DoesNotMergeAbuttingSpecies()
        {
            // Two peptides whose ladders abut, each starting at its own monoisotopic peak (scan 26 of TaGe_SA_A549_3_snip)
            int charge = 2;
            double spacing = Constants.C13MinusC12 / charge;
            var left = GetLadderPrecursor(446.7105, 446.7105, 2, charge, 100);
            double rightMono = 446.7105 + 2 * spacing;
            var right = GetLadderPrecursor(rightMono, rightMono, 3, charge, 200);

            Assert.That(PrecursorSet.AreSequentialAtSpacing(left, right, spacing, new PpmTolerance(5)), Is.False);

            var precursors = new List<Precursor> { left, right };
            PrecursorSet.MergeSplitEnvelopes(in precursors, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(precursors, Has.Count.EqualTo(2));
            Assert.That(precursors, Does.Contain(left));
            Assert.That(precursors, Does.Contain(right));
        }

        [Test]
        public static void MergeSplitEnvelopes_LooksPastOverlappingNeighbour()
        {
            // L and R are halves of one split envelope; X is another species that overlaps L and sorts between them
            int charge = 2;
            double spacing = Constants.C13MinusC12 / charge;
            double m0 = 500.0;
            var left = GetLadderPrecursor(m0, m0, 6, charge, 300);
            var overlapping = GetLadderPrecursor(m0 + 0.5, m0 + 0.5, 2, charge, 50);
            var right = GetLadderPrecursor(m0 + 1.0, m0 + 6 * spacing, 3, charge, 100);

            var precursors = new List<Precursor> { left, overlapping, right };
            PrecursorSet.MergeSplitEnvelopes(in precursors, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(precursors, Has.Count.EqualTo(2));
            Assert.That(precursors, Does.Contain(overlapping));
            var merged = precursors.Single(p => !ReferenceEquals(p, overlapping));
            Assert.That(merged.EnvelopePeakCount, Is.EqualTo(9));
            Assert.That(merged.MonoisotopicPeakMz, Is.EqualTo(m0).Within(1e-6));
        }

        [Test]
        public static void MergeSplitEnvelopes_MergesEachRightOnce()
        {
            // Two rights share a first peak; after L takes one, the merged envelope no longer abuts the other
            int charge = 2;
            double spacing = Constants.C13MinusC12 / charge;
            double m0 = 500.0;
            var left = GetLadderPrecursor(m0, m0, 3, charge, 300);
            var right1 = GetLadderPrecursor(m0 + 0.9, m0 + 3 * spacing, 2, charge, 100);
            var right2 = GetLadderPrecursor(m0 + 1.1, m0 + 3 * spacing, 3, charge, 90);

            var precursors = new List<Precursor> { left, right1, right2 };
            PrecursorSet.MergeSplitEnvelopes(in precursors, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(precursors, Has.Count.EqualTo(2), "Only one of the two rights should be merged into left");
            Assert.That(precursors.Count(p => Math.Abs(p.MonoisotopicPeakMz - m0) < 1e-6), Is.EqualTo(1));
        }

        [Test]
        public static void MergePrecursors_CarriesOverBestDeconvolutionScore()
        {
            int charge = 3;
            double spacing = Constants.C13MinusC12 / charge;
            var left = GetLadderPrecursor(700, 700, 3, charge, 300, score: 0.8);
            var right = GetLadderPrecursor(700, 700 + 3 * spacing, 3, charge, 100, score: 0.4);

            var merged = PrecursorSet.MergePrecursors(left, right, charge, new PpmTolerance(5));

            Assert.That(merged.DeconvolutionScore, Is.EqualTo(0.8).Within(1e-9));
        }

        [Test]
        [TestCase(4, 2)]
        [TestCase(12, 4)]
        public static void FilterLowHarmonics_PeaksInApexOrder(int baseCharge, int harmonicCharge)
        {
            // Classic lists envelope peaks starting from the apex rather than in m/z order
            var tol = new PpmTolerance(5);
            string sequence = "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE";
            var mzOrdered = GetEnvelope(sequence, baseCharge, 100, 0.05);
            var harmonicMzOrdered = GetHarmonicEnvelope(mzOrdered, harmonicCharge);

            var apexOrdered = new IsotopicEnvelope(mzOrdered.Peaks.OrderByDescending(p => p.intensity).ToList(),
                mzOrdered.MonoisotopicMass, mzOrdered.Charge, mzOrdered.TotalIntensity, 0.5);
            var harmonicApexOrdered = new IsotopicEnvelope(harmonicMzOrdered.Peaks.OrderByDescending(p => p.intensity).ToList(),
                harmonicMzOrdered.MonoisotopicMass, harmonicMzOrdered.Charge, harmonicMzOrdered.TotalIntensity, 0.5);

            var basePrecursor = new Precursor(apexOrdered);
            var harmonicPrecursor = new Precursor(harmonicApexOrdered);
            var precursors = new List<Precursor> { basePrecursor, harmonicPrecursor };
            PrecursorSet.RemoveLowHarmonics(in precursors, tol);

            Assert.That(precursors, Does.Contain(basePrecursor));
            Assert.That(precursors, Does.Not.Contain(harmonicPrecursor));
        }

        [Test]
        public static void FilterLowHarmonics_RequiresMajorityOfLowPeaks()
        {
            // A z=2 species with 5 peaks, only 2 of which land on a 6-peak z=4 ladder, is not a harmonic
            var tol = new PpmTolerance(5);
            double spacing4 = Constants.C13MinusC12 / 4;
            double h0 = 800.0;
            var high = GetLadderPrecursor(h0, h0, 6, 4, 100);
            var low = GetLadderPrecursor(h0 + 2 * spacing4, h0 + 2 * spacing4, 5, 2, 100);

            var precursors = new List<Precursor> { high, low };
            PrecursorSet.RemoveLowHarmonics(in precursors, tol);

            Assert.That(precursors, Does.Contain(low));
        }

        private static Precursor GetPrecursorFromPeaks(List<(double mz, double intensity)> peaks, double monoMz, int charge,
            double? intensity = null, double? fractionalIntensity = null)
        {
            var envelope = new IsotopicEnvelope(peaks, monoMz.ToMass(charge), charge, peaks.Sum(p => p.intensity), 0.5);
            return new Precursor(envelope, intensity, fractionalIntensity);
        }

        [Test]
        public static void NonGenericEnumerator_ReturnsAllPrecursors()
        {
            var set = new PrecursorSet(new PpmTolerance(5))
            {
                new Precursor(500.0, 2, 100, 3),
                new Precursor(600.0, 2, 100, 3)
            };

            var enumerated = new List<object>();
            var enumerator = ((System.Collections.IEnumerable)set).GetEnumerator();
            while (enumerator.MoveNext())
                enumerated.Add(enumerator.Current);

            Assert.That(enumerated, Has.Count.EqualTo(2));
        }

        [Test]
        public static void RemoveDuplicates_KeepsTheBetterRepresentative()
        {
            var set = new PrecursorSet(new PpmTolerance(5));

            // One with an envelope beats one without, whatever the order
            var withEnvelope = GetLadderPrecursor(500.0, 500.0, 3, 2, 100);
            var withoutEnvelope = new Precursor(500.0, 2, 1000, 5);
            var list = new List<Precursor> { withoutEnvelope, withEnvelope };
            set.RemoveDuplicates(in list);
            Assert.That(list, Is.EquivalentTo(new[] { withEnvelope }));

            list = new List<Precursor> { withEnvelope, withoutEnvelope };
            set.RemoveDuplicates(in list);
            Assert.That(list, Is.EquivalentTo(new[] { withEnvelope }));

            // Same peak count: the more intense one wins
            var dim = new Precursor(500.0, 2, 100, 3);
            var bright = new Precursor(500.0, 2, 200, 3);
            list = new List<Precursor> { dim, bright };
            set.RemoveDuplicates(in list);
            Assert.That(list, Is.EquivalentTo(new[] { bright }));

            list = new List<Precursor> { bright, dim };
            set.RemoveDuplicates(in list);
            Assert.That(list, Is.EquivalentTo(new[] { bright }));

            // Same intensity: the larger fractional intensity wins
            var lowFraction = new Precursor(500.0, 2, 100, 3, 0.2);
            var highFraction = new Precursor(500.0, 2, 100, 3, 0.6);
            list = new List<Precursor> { lowFraction, highFraction };
            set.RemoveDuplicates(in list);
            Assert.That(list, Is.EquivalentTo(new[] { highFraction }));

            // Everything else equal: the lower m/z wins
            var lowerMz = new Precursor(500.0, 2, 100, 3, 0.5);
            var higherMz = new Precursor(500.001, 2, 100, 3, 0.5);
            list = new List<Precursor> { higherMz, lowerMz };
            set.RemoveDuplicates(in list);
            Assert.That(list, Is.EquivalentTo(new[] { lowerMz }));
        }

        [Test]
        public static void MergeSplitEnvelopes_ReusesABetterExistingPrecursor()
        {
            // An existing full envelope at the merged mass with more peaks is kept instead of the synthetic merge
            int charge = 2;
            double spacing = Constants.C13MinusC12 / charge;
            double m0 = 500.0;
            var left = GetLadderPrecursor(m0, m0, 3, charge, 300);
            var right = GetLadderPrecursor(m0 + 1.0, m0 + 3 * spacing, 2, charge, 100);
            var existing = GetLadderPrecursor(m0, m0, 7, charge, 300);

            var precursors = new List<Precursor> { left, existing, right };
            PrecursorSet.MergeSplitEnvelopes(in precursors, new PpmTolerance(5), Constants.C13MinusC12);

            Assert.That(precursors, Is.EquivalentTo(new[] { existing }));
        }

        [Test]
        public static void MergeAndHarmonics_IgnoreNullAndSingleItemLists()
        {
            var tol = new PpmTolerance(5);
            List<Precursor> nullList = null;
            Assert.DoesNotThrow(() => PrecursorSet.MergeSplitEnvelopes(in nullList, tol, Constants.C13MinusC12));
            Assert.DoesNotThrow(() => PrecursorSet.RemoveLowHarmonics(in nullList, tol));

            var single = new List<Precursor> { GetLadderPrecursor(500, 500, 3, 2, 100) };
            PrecursorSet.MergeSplitEnvelopes(in single, tol, Constants.C13MinusC12);
            PrecursorSet.RemoveLowHarmonics(in single, tol);
            Assert.That(single, Has.Count.EqualTo(1));
        }

        [Test]
        public static void MergePrecursors_OverlappingBoundaryPeakKeepsTheMoreIntenseOne()
        {
            var tol = new PpmTolerance(5);
            int charge = 2;
            double s = Constants.C13MinusC12 / charge;
            double a = 500.0;
            var left = GetPrecursorFromPeaks(new List<(double, double)> { (a, 100), (a + s, 100) }, a, charge);

            var brighterRight = GetPrecursorFromPeaks(new List<(double, double)> { (a + s, 200), (a + 2 * s, 50) }, a, charge);
            var merged = PrecursorSet.MergePrecursors(left, brighterRight, charge, tol);
            Assert.That(merged.EnvelopePeakCount, Is.EqualTo(3));
            Assert.That(merged.Envelope!.Peaks.Single(p => tol.Within(p.mz, a + s)).intensity, Is.EqualTo(200));

            var dimmerRight = GetPrecursorFromPeaks(new List<(double, double)> { (a + s, 50), (a + 2 * s, 50) }, a, charge);
            merged = PrecursorSet.MergePrecursors(left, dimmerRight, charge, tol);
            Assert.That(merged.EnvelopePeakCount, Is.EqualTo(3));
            Assert.That(merged.Envelope!.Peaks.Single(p => tol.Within(p.mz, a + s)).intensity, Is.EqualTo(100));
        }

        [Test]
        public static void MergePrecursors_TakesMassFromTheMoreIntenseSide()
        {
            var tol = new PpmTolerance(5);
            int charge = 2;
            double s = Constants.C13MinusC12 / charge;
            var left = GetLadderPrecursor(500, 500, 2, charge, 100);
            var right = GetLadderPrecursor(500.7, 500 + 2 * s, 2, charge, 400);

            var merged = PrecursorSet.MergePrecursors(left, right, charge, tol);

            Assert.That(merged.MonoisotopicPeakMz, Is.EqualTo(500.7).Within(1e-6));
            Assert.That(merged.Intensity, Is.EqualTo(1000).Within(1e-6), "Summed-intensity precursors keep summing");
        }

        [Test]
        public static void MergePrecursors_MostAbundantIntensityStaysMostAbundant()
        {
            var tol = new PpmTolerance(5);
            int charge = 2;
            double s = Constants.C13MinusC12 / charge;
            double a = 500.0;
            var leftPeaks = new List<(double, double)> { (a, 100), (a + s, 300) };
            var rightPeaks = new List<(double, double)> { (a + 2 * s, 200), (a + 3 * s, 50) };
            var left = GetPrecursorFromPeaks(leftPeaks, a, charge, intensity: 300);
            var right = GetPrecursorFromPeaks(rightPeaks, a, charge, intensity: 200);

            var merged = PrecursorSet.MergePrecursors(left, right, charge, tol);

            Assert.That(merged.Intensity, Is.EqualTo(300).Within(1e-6));
        }

        [Test]
        [TestCase(0.2, 0.3, 0.5)]
        [TestCase(0.2, null, 0.2)]
        [TestCase(null, 0.3, 0.3)]
        [TestCase(null, null, null)]
        public static void MergePrecursors_CombinesFractionalIntensity(double? leftFraction, double? rightFraction, double? expected)
        {
            var tol = new PpmTolerance(5);
            int charge = 2;
            double s = Constants.C13MinusC12 / charge;
            double a = 500.0;
            var left = GetPrecursorFromPeaks(new List<(double, double)> { (a, 100), (a + s, 100) }, a, charge, fractionalIntensity: leftFraction);
            var right = GetPrecursorFromPeaks(new List<(double, double)> { (a + 2 * s, 50), (a + 3 * s, 50) }, a, charge, fractionalIntensity: rightFraction);

            var merged = PrecursorSet.MergePrecursors(left, right, charge, tol);

            if (expected is null)
                Assert.That(merged.FractionalIntensity, Is.Null);
            else
                Assert.That(merged.FractionalIntensity, Is.EqualTo(expected.Value).Within(1e-9));
        }

        #endregion
    }
}
