using Easy.Common.Extensions;
using EngineLayer.Util;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;

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
            Assert.That(_precursorSet.Tolerance, Is.EqualTo(_tolerance));
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
            Assert.That(p1.Mass, Is.EqualTo(p2.Mass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p2.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p2.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p2.FractionalIntensity));

            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p3.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p3.Charge));
            Assert.That(p1.Mass, Is.EqualTo(p3.Mass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p3.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p3.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p3.FractionalIntensity));

            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p4.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p4.Charge));
            Assert.That(p1.Mass, Is.EqualTo(p4.Mass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p4.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p4.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p4.FractionalIntensity));

            Assert.That(p1.MonoisotopicPeakMz, Is.EqualTo(p5.MonoisotopicPeakMz));
            Assert.That(p1.Charge, Is.EqualTo(p5.Charge));
            Assert.That(p1.Mass, Is.EqualTo(p5.Mass).Within(1e-6));
            Assert.That(p1.Intensity, Is.EqualTo(p5.Intensity).Within(1e-6));
            Assert.That(p1.EnvelopePeakCount, Is.EqualTo(p5.EnvelopePeakCount));
            Assert.That(p1.FractionalIntensity, Is.EqualTo(p5.FractionalIntensity));
        }
    }
}