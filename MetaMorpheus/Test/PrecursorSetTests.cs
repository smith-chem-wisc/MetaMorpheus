using Easy.Common.Extensions;
using EngineLayer.Util;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

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
        public void Add_DuplicatePrecursor_ReturnsFalse()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5, 2, 999.0, 2000.0, 4); // Same m/z and charge, different intensity

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.False);
            Assert.That(_precursorSet.Count, Is.EqualTo(1));
        }

        [Test]
        public void Add_EquivalentPrecursorWithinTolerance_ReturnsFalse()
        {
            // Arrange
            var precursor1 = new Precursor(500.5, 2, 999.0, 1000.0, 3);
            var precursor2 = new Precursor(500.5 + (_tolerance.Value * 500.5 / 1e6 * 0.5), 2, 999.0, 2000.0, 4); // Within tolerance

            // Act
            bool result1 = _precursorSet.Add(precursor1);
            bool result2 = _precursorSet.Add(precursor2);

            // Assert
            Assert.That(result1, Is.True);
            Assert.That(result2, Is.False);
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
            Assert.That(_precursorSet.Count, Is.EqualTo(3), "Count should not change when adding duplicate");
        }
    }
}