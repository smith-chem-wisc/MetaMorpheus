using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public class DictionaryExtensionsTests
    {
        [Test]
        public void AddOrCreate_AddsValueToExistingKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, IList<int>>
            {
                { "key1", new List<int> { 1, 2 } }
            };

            // Act
            dictionary.AddOrCreate("key1", 3);

            // Assert
            Assert.That(dictionary["key1"], Is.EquivalentTo(new List<int> { 1, 2, 3 }));
        }

        [Test]
        public void AddOrCreate_CreatesNewKeyWithList()
        {
            // Arrange
            var dictionary = new Dictionary<string, IList<int>>();

            // Act
            dictionary.AddOrCreate("key1", 1);

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EquivalentTo(new List<int> { 1 }));
        }

        [Test]
        public void AddOrCreate_AddsMultipleValuesToNewKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, IList<int>>();

            // Act
            dictionary.AddOrCreate("key1", 1);
            dictionary.AddOrCreate("key1", 2);
            dictionary.AddOrCreate("key1", 3);

            // Assert
            Assert.That(dictionary["key1"], Is.EquivalentTo(new List<int> { 1, 2, 3 }));
        }

        [Test]
        public void AddOrCreate_AddsMultipleValuesToExistingKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, IList<int>>
            {
                { "key1", new List<int> { 1 } }
            };

            // Act
            dictionary.AddOrCreate("key1", 2);
            dictionary.AddOrCreate("key1", 3);

            // Assert
            Assert.That(dictionary["key1"], Is.EquivalentTo(new List<int> { 1, 2, 3 }));
        }
    }
}
