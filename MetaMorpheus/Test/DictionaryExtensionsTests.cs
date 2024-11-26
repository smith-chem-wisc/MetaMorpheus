using NUnit.Framework;
using System.Collections.Generic;
using EngineLayer;

namespace Test
{
    [TestFixture]
    public class DictionaryExtensionsTests
    {
        [Test]
        public void AddOrCreate_AddsValueToExistingKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, List<int>>
            {
                { "key1", new List<int> { 1, 2 } }
            };

            // Act
            dictionary.AddOrCreate("key1", 3);

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(new List<int> { 1, 2, 3 }));
        }

        [Test]
        public void AddOrCreate_CreatesNewKeyWithValue()
        {
            // Arrange
            var dictionary = new Dictionary<string, List<int>>();

            // Act
            dictionary.AddOrCreate("key1", 1);

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(new List<int> { 1 }));
        }

        [Test]
        public void AddOrCreate_AddsMultipleValuesToExistingKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, List<int>>
            {
                { "key1", new List<int> { 1 } }
            };

            // Act
            dictionary.AddOrCreate("key1", 2);
            dictionary.AddOrCreate("key1", 3);

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(new List<int> { 1, 2, 3 }));
        }

        [Test]
        public void AddOrCreate_CreatesMultipleKeysWithValues()
        {
            // Arrange
            var dictionary = new Dictionary<string, List<int>>();

            // Act
            dictionary.AddOrCreate("key1", 1);
            dictionary.AddOrCreate("key2", 2);

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(new List<int> { 1 }));
            Assert.That(dictionary.ContainsKey("key2"));
            Assert.That(dictionary["key2"], Is.EqualTo(new List<int> { 2 }));
        }
    }
}
