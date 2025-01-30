using EngineLayer;
using NUnit.Framework;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;

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

        [Test]
        public void Increment_IncrementsExistingKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(2));
        }

        [Test]
        public void Increment_AddsNewKeyWithInitialValue()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(1));
        }

        [Test]
        public void Increment_IncrementsMultipleTimes()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            dictionary.Increment("key1");
            dictionary.Increment("key1");
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(4));
        }

        [Test]
        public void Increment_AddsAndIncrementsNewKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            dictionary.Increment("key1");
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(2));
        }

        [Test]
        public void Increment_ThreadSafeWithConcurrentDictionary()
        {
            // Arrange
            var dictionary = new ConcurrentDictionary<string, int>();
            var tasks = new List<Task>();

            // Act
            for (int i = 0; i < 1000; i++)
            {
                tasks.Add(Task.Run(() => dictionary.Increment("key1")));
            }
            Task.WaitAll(tasks.ToArray());

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(1000));
        }

        [Test]
        public void AddOrCreate_ThreadSafeWithConcurrentDictionary()
        {
            // Arrange
            var dictionary = new ConcurrentDictionary<string, IList<int>>();
            var tasks = new List<Task>();
            var lockObject = new object();

            // Act
            for (int i = 0; i < 1000; i++)
            {
                int value = i;
                tasks.Add(Task.Run(() => dictionary.AddOrCreateThreadSafe("key1", value, lockObject)));
            }
            Task.WaitAll(tasks.ToArray());

            // Assert
            Assert.That(dictionary["key1"].Count, Is.EqualTo(1000));
            for (int i = 0; i < 1000; i++)
            {
                Assert.That(dictionary["key1"], Contains.Item(i));
            }
        }

        [Test]
        public void AddOrCreate_ThreadSafeWithDictionary()
        {
            // Arrange
            var dictionary = new Dictionary<string, IList<int>>();
            var tasks = new List<Task>();
            var lockObject = new object();

            // Act
            for (int i = 0; i < 1000; i++)
            {
                int value = i;
                tasks.Add(Task.Run(() => dictionary.AddOrCreateThreadSafe("key1", value, lockObject)));
            }
            Task.WaitAll(tasks.ToArray());

            // Assert
            Assert.That(dictionary["key1"].Count, Is.EqualTo(1000));
            for (int i = 0; i < 1000; i++)
            {
                Assert.That(dictionary["key1"], Contains.Item(i));
            }
        }

        [Test]
        public void IsNullOrEmpty_ReturnsTrueForNullDictionary()
        {
            // Arrange
            Dictionary<string, int> dictionary = null;

            // Act
            var result = dictionary.IsNullOrEmpty();

            // Assert
            Assert.That(result, Is.True);
        }

        [Test]
        public void IsNullOrEmpty_ReturnsTrueForEmptyDictionary()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            var result = dictionary.IsNullOrEmpty();

            // Assert
            Assert.That(result, Is.True);
        }

        [Test]
        public void IsNullOrEmpty_ReturnsFalseForNonEmptyDictionary()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            var result = dictionary.IsNullOrEmpty();

            // Assert
            Assert.That(result, Is.False);
        }
    }
}
