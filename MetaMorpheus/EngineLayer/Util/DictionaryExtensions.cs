using Nett;
using Newtonsoft.Json.Linq;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Numerics;

namespace EngineLayer;

public static class DictionaryExtensions
{
    private static readonly object AddOrCreateLock = new object();

    /// <summary>
    /// Adds a value to the list associated with the specified key in the dictionary.
    /// If the key does not exist, a new list is created with the value and added to the dictionary.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValues">The type of the values in the lists.</typeparam>
    /// <param name="dictionary">The dictionary to operate on.</param>
    /// <param name="key">The key whose value list to add to or create.</param>
    /// <param name="value">The value to add to the list associated with the specified key.</param>
    /// <remarks>
    /// This is not thread-safe!
    /// </remarks>
    public static void AddOrCreate<TKey, TValues>(this IDictionary<TKey, IList<TValues>> dictionary, TKey key, TValues value)
    {
        if (dictionary.TryGetValue(key, out IList<TValues> values))
        {
            values.Add(value);
        }
        else
        {
            dictionary.Add(key, new List<TValues> { value });
        }
    }

    /// <summary>
    /// Adds a value to the list associated with the specified key in the dictionary.
    /// If the key does not exist, a new list is created with the value and added to the dictionary.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValues">The type of the values in the lists.</typeparam>
    /// <param name="dictionary">The dictionary to operate on.</param>
    /// <param name="key">The key whose value list to add to or create.</param>
    /// <param name="value">The value to add to the list associated with the specified key.</param>
    /// <remarks>
    /// This is thread safe for all dictionary types. 
    /// </remarks>
    public static void AddOrCreateThreadSafe<TKey, TValues>(this IDictionary<TKey, IList<TValues>> dictionary, TKey key, TValues value)
    {
        if (dictionary is ConcurrentDictionary<TKey, IList<TValues>> concurrentDictionary)
        {
            concurrentDictionary.AddOrUpdate(key, new List<TValues> { value }, (k, v) =>
            {
                lock (AddOrCreateLock)
                {
                    v.Add(value);
                    return v;
                }
            });
        }
        else
        {
            lock (AddOrCreateLock)
            {
                if (dictionary.TryGetValue(key, out IList<TValues> values))
                {
                    values.Add(value);
                }
                else
                {
                    dictionary.Add(key, new List<TValues> { value });
                }
            }
        }
    }

    /// <summary>
    /// Increments the value associated with the specified key in the dictionary.
    /// If the key does not exist, a new entry is created with the value set to one.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValue">The type of the values in the dictionary, which must implement <see cref="INumber{TValue}"/>.</typeparam>
    /// <param name="dictionary">The dictionary to operate on.</param>
    /// <param name="key">The key whose value to increment or create.</param>
    public static void Increment<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key)
        where TValue : INumber<TValue>
    {
        if (dictionary is ConcurrentDictionary<TKey, TValue> concurrentDictionary)
        {
            concurrentDictionary.AddOrUpdate(key, TValue.One, (k, v) => v + TValue.One);
        }
        else
        {
            if (dictionary.TryGetValue(key, out TValue value))
            {
                dictionary[key] = value + TValue.One;
            }
            else
            {
                dictionary.Add(key, TValue.One);
            }
        }
    }

    /// <summary>
    /// Determines whether the dictionary is null or has no elements.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValue">The type of the values in the dictionary.</typeparam>
    /// <param name="dictionary">The dictionary to check.</param>
    /// <returns>
    /// <c>true</c> if the dictionary is null or has no elements; otherwise, <c>false</c>.
    /// </returns>
    public static bool IsNullOrEmpty<TKey, TValue>(this IDictionary<TKey, TValue> dictionary)
    {
        return dictionary == null || dictionary.Count == 0;
    }
}
