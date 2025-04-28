using Nett;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Numerics;

namespace EngineLayer;

public static class DictionaryExtensions
{
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
    /// Adds a value to the HashSet associated with the specified key in the dictionary.
    /// If the key does not exist, a new HashSet is created with the value and added to the dictionary.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValues">The type of the values in the hash set.</typeparam>
    /// <param name="dictionary">The dictionary to operate on.</param>
    /// <param name="key">The key whose value list to add to or create.</param>
    /// <param name="value">The value to add to the list associated with the specified key.</param>
    /// <remarks>
    /// This is not thread-safe!
    /// </remarks>
    public static void AddOrCreate<TKey, TValues>(this IDictionary<TKey, HashSet<TValues>> dictionary, TKey key, TValues value)
    {
        if (dictionary.TryGetValue(key, out HashSet<TValues> values))
        {
            values.Add(value);
        }
        else
        {
            dictionary.Add(key, new HashSet<TValues> { value });
        }
    }

    /// <summary>
    /// Merges two HashSets associated with the specified key in the dictionary.
    /// If the key does not exist, a new HashSet is created with the value and added to the dictionary.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValues">The type of the values in the HashSets.</typeparam>
    /// <param name="dictionary">The dictionary to operate on.</param>
    /// <param name="key">The key whose value HashSet to add to or create.</param>
    /// <param name="value">The HashSet of values associated with the specified key.</param>
    /// <remarks>
    /// This is not thread-safe!
    /// </remarks>
    public static void MergeOrCreate<TKey, TValues>(this IDictionary<TKey, HashSet<TValues>> dictionary, TKey key, HashSet<TValues> value)
    {
        if (dictionary.TryGetValue(key, out HashSet<TValues> values))
        {
            values.UnionWith(value);
        }
        else
        {
            dictionary.Add(key, value);
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
    /// <param name="lockObject">Object used to lock this specific dictionary</param>
    /// <remarks>
    /// This is thread safe for all dictionary types. 
    /// </remarks>
    public static void AddOrCreateThreadSafe<TKey, TValues>(this IDictionary<TKey, IList<TValues>> dictionary, TKey key, TValues value, object lockObject)
    {
        if (dictionary is ConcurrentDictionary<TKey, IList<TValues>> concurrentDictionary)
        {
            concurrentDictionary.AddOrUpdate(key, new List<TValues> { value }, (k, v) =>
            {
                // must lock inside the add or update as the List.Add method is not thread safe. 
                lock (lockObject)
                {
                    v.Add(value);
                    return v;
                }
            });
        }
        else
        {
            lock (lockObject)
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
    /// Increments the value associated with the specified key in the dictionary by a defined number with a default of one.
    /// If the key does not exist, a new entry is created with the value set to one.
    /// </summary>
    /// <typeparam name="TKey">The type of the keys in the dictionary.</typeparam>
    /// <typeparam name="TValue">The type of the values in the dictionary, which must implement <see cref="INumber{TValue}"/>.</typeparam>
    /// <param name="dictionary">The dictionary to operate on.</param>
    /// <param name="key">The key whose value to increment or create.</param>
    /// <param name="incrementBy">The amount to increment by with a default of one</param>
    public static void Increment<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue incrementBy = default)
    where TValue : INumber<TValue>
    {
        TValue incrementValue = incrementBy.Equals(default) ? TValue.One : incrementBy;

        if (dictionary is ConcurrentDictionary<TKey, TValue> concurrentDictionary)
        {
            concurrentDictionary.AddOrUpdate(key, incrementValue, (k, v) => v + incrementValue);
        }
        else
        {
            if (dictionary.TryGetValue(key, out TValue value))
            {
                dictionary[key] = value + incrementValue;
            }
            else
            {
                dictionary.Add(key, incrementValue);
            }
        }
    }
}
