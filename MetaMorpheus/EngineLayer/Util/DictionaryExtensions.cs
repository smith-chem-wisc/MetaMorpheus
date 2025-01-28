using System.Collections.Generic;

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
}