using System.Collections.Generic;

namespace EngineLayer;

public static class DictionaryExtensions
{
    public static void AddOrCreate<TKey, TValue>(this IDictionary<TKey, List<TValue>> dictionary, TKey key, TValue value)
    {
        if (dictionary.TryGetValue(key, out var list))
        {
            list.Add(value);
        }
        else
        {
            dictionary[key] = [value];
        }
    }
}