using System;
using System.Collections.Generic;

namespace GuiFunctions.Util;

public class CyclicalQueue<T>
{
    private readonly List<T> _items;
    private int _index;

    public CyclicalQueue(IEnumerable<T> items)
    {
        if (items == null) throw new ArgumentNullException(nameof(items));
        _items = new List<T>(items);
        if (_items.Count == 0)
            throw new ArgumentException("CyclicalQueue must be initialized with at least one item.");
        _index = 0;
    }

    public T Dequeue()
    {
        var item = _items[_index];
        _index = (_index + 1) % _items.Count;
        return item;
    }

    public void Reset()
    {
        _index = 0;
    }

    public int Count => _items.Count;
}
