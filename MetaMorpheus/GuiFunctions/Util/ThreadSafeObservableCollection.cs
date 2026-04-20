using System.Collections;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Collections.Specialized;
using System.Linq;
using System.Windows;
using System.Windows.Threading;

namespace GuiFunctions.Util;

/// <summary>
/// An ObservableCollection that is safe for concurrent reads and writes.
/// All mutations are serialized with a lock, and GetEnumerator returns a
/// snapshot so callers never see a "collection was modified" exception.
/// CollectionChanged notifications are automatically marshalled to the UI
/// dispatcher when raised from a background thread.
/// </summary>
public class ThreadSafeObservableCollection<T> : ObservableCollection<T>, IEnumerable<T>, IEnumerable
{
    private readonly object _lock = new();

    public ThreadSafeObservableCollection() { }

    public ThreadSafeObservableCollection(IEnumerable<T> collection) : base(collection) { }

    // ── Mutations ────────────────────────────────────────────────────────────

    protected override void InsertItem(int index, T item)
    {
        lock (_lock) base.InsertItem(index, item);
    }

    protected override void RemoveItem(int index)
    {
        lock (_lock) base.RemoveItem(index);
    }

    protected override void SetItem(int index, T item)
    {
        lock (_lock) base.SetItem(index, item);
    }

    protected override void MoveItem(int oldIndex, int newIndex)
    {
        lock (_lock) base.MoveItem(oldIndex, newIndex);
    }

    protected override void ClearItems()
    {
        lock (_lock) base.ClearItems();
    }

    // ── Snapshot enumeration ─────────────────────────────────────────────────
    // Explicit interface implementation ensures LINQ (which dispatches through
    // IEnumerable<T>) gets the snapshot enumerator, not the base class one.

    IEnumerator<T> IEnumerable<T>.GetEnumerator()
    {
        T[] snapshot;
        lock (_lock)
            snapshot = Items.ToArray(); // Items is the protected List<T> from Collection<T>
        return ((IEnumerable<T>)snapshot).GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator() =>
        ((IEnumerable<T>)this).GetEnumerator();

    // ── UI-thread marshalling ─────────────────────────────────────────────────

    protected override void OnCollectionChanged(NotifyCollectionChangedEventArgs e)
    {
        Dispatcher dispatcher = Application.Current?.Dispatcher;
        if (dispatcher != null && !dispatcher.CheckAccess())
            dispatcher.Invoke(() => base.OnCollectionChanged(e));
        else
            base.OnCollectionChanged(e);
    }
}
