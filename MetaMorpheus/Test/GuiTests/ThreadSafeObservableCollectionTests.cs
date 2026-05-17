using GuiFunctions.Util;
using NUnit.Framework;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Test.GuiTests;

[TestFixture]
public class ThreadSafeObservableCollectionTests
{
    [Test]
    public void DefaultConstructor_StartsEmpty()
    {
        var c = new ThreadSafeObservableCollection<int>();
        Assert.That(c.Count, Is.EqualTo(0));
    }

    [Test]
    public void Constructor_WithInitialCollection_CopiesItems()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2, 3 });
        Assert.That(c.Count, Is.EqualTo(3));
        Assert.That(c.ToList(), Is.EqualTo(new[] { 1, 2, 3 }));
    }

    // Plain foreach on the concrete type uses the snapshot enumerator (via the `new` GetEnumerator),
    // so mutation mid-iteration does not throw and is not observed by the in-flight enumeration.
    [Test]
    public void SnapshotEnumeration_MutationDuringEnumerationDoesNotThrow()
    {
        var c = new ThreadSafeObservableCollection<int>(Enumerable.Range(0, 10));

        var enumerated = new List<int>();
        foreach (var item in c)
        {
            enumerated.Add(item);
            if (item == 5) c.Add(999);
        }

        Assert.That(enumerated, Is.EqualTo(Enumerable.Range(0, 10).ToList()));
        Assert.That(c.Count, Is.EqualTo(11));
    }

    // Pins the snapshot guarantee at the GetEnumerator call site on the concrete type.
    [Test]
    public void GetEnumerator_OnConcreteType_ReturnsSnapshot()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2, 3 });
        var e = c.GetEnumerator();
        c.Add(4);

        var enumerated = new List<int>();
        while (e.MoveNext()) enumerated.Add(e.Current);

        Assert.That(enumerated, Is.EqualTo(new[] { 1, 2, 3 }));
        Assert.That(c.Count, Is.EqualTo(4));
    }

    // The explicit IEnumerable<T>.GetEnumerator is what LINQ uses; verify a representative LINQ op
    // returns the snapshot view rather than a live one.
    [Test]
    public void Linq_ToListReturnsSnapshot()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2, 3 });
        var snapshot = c.ToList();
        c.Add(4);
        Assert.That(snapshot, Is.EqualTo(new[] { 1, 2, 3 }));
        Assert.That(c.ToList(), Is.EqualTo(new[] { 1, 2, 3, 4 }));
    }

    // Concurrent Adds must all be preserved (no lost writes from races on the base list).
    [Test]
    public void ConcurrentAdds_AllItemsPreserved()
    {
        const int threads = 8;
        const int itemsPerThread = 500;
        var c = new ThreadSafeObservableCollection<int>();

        Parallel.For(0, threads, t =>
        {
            for (int i = 0; i < itemsPerThread; i++)
                c.Add(t * itemsPerThread + i);
        });

        Assert.That(c.Count, Is.EqualTo(threads * itemsPerThread));
        // Every produced value should appear exactly once.
        Assert.That(c.ToList().Distinct().Count(), Is.EqualTo(threads * itemsPerThread));
    }

    // Enumeration on one thread while another thread mutates must not throw.
    [Test]
    public void ConcurrentEnumerationDuringMutation_DoesNotThrow()
    {
        var c = new ThreadSafeObservableCollection<int>(Enumerable.Range(0, 100));
        var cts = new CancellationTokenSource();

        var mutator = Task.Run(() =>
        {
            int i = 1000;
            while (!cts.IsCancellationRequested)
            {
                c.Add(i++);
                if (c.Count > 200) c.RemoveAt(0);
            }
        });

        var enumerator = Task.Run(() =>
        {
            for (int i = 0; i < 1000; i++)
                _ = c.Sum();
        });

        Assert.DoesNotThrow(() => enumerator.Wait());
        cts.Cancel();
        mutator.Wait();
    }

    [Test]
    public void CollectionChanged_FiresOnAdd()
    {
        var c = new ThreadSafeObservableCollection<int>();
        NotifyCollectionChangedAction? lastAction = null;
        c.CollectionChanged += (_, e) => lastAction = e.Action;

        c.Add(1);

        Assert.That(lastAction, Is.EqualTo(NotifyCollectionChangedAction.Add));
    }

    [Test]
    public void CollectionChanged_FiresOnRemove()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2, 3 });
        NotifyCollectionChangedAction? lastAction = null;
        c.CollectionChanged += (_, e) => lastAction = e.Action;

        c.Remove(2);

        Assert.That(lastAction, Is.EqualTo(NotifyCollectionChangedAction.Remove));
    }

    [Test]
    public void CollectionChanged_FiresOnClear()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2 });
        NotifyCollectionChangedAction? lastAction = null;
        c.CollectionChanged += (_, e) => lastAction = e.Action;

        c.Clear();

        Assert.That(lastAction, Is.EqualTo(NotifyCollectionChangedAction.Reset));
    }

    [Test]
    public void CollectionChanged_FiresOnReplace()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2, 3 });
        NotifyCollectionChangedAction? lastAction = null;
        c.CollectionChanged += (_, e) => lastAction = e.Action;

        c[1] = 99;

        Assert.That(lastAction, Is.EqualTo(NotifyCollectionChangedAction.Replace));
        Assert.That(c[1], Is.EqualTo(99));
    }

    [Test]
    public void CollectionChanged_FiresOnMove()
    {
        var c = new ThreadSafeObservableCollection<int>(new[] { 1, 2, 3 });
        NotifyCollectionChangedAction? lastAction = null;
        c.CollectionChanged += (_, e) => lastAction = e.Action;

        c.Move(0, 2);

        Assert.That(lastAction, Is.EqualTo(NotifyCollectionChangedAction.Move));
        Assert.That(c.ToList(), Is.EqualTo(new[] { 2, 3, 1 }));
    }
}
