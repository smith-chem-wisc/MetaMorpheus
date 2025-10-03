using GuiFunctions.Util;
using NUnit.Framework;

namespace Test.GuiTests;

[TestFixture]
public class CyclicalQueueTests
{
    [Test]
    public void Constructor_ThrowsOnNull()
    {
        Assert.That(() => new CyclicalQueue<int>(null), Throws.ArgumentNullException);
    }

    [Test]
    public void Constructor_ThrowsOnEmpty()
    {
        Assert.That(() => new CyclicalQueue<int>(new int[0]), Throws.ArgumentException);
    }

    [Test]
    public void Constructor_WorksWithValidInput()
    {
        var queue = new CyclicalQueue<string>(new[] { "a", "b" });
        Assert.That(queue.Count, Is.EqualTo(2));
    }

    [Test]
    public void Dequeue_CyclesThroughAllItems()
    {
        var items = new[] { 1, 2, 3 };
        var queue = new CyclicalQueue<int>(items);

        Assert.That(queue.Dequeue(), Is.EqualTo(1));
        Assert.That(queue.Dequeue(), Is.EqualTo(2));
        Assert.That(queue.Dequeue(), Is.EqualTo(3));
        // Should cycle back to start
        Assert.That(queue.Dequeue(), Is.EqualTo(1));
        Assert.That(queue.Dequeue(), Is.EqualTo(2));
    }

    [Test]
    public void Reset_SetsIndexToZero()
    {
        var queue = new CyclicalQueue<string>(new[] { "x", "y" });
        queue.Dequeue(); // "x"
        queue.Dequeue(); // "y"
        queue.Reset();
        Assert.That(queue.Dequeue(), Is.EqualTo("x"));
    }

    [Test]
    public void Count_ReturnsCorrectCount()
    {
        var queue = new CyclicalQueue<double>(new[] { 1.1, 2.2, 3.3 });
        Assert.That(queue.Count, Is.EqualTo(3));
    }

    [Test]
    public void Dequeue_OneItem_AlwaysReturnsThatItem()
    {
        var queue = new CyclicalQueue<int>(new[] { 42 });
        Assert.That(queue.Dequeue(), Is.EqualTo(42));
        Assert.That(queue.Dequeue(), Is.EqualTo(42));
        Assert.That(queue.Dequeue(), Is.EqualTo(42));
    }

    [Test]
    public static void Dequeue_ResetsOnCollectionBeingEmpty()
    {
        var queue = new CyclicalQueue<int>(new[] { 1, 2, 3 });
        Assert.That(queue.Dequeue(), Is.EqualTo(1));
        Assert.That(queue.Dequeue(), Is.EqualTo(2));
        Assert.That(queue.Dequeue(), Is.EqualTo(3));

        // Collection is now "empty" in the sense that we've cycled through all items
        Assert.That(queue.Dequeue(), Is.EqualTo(1)); // Should reset and return first item
        Assert.That(queue.Dequeue(), Is.EqualTo(2));
        Assert.That(queue.Dequeue(), Is.EqualTo(3));


        Assert.That(queue.Dequeue(), Is.EqualTo(1));
        Assert.That(queue.Dequeue(), Is.EqualTo(2));
        Assert.That(queue.Dequeue(), Is.EqualTo(3));
    }   
}
