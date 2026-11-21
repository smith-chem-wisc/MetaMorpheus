using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using NUnit.Framework;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;

namespace Test.ParallelSearchTask.Analysis;

[TestFixture]
public class CollectorColumnAndMetricsTests
{
    private static IEnumerable<IMetricCollector> AllCollectors()
    {
        yield return new BasicMetricCollector();
        yield return new PsmPeptideSearchCollector("Homo sapiens");
        yield return new FragmentIonCollector();
        yield return new RetentionTimeCollector();
        yield return new ProteinGroupCollector("Homo sapiens");
    }

    [TestCaseSource(nameof(AllCollectors))]
    public void GetOutputColumns_AreNamedAndUniqueWithinTheCollector(IMetricCollector collector)
    {
        var columns = collector.GetOutputColumns().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(collector.CollectorName, Is.Not.Empty);
            Assert.That(columns, Is.Not.Empty);
            Assert.That(columns, Has.None.Empty);
            Assert.That(columns, Is.Unique);
        });
    }

    [Test]
    public void GetOutputColumns_DoNotCollideAcrossCollectors()
    {
        var columns = AllCollectors().SelectMany(c => c.GetOutputColumns()).ToList();

        Assert.That(columns, Is.Unique);
    }

    [Test]
    public void GetValue_ReturnsTypedConvertedOrDefaultValues()
    {
        var metrics = new TransientDatabaseMetrics("Db");
        metrics.Results["int"] = 3;
        metrics.Results["text"] = "not a number";

        Assert.Multiple(() =>
        {
            Assert.That(metrics.GetValue<int>("int"), Is.EqualTo(3));
            Assert.That(metrics.GetValue<double>("int"), Is.EqualTo(3.0));
            Assert.That(metrics.GetValue("text", -1), Is.EqualTo(-1));
            Assert.That(metrics.GetValue("missing", 7), Is.EqualTo(7));
        });
    }

    [Test]
    public void Equals_ComparesDatabaseNameAndResultCount()
    {
        var metrics = new TransientDatabaseMetrics("Db");
        metrics.Results["a"] = 1;
        var same = new TransientDatabaseMetrics("Db");
        same.Results["a"] = 2;
        var otherName = new TransientDatabaseMetrics("Other");
        otherName.Results["a"] = 1;
        var otherCount = new TransientDatabaseMetrics("Db");

        Assert.Multiple(() =>
        {
            Assert.That(metrics.Equals(same), Is.True);
            Assert.That(metrics.GetHashCode(), Is.EqualTo(same.GetHashCode()));
            Assert.That(metrics.Equals(otherName), Is.False);
            Assert.That(metrics.Equals(otherCount), Is.False);
            Assert.That(metrics.Equals(null), Is.False);
        });
    }

    [Test]
    public void Add_SumsNumbersConcatenatesArraysAndKeepsUnmatchedKeys()
    {
        var left = new TransientDatabaseMetrics("Db");
        left.Results["count"] = 2;
        left.Results["mean"] = 1.5;
        left.Results["values"] = new[] { 1.0, 2.0 };
        left.Results["onlyLeft"] = 9;

        var right = new TransientDatabaseMetrics("Db");
        right.Results["count"] = 3;
        right.Results["mean"] = 2.5;
        right.Results["values"] = new[] { 3.0 };

        var sum = left.Add(right);

        Assert.Multiple(() =>
        {
            Assert.That(sum.DatabaseName, Is.EqualTo("Db"));
            Assert.That(sum.Results["count"], Is.EqualTo(5));
            Assert.That(sum.Results["mean"], Is.EqualTo(4.0));
            Assert.That(sum.Results["values"], Is.EqualTo(new[] { 1.0, 2.0, 3.0 }));
            Assert.That(sum.Results["onlyLeft"], Is.EqualTo(9));
            Assert.That(left.Add(null!), Is.SameAs(left));
        });
    }

    [Test]
    public void Add_RejectsADifferentDatabaseOrAnUnsupportedType()
    {
        var left = new TransientDatabaseMetrics("Db");
        left.Results["text"] = "a";
        var sameNameText = new TransientDatabaseMetrics("Db");
        sameNameText.Results["text"] = "b";

        Assert.Multiple(() =>
        {
            Assert.That(() => left.Add(new TransientDatabaseMetrics("Other")), Throws.TypeOf<MetaMorpheusException>());
            Assert.That(() => left.Add(sameNameText), Throws.TypeOf<MetaMorpheusException>());
        });
    }
}
