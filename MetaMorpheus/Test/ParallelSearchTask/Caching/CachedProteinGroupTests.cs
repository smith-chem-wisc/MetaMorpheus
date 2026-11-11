using System;
using System.Reflection;
using NUnit.Framework;
using Test.ParallelSearchTask.Utility;

namespace Test.ParallelSearchTask.Caching;

[TestFixture]
public class CachedProteinGroupTests
{
    [Test]
    public void Constructor_CreatesDefensiveCopyFromOriginalProteinGroup()
    {
        var original = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.02, peptideCount: 2);
        original.ProteinGroupScore = 123.45;
        original.BestPeptideScore = 88.9;
        original.CumulativeTarget = 10;
        original.CumulativeDecoy = 1;

        object cachedProteinGroup = CreateCachedProteinGroup(original);

        original.Proteins.Clear();
        original.AllPeptides.Clear();
        original.UniquePeptides.Clear();
        original.ProteinGroupScore = 1.0;

        var runtimeCopy = CreateRuntimeCopy(cachedProteinGroup);

        Assert.Multiple(() =>
        {
            Assert.That(runtimeCopy.Proteins.Count, Is.EqualTo(2));
            Assert.That(runtimeCopy.AllPeptides.Count, Is.EqualTo(2));
            Assert.That(runtimeCopy.UniquePeptides.Count, Is.EqualTo(2));
            Assert.That(runtimeCopy.ProteinGroupScore, Is.EqualTo(123.45).Within(1e-10));
            Assert.That(runtimeCopy.BestPeptideScore, Is.EqualTo(88.9).Within(1e-10));
            Assert.That(runtimeCopy.CumulativeTarget, Is.EqualTo(10));
            Assert.That(runtimeCopy.CumulativeDecoy, Is.EqualTo(1));
            Assert.That(runtimeCopy.QValue, Is.EqualTo(0.02).Within(1e-10));
        });
    }

    [Test]
    public void CreateRuntimeCopy_ReturnsIndependentInstances()
    {
        var original = ParallelSearchTestContextFactory.CreateProteinGroup(isDecoy: false, qValue: 0.001, peptideCount: 2);
        object cachedProteinGroup = CreateCachedProteinGroup(original);

        var firstCopy = CreateRuntimeCopy(cachedProteinGroup);
        var secondCopy = CreateRuntimeCopy(cachedProteinGroup);

        firstCopy.Proteins.Clear();
        firstCopy.AllPeptides.Clear();
        firstCopy.UniquePeptides.Clear();

        Assert.Multiple(() =>
        {
            Assert.That(secondCopy.Proteins.Count, Is.EqualTo(2));
            Assert.That(secondCopy.AllPeptides.Count, Is.EqualTo(2));
            Assert.That(secondCopy.UniquePeptides.Count, Is.EqualTo(2));
        });
    }

    private static object CreateCachedProteinGroup(EngineLayer.ProteinGroup proteinGroup)
    {
        var type = typeof(TaskLayer.ParallelSearch.ParallelSearchTask).Assembly.GetType("TaskLayer.ParallelSearch.CachedProteinGroup", throwOnError: true)!;
        return Activator.CreateInstance(type, proteinGroup)!;
    }

    private static EngineLayer.ProteinGroup CreateRuntimeCopy(object cachedProteinGroup)
    {
        MethodInfo method = cachedProteinGroup.GetType().GetMethod("CreateRuntimeCopy", BindingFlags.Instance | BindingFlags.Public)!;
        return (EngineLayer.ProteinGroup)method.Invoke(cachedProteinGroup, null)!;
    }
}
