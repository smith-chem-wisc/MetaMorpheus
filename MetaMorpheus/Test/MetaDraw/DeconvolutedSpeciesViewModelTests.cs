using System.Collections.Generic;
using System.Linq;
using Chemistry;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot;

namespace Test.MetaDraw;

[TestFixture]
public class DeconvolutedSpeciesViewModelTests
{
    [Test]
    public void Properties_AreCorrect()
    {
        // Arrange: create a fake isotopic envelope
        var peaks = new List<(double mz, double intensity)>
        {
            (500.0, 100.0),
            (501.0, 200.0),
            (502.0, 50.0)
        };
        var envelope = new IsotopicEnvelope(0, peaks, 1000, 2, peaks.Sum(p => p.intensity), 1);

        // Act
        var vm = new DeconvolutedSpeciesViewModel(envelope);

        // Assert
        Assert.That(vm.Envelope, Is.EqualTo(envelope));
        Assert.That(vm.Color, Is.EqualTo(OxyColors.Automatic));
        Assert.That(vm.MonoisotopicMass, Is.EqualTo(1000.0));
        Assert.That(vm.Charge, Is.EqualTo(2));
        Assert.That(vm.PeakCount, Is.EqualTo(3));
        Assert.That(vm.Intensity, Is.EqualTo(350.0));
        Assert.That(vm.Annotation, Does.Contain("M=1000"));
        Assert.That(vm.Annotation, Does.Contain("z=2"));
        Assert.That(vm.MostAbundantMz, Is.EqualTo(501.0));
        Assert.That(vm.PeakMzs, Is.EqualTo("500, 501, 502"));
        Assert.That(vm.MonoMz, Is.EqualTo(1000.0.ToMz(2)));

    }

    [Test]
    public void Properties_CacheValues()
    {
        // Arrange
        var peaks = new List<(double mz, double intensity)>
        {
            (100, 1), (200, 2)
        };
        var envelope = new IsotopicEnvelope(0, peaks, 1000, 3, peaks.Sum(p => p.intensity), 1);
        var vm = new DeconvolutedSpeciesViewModel(envelope);

        // Act
        var intensity1 = vm.Intensity;
        var intensity2 = vm.Intensity;
        var annotation1 = vm.Annotation;
        var annotation2 = vm.Annotation;
        var mostAbundantMz1 = vm.MostAbundantMz;
        var mostAbundantMz2 = vm.MostAbundantMz;
        var peakMzs1 = vm.PeakMzs;
        var peakMzs2 = vm.PeakMzs;

        // Assert: values are cached and consistent
        Assert.That(intensity1, Is.EqualTo(intensity2));
        Assert.That(annotation1, Is.EqualTo(annotation2));
        Assert.That(mostAbundantMz1, Is.EqualTo(mostAbundantMz2));
        Assert.That(peakMzs1, Is.EqualTo(peakMzs2));
    }
}