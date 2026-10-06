#nullable enable
using System.Collections.Generic;
using System.IO;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.GuiTests.Deconvolution;

[TestFixture]
public class FromFileDeconParamsViewModelTests
{
    private static string Ms1FeatureFixturePath =>
        Path.Combine(TestContext.CurrentContext.TestDirectory,
            "TestData", "file-specific-decon", "TaGe_SA_A549_3_snip_ms1.feature");

    private static string Ms2FeatureFixturePath =>
        Path.Combine(TestContext.CurrentContext.TestDirectory,
            "TestData", "file-specific-decon", "TaGe_SA_A549_3_snip_ms2.feature");

    [Test]
    public void Constructor_SetsParametersAndDeconvolutionType()
    {
        var parameters = new FromFileDeconvolutionParameters("dummy.feature", 1, 12);
        var vm = new FromFileDeconParamsViewModel(parameters);

        Assert.That(vm.Parameters, Is.SameAs(parameters));
        Assert.That(vm.DeconvolutionType, Is.EqualTo(DeconvolutionType.FromFile));
    }

    [Test]
    public void ToString_ReturnsFromFile()
    {
        var vm = new FromFileDeconParamsViewModel(new FromFileDeconvolutionParameters("x.feature", 1, 12));
        Assert.That(vm.ToString(), Is.EqualTo("FromFile"));
    }

    [Test]
    public void FilePath_ReflectsWrappedParameters()
    {
        var parameters = new FromFileDeconvolutionParameters("/abs/path_ms1.feature", 1, 12);
        var vm = new FromFileDeconParamsViewModel(parameters);

        Assert.That(vm.FilePath, Is.EqualTo("/abs/path_ms1.feature"));
    }

    [Test]
    public void FilePath_Set_UpdatesParametersAndRaisesPropertyChanged()
    {
        var parameters = new FromFileDeconvolutionParameters("initial_ms1.feature", 1, 12);
        var vm = new FromFileDeconParamsViewModel(parameters);
        var raised = new List<string?>();
        vm.PropertyChanged += (_, e) => raised.Add(e.PropertyName);

        vm.FilePath = "updated_ms1.feature";

        Assert.That(parameters.FilePath, Is.EqualTo("updated_ms1.feature"));
        Assert.That(raised, Contains.Item(nameof(vm.FilePath)));
    }

    [Test]
    public void FilePath_SetSameValue_DoesNotRaisePropertyChanged()
    {
        var parameters = new FromFileDeconvolutionParameters("same_ms1.feature", 1, 12);
        var vm = new FromFileDeconParamsViewModel(parameters);
        var raised = new List<string?>();
        vm.PropertyChanged += (_, e) => raised.Add(e.PropertyName);

        vm.FilePath = "same_ms1.feature";

        Assert.That(raised, Does.Not.Contain(nameof(vm.FilePath)));
    }

    [Test]
    public void UseGenericScore_ReflectsAndUpdatesParameters()
    {
        var parameters = new FromFileDeconvolutionParameters("x.feature", 1, 12);
        parameters.UseGenericScore = false;
        var vm = new FromFileDeconParamsViewModel(parameters);

        Assert.That(vm.UseGenericScore, Is.False);

        var raised = new List<string?>();
        vm.PropertyChanged += (_, e) => raised.Add(e.PropertyName);
        vm.UseGenericScore = true;

        Assert.That(parameters.UseGenericScore, Is.True);
        Assert.That(raised, Contains.Item(nameof(vm.UseGenericScore)));
    }

    [Test]
    public void UseGenericScore_SetSameValue_DoesNotRaisePropertyChanged()
    {
        var parameters = new FromFileDeconvolutionParameters("x.feature", 1, 12);
        parameters.UseGenericScore = true;
        var vm = new FromFileDeconParamsViewModel(parameters);
        var raised = new List<string?>();
        vm.PropertyChanged += (_, e) => raised.Add(e.PropertyName);

        vm.UseGenericScore = true;

        Assert.That(raised, Does.Not.Contain(nameof(vm.UseGenericScore)));
    }

    [Test]
    public void ChargeState_InheritedFromBaseViewModel()
    {
        var parameters = new FromFileDeconvolutionParameters("x.feature", 2, 18);
        var vm = new FromFileDeconParamsViewModel(parameters);

        Assert.That(vm.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(vm.MaxAssumedChargeState, Is.EqualTo(18));
        Assert.That(vm.Polarity, Is.EqualTo(Polarity.Positive));
    }

    [Test]
    public void Polarity_NegativeMode_InheritedFromBaseViewModel()
    {
        var parameters = new FromFileDeconvolutionParameters("x.feature", -20, -1, Polarity.Negative);
        var vm = new FromFileDeconParamsViewModel(parameters);

        Assert.That(vm.Polarity, Is.EqualTo(Polarity.Negative));
    }

    [Test]
    [TestCase(null)]
    [TestCase("")]
    [TestCase("   ")]
    public void TryValidateFilePath_NullOrWhitespace_ReturnsFalseWithMessage(string? path)
    {
        var result = FromFileDeconParamsViewModel.TryValidateFilePath(path, out var msg);

        Assert.That(result, Is.False);
        Assert.That(msg, Is.Not.Null.And.Not.Empty);
    }

    [Test]
    public void TryValidateFilePath_NonExistentPath_ReturnsFalseWithNotFoundMessage()
    {
        var fakePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "nonexistent_ms1.feature");
        var result = FromFileDeconParamsViewModel.TryValidateFilePath(fakePath, out var msg);

        Assert.That(result, Is.False);
        Assert.That(msg.ToLowerInvariant(), Does.Contain("not found").Or.Contain("not exist"));
    }

    [Test]
    public void TryValidateFilePath_ValidMs1FeatureFile_ReturnsTrue()
    {
        Assume.That(File.Exists(Ms1FeatureFixturePath),
            $"MS1 fixture must exist at: {Ms1FeatureFixturePath}");

        var result = FromFileDeconParamsViewModel.TryValidateFilePath(Ms1FeatureFixturePath, out var msg);

        Assert.That(result, Is.True, $"Validation message: {msg}");
        Assert.That(msg, Is.Empty);
    }

    [Test]
    public void TryValidateFilePath_ValidMs1FeatureFile_DoesNotThrow()
    {
        Assume.That(File.Exists(Ms1FeatureFixturePath),
            $"MS1 fixture must exist at: {Ms1FeatureFixturePath}");

        Assert.DoesNotThrow(() =>
            FromFileDeconParamsViewModel.TryValidateFilePath(Ms1FeatureFixturePath, out _));

        var vm = new FromFileDeconParamsViewModel(
            new FromFileDeconvolutionParameters(Ms1FeatureFixturePath, 1, 12));

        Assert.That((vm.Parameters as FromFileDeconvolutionParameters)!.FilePath,
            Is.EqualTo(Ms1FeatureFixturePath));
    }

    [Test]
    public void TryValidateFilePath_Ms2FeatureFile_ReturnsFalseAsNonMs1()
    {
        Assume.That(File.Exists(Ms2FeatureFixturePath),
            $"MS2 fixture must exist at: {Ms2FeatureFixturePath}");

        var result = FromFileDeconParamsViewModel.TryValidateFilePath(Ms2FeatureFixturePath, out var msg);

        Assert.That(result, Is.False, "Non-MS1 result file must be rejected");
        Assert.That(msg, Is.Not.Null.And.Not.Empty);
    }

    [Test]
    public void TryValidateFilePath_RelativePath_ReturnsFalseWithAbsolutePathMessage()
    {
        var result = FromFileDeconParamsViewModel.TryValidateFilePath(
            "relative/path_ms1.feature", out var msg);

        Assert.That(result, Is.False);
        Assert.That(msg.ToLowerInvariant(), Does.Contain("absolute"),
            "Validation message must explain that an absolute path is required");
    }

    [Test]
    public void TryValidateFilePath_DriveRelativePath_ReturnsFalseWithAbsolutePathMessage()
    {
        if (!System.OperatingSystem.IsWindows())
            Assert.Ignore("Drive-relative paths are Windows-specific");

        var result = FromFileDeconParamsViewModel.TryValidateFilePath("C:relative.feature", out var msg);

        Assert.That(result, Is.False,
            "A drive-relative path such as 'C:relative.feature' is not fully qualified and must be rejected");
        Assert.That(msg.ToLowerInvariant(), Does.Contain("fully qualified").Or.Contain("absolute"),
            "Validation message must explain that a fully qualified absolute path is required");
    }
}
