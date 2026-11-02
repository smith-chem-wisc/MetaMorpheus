using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer.DatabaseLoading;
using GuiFunctions;
using NUnit.Framework;

namespace Test.GuiTests;

[TestFixture]
public class FastaHeaderParsingViewModelTests
{
    [Test]
    public void TheSavedSettingRoundTripsThroughTheViewModel()
    {
        var saved = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, "acc(1)", "full(2)", "name(3)", "gene(4)", "org(5)", "orgId(6)");

        var roundTripped = new FastaHeaderParsingViewModel(saved).ToParameters();

        Assert.That(roundTripped.HeaderFormat, Is.EqualTo(FastaHeaderFormat.Custom));
        Assert.That(new[]
        {
            roundTripped.CustomAccessionRegex, roundTripped.CustomFullNameRegex, roundTripped.CustomNameRegex,
            roundTripped.CustomGeneNameRegex, roundTripped.CustomOrganismRegex, roundTripped.CustomOrganismIdRegex
        }, Is.EqualTo(new[] { "acc(1)", "full(2)", "name(3)", "gene(4)", "org(5)", "orgId(6)" }));
    }

    [Test]
    public void ATaskWithNoSettingShowsTheDefaultAndEveryFormat()
    {
        var viewModel = new FastaHeaderParsingViewModel(null);

        Assert.That(viewModel.HeaderFormat, Is.EqualTo(FastaHeaderFormat.UniProt));
        Assert.That(viewModel.IsCustom, Is.False);
        Assert.That(viewModel.ToParameters().CustomAccessionRegex, Is.Empty);
        Assert.That(viewModel.AvailableFormats, Is.EqualTo(Enum.GetValues<FastaHeaderFormat>()));
    }

    [Test]
    public void ChoosingCustomEnablesTheRegexBoxesAndClearsTheSummary()
    {
        var viewModel = new FastaHeaderParsingViewModel(new FastaHeaderParsingParameters(FastaHeaderFormat.Gencode));
        var changed = new List<string>();
        viewModel.PropertyChanged += (_, e) => changed.Add(e.PropertyName);
        Assert.That(viewModel.PresetSummary, Is.Not.Empty);

        viewModel.HeaderFormat = FastaHeaderFormat.Custom;

        // The panel binds IsEnabled to IsCustom and the summary text to PresetSummary; both must be told.
        Assert.That(changed, Is.EquivalentTo(new[] { nameof(viewModel.HeaderFormat), nameof(viewModel.IsCustom), nameof(viewModel.PresetSummary) }));
        Assert.That(viewModel.IsCustom, Is.True);
        Assert.That(viewModel.PresetSummary, Is.Empty);
        Assert.That(viewModel.ToParameters().HeaderFormat, Is.EqualTo(FastaHeaderFormat.Custom));
    }

    [Test]
    public void EditingARegexIsSavedAndNotified()
    {
        var viewModel = new FastaHeaderParsingViewModel(new FastaHeaderParsingParameters(FastaHeaderFormat.Custom));
        var changed = new List<string>();
        viewModel.PropertyChanged += (_, e) => changed.Add(e.PropertyName);

        viewModel.CustomAccessionRegex = "a(1)";
        viewModel.CustomFullNameRegex = "b(2)";
        viewModel.CustomNameRegex = "c(3)";
        viewModel.CustomGeneNameRegex = "d(4)";
        viewModel.CustomOrganismRegex = "e(5)";
        viewModel.CustomOrganismIdRegex = "f(6)";

        Assert.That(changed, Is.EqualTo(new[]
        {
            nameof(viewModel.CustomAccessionRegex), nameof(viewModel.CustomFullNameRegex), nameof(viewModel.CustomNameRegex),
            nameof(viewModel.CustomGeneNameRegex), nameof(viewModel.CustomOrganismRegex), nameof(viewModel.CustomOrganismIdRegex)
        }));
        var saved = viewModel.ToParameters();
        Assert.That(new[] { saved.CustomAccessionRegex, saved.CustomFullNameRegex, saved.CustomNameRegex, saved.CustomGeneNameRegex, saved.CustomOrganismRegex, saved.CustomOrganismIdRegex },
            Is.EqualTo(new[] { "a(1)", "b(2)", "c(3)", "d(4)", "e(5)", "f(6)" }));
    }

    [Test]
    public void EveryFormatHasASummaryAndAnUnknownOneIsLoud()
    {
        foreach (var format in Enum.GetValues<FastaHeaderFormat>().Where(f => f != FastaHeaderFormat.Custom))
            Assert.That(FastaHeaderParsingViewModel.DescribePreset(format), Is.Not.Empty, format.ToString());

        Assert.That(() => FastaHeaderParsingViewModel.DescribePreset((FastaHeaderFormat)99), Throws.TypeOf<ArgumentOutOfRangeException>());
    }
}
