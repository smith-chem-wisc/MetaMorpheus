using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using GuiFunctions;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics;

[TestFixture]
public class RnaTerminalCapsTests
{
    private const string ThreePrimeCap = "H2O4P";
    private const string FivePrimeCap = "H1O1";

    [Test]
    public void DefaultsAreNull()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams());
        Assert.That(common.FivePrimeTerminus, Is.Null);
        Assert.That(common.ThreePrimeTerminus, Is.Null);
    }

    [Test]
    public void CapsAreNullInProteinMode()
    {
        var common = new CommonParameters(digestionParams: new DigestionParams(), fivePrimeTerminus: FivePrimeCap, threePrimeTerminus: ThreePrimeCap);
        Assert.That(common.FivePrimeTerminus, Is.Null);
        Assert.That(common.ThreePrimeTerminus, Is.Null);
    }

    [Test]
    public void CapsAreKeptInRnaModeAndSurviveClones()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: FivePrimeCap, threePrimeTerminus: ThreePrimeCap);
        Assert.That(common.FivePrimeTerminus!.Formula, Is.EqualTo(ChemicalFormula.ParseFormula(FivePrimeCap).Formula));
        Assert.That(common.ThreePrimeTerminus!.Formula, Is.EqualTo(ChemicalFormula.ParseFormula(ThreePrimeCap).Formula));

        foreach (var clone in new[] { common.Clone(), common.CloneWithNewTerminus(), common.CloneWithNewTotalPartitions(2), common.CloneWithNewMaxThreadsToUsePerFile(1) })
        {
            Assert.That(clone.FivePrimeTerminus, Is.EqualTo(common.FivePrimeTerminus));
            Assert.That(clone.ThreePrimeTerminus, Is.EqualTo(common.ThreePrimeTerminus));
        }
    }

    [Test]
    public void UnparseableFormulaBecomesNull()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: "not a formula!", threePrimeTerminus: "   ");
        Assert.That(common.FivePrimeTerminus, Is.Null);
        Assert.That(common.ThreePrimeTerminus, Is.Null);
    }

    [Test]
    public void FileSpecificParametersKeepTaskCaps()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: FivePrimeCap, threePrimeTerminus: ThreePrimeCap)
        };
        var fileSpecific = new FileSpecificParameters { MaxMissedCleavages = 1 };
        var result = MetaMorpheusTask.SetAllFileSpecificCommonParams(task.CommonParameters, fileSpecific);
        Assert.That(result.FivePrimeTerminus, Is.EqualTo(task.CommonParameters.FivePrimeTerminus));
        Assert.That(result.ThreePrimeTerminus, Is.EqualTo(task.CommonParameters.ThreePrimeTerminus));
        Assert.That(result.ThreePrimeTerminus, Is.Not.Null);
    }

    [Test]
    public void TomlRoundTrip()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "RnaCapsTask.toml");
        try
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: FivePrimeCap, threePrimeTerminus: ThreePrimeCap)
            };
            Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
            var loaded = Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig);
            Assert.That(loaded.CommonParameters.FivePrimeTerminus, Is.EqualTo(task.CommonParameters.FivePrimeTerminus));
            Assert.That(loaded.CommonParameters.ThreePrimeTerminus, Is.EqualTo(task.CommonParameters.ThreePrimeTerminus));
            Assert.That(loaded.CommonParameters.ThreePrimeTerminus, Is.Not.Null);

            // default (null) caps stay null after a round trip
            var defaults = new SearchTask { CommonParameters = new CommonParameters(digestionParams: new RnaDigestionParams()) };
            Toml.WriteFile(defaults, path, MetaMorpheusTask.tomlConfig);
            var loadedDefaults = Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig);
            Assert.That(loadedDefaults.CommonParameters.FivePrimeTerminus, Is.Null);
            Assert.That(loadedDefaults.CommonParameters.ThreePrimeTerminus, Is.Null);
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    [Test]
    [TestCase("20mer1.fasta", true, true)]
    [TestCase("20mer1.fasta", true, false)]
    [TestCase("20mer1.fasta", false, true)]
    [TestCase("20mer1.xml", true, true)]
    [TestCase("20mer1.xml", true, false)]
    [TestCase("20mer1.xml", false, true)]
    public void DatabaseLoadingAppliesWhicheverCapsAreSet(string databaseFileName, bool setFivePrime, bool setThreePrime)
    {
        GlobalVariables.AnalyteType = AnalyteType.Oligo;
        var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", databaseFileName);

        RNA Load(CommonParameters common)
        {
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };
            var results = (DatabaseLoadingEngineResults)new DatabaseLoadingEngine(common, [], [], dbForTask, "TestTaskId", DecoyType.None).Run()!;
            return (RNA)results.BioPolymers.Single();
        }

        var defaultRna = Load(new CommonParameters(digestionParams: new RnaDigestionParams()));
        var cappedRna = Load(new CommonParameters(digestionParams: new RnaDigestionParams(),
            fivePrimeTerminus: setFivePrime ? FivePrimeCap : null, threePrimeTerminus: setThreePrime ? ThreePrimeCap : null));

        // a cap that was not set keeps whatever the database gave it
        var expectedFive = setFivePrime ? ChemicalFormula.ParseFormula(FivePrimeCap) : defaultRna.FivePrimeTerminus.ThisChemicalFormula;
        var expectedThree = setThreePrime ? ChemicalFormula.ParseFormula(ThreePrimeCap) : defaultRna.ThreePrimeTerminus.ThisChemicalFormula;
        Assert.That(cappedRna.FivePrimeTerminus.ThisChemicalFormula, Is.EqualTo(expectedFive));
        Assert.That(cappedRna.ThreePrimeTerminus.ThisChemicalFormula, Is.EqualTo(expectedThree));

        double expectedShift = expectedFive.MonoisotopicMass + expectedThree.MonoisotopicMass
            - defaultRna.FivePrimeTerminus.ThisChemicalFormula.MonoisotopicMass - defaultRna.ThreePrimeTerminus.ThisChemicalFormula.MonoisotopicMass;
        Assert.That(cappedRna.MonoisotopicMass - defaultRna.MonoisotopicMass, Is.EqualTo(expectedShift).Within(1e-6));
    }

    [Test]
    public void LoadOligoDbWithoutCommonParametersLeavesCapsAtDefault()
    {
        GlobalVariables.AnalyteType = AnalyteType.Oligo;
        var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "20mer1.fasta");
        var rna = DatabaseLoadingEngine.LoadOligoDb(dbPath, true, DecoyType.None, new List<string>(), false, out _, out _, null).Single();
        Assert.That(rna.FivePrimeTerminus.ThisChemicalFormula, Is.EqualTo(NucleicAcid.DefaultFivePrimeTerminus));
        Assert.That(rna.ThreePrimeTerminus.ThisChemicalFormula, Is.EqualTo(NucleicAcid.DefaultThreePrimeTerminus));
    }

    [Test]
    public void TomlWithUnreadableCapFormulaFallsBackToNull()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "RnaCapsBadFormula.toml");
        try
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: FivePrimeCap, threePrimeTerminus: ThreePrimeCap)
            };
            Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
            string text = File.ReadAllText(path);
            string written = $"FivePrimeTerminus = \"{task.CommonParameters.FivePrimeTerminus.Formula}\"";
            Assert.That(text, Does.Contain(written));
            File.WriteAllText(path, text.Replace(written, "FivePrimeTerminus = \"zzz!\""));

            var loaded = Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig);
            Assert.That(loaded.CommonParameters.FivePrimeTerminus, Is.Null);
            Assert.That(loaded.CommonParameters.ThreePrimeTerminus, Is.EqualTo(task.CommonParameters.ThreePrimeTerminus));
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    [Test]
    [TestCase(null, null)]
    [TestCase("", null)]
    [TestCase("   ", null)]
    [TestCase("not a formula!", null)]
    [TestCase("OH", "OH")]
    [TestCase("  H2O4P  ", "H2O4P")]
    public void TryParseFormulaHandlesBlankBadAndPaddedText(string text, string expected)
    {
        var parsed = CommonParameters.TryParseFormula(text);
        if (expected == null)
            Assert.That(parsed, Is.Null);
        else
            Assert.That(parsed, Is.EqualTo(ChemicalFormula.ParseFormula(expected)));
    }

    [Test]
    public void EachCapIsIndependentInRnaMode()
    {
        var fiveOnly = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: FivePrimeCap);
        Assert.That(fiveOnly.FivePrimeTerminus, Is.Not.Null);
        Assert.That(fiveOnly.ThreePrimeTerminus, Is.Null);

        var threeOnly = new CommonParameters(digestionParams: new RnaDigestionParams(), threePrimeTerminus: ThreePrimeCap);
        Assert.That(threeOnly.FivePrimeTerminus, Is.Null);
        Assert.That(threeOnly.ThreePrimeTerminus, Is.Not.Null);
    }

    [Test]
    public void ClonesOfParametersWithoutCapsKeepThemNull()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams());
        foreach (var clone in new[] { common.Clone(), common.CloneWithNewTerminus(), common.CloneWithNewTotalPartitions(2), common.CloneWithNewMaxThreadsToUsePerFile(1) })
        {
            Assert.That(clone.FivePrimeTerminus, Is.Null);
            Assert.That(clone.ThreePrimeTerminus, Is.Null);
        }
    }

    [Test]
    public void FileSpecificParametersWithoutCapsStayNull()
    {
        var fileSpecific = new FileSpecificParameters { MaxMissedCleavages = 1 };
        foreach (var common in new[]
                 {
                     new CommonParameters(digestionParams: new RnaDigestionParams()),
                     new CommonParameters(digestionParams: new DigestionParams())
                 })
        {
            var result = MetaMorpheusTask.SetAllFileSpecificCommonParams(common, fileSpecific);
            Assert.That(result.FivePrimeTerminus, Is.Null);
            Assert.That(result.ThreePrimeTerminus, Is.Null);
        }
    }

    [Test]
    public void CatalogShiftsMatchPublishedValues()
    {
        var vm = new RnaTerminalCapsViewModel();
        double Shift(IEnumerable<TerminalCapOption> options, string name) => options.Single(o => o.Name.StartsWith(name)).MassShift;

        Assert.That(Shift(vm.FivePrimeOptions, "Hydroxyl"), Is.EqualTo(0).Within(1e-9));
        Assert.That(Shift(vm.ThreePrimeOptions, "Hydroxyl"), Is.EqualTo(0).Within(1e-9));

        double hpo3 = ChemicalFormula.ParseFormula("HPO3").MonoisotopicMass; // 79.9663
        Assert.That(Shift(vm.FivePrimeOptions, "Monophosphate"), Is.EqualTo(hpo3).Within(1e-6));
        Assert.That(Shift(vm.FivePrimeOptions, "Diphosphate"), Is.EqualTo(2 * hpo3).Within(1e-6));
        Assert.That(Shift(vm.FivePrimeOptions, "Triphosphate"), Is.EqualTo(3 * hpo3).Within(1e-6));
        Assert.That(Shift(vm.ThreePrimeOptions, "Phosphate"), Is.EqualTo(hpo3).Within(1e-6));

        // published: 5'/3'-phosphate +79.98, 2',3'-cyclic phosphate +61.96
        Assert.That(Shift(vm.FivePrimeOptions, "Monophosphate"), Is.EqualTo(79.98).Within(0.02));
        Assert.That(Shift(vm.ThreePrimeOptions, "Cyclic phosphate"), Is.EqualTo(61.96).Within(0.01));
    }

    [Test]
    public void CatalogDefaultsAreMzLibDefaults()
    {
        var vm = new RnaTerminalCapsViewModel();
        Assert.That(vm.FivePrimeOptions.Single(o => o.IsDefault).CapFormula, Is.EqualTo(NucleicAcid.DefaultFivePrimeTerminus));
        Assert.That(vm.ThreePrimeOptions.Single(o => o.IsDefault).CapFormula, Is.EqualTo(NucleicAcid.DefaultThreePrimeTerminus));
        Assert.That(vm.FivePrimeOptions.Single(o => o.IsDefault).Formula, Is.Null);
    }

    [Test]
    public void ViewModelOffersCatalogNamesThenCustom()
    {
        var vm = new RnaTerminalCapsViewModel();
        Assert.That(vm.FivePrimeOptions.Select(o => o.Name), Is.EqualTo(new[] { "Hydroxyl (default)", "Monophosphate", "Diphosphate", "Triphosphate", "Custom formula..." }));
        Assert.That(vm.ThreePrimeOptions.Select(o => o.Name), Is.EqualTo(new[] { "Hydroxyl (default)", "Phosphate", "Cyclic phosphate", "Custom formula..." }));

        vm.SelectedFivePrimeOption = vm.FivePrimeOptions.Single(o => o.Name == "Triphosphate");
        vm.SelectedThreePrimeOption = vm.ThreePrimeOptions.Single(o => o.Name == "Cyclic phosphate");
        var common = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: vm.FivePrimeTerminus, threePrimeTerminus: vm.ThreePrimeTerminus);
        Assert.That(common.FivePrimeTerminus, Is.EqualTo(ChemicalFormula.ParseFormula("H3O6P2")));
        Assert.That(common.ThreePrimeTerminus, Is.EqualTo(ChemicalFormula.ParseFormula("O3P")));

        // and back again: a loaded task shows the name, not a formula
        var restored = new RnaTerminalCapsViewModel(common);
        Assert.That(restored.SelectedFivePrimeOption.Name, Is.EqualTo("Triphosphate"));
        Assert.That(restored.SelectedThreePrimeOption.Name, Is.EqualTo("Cyclic phosphate"));
        Assert.That(restored.IsFivePrimeCustom, Is.False);
    }

    [Test]
    public void ViewModelDefaultsMapToNull()
    {
        var vm = new RnaTerminalCapsViewModel(new CommonParameters(digestionParams: new RnaDigestionParams()));
        Assert.That(vm.FivePrimeTerminus, Is.Null);
        Assert.That(vm.ThreePrimeTerminus, Is.Null);
        Assert.That(vm.IsFivePrimeCustom, Is.False);
        Assert.That(vm.Validate(out _), Is.True);
    }

    [Test]
    public void ViewModelRestoresPresetAndCustom()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams(), fivePrimeTerminus: FivePrimeCap, threePrimeTerminus: ThreePrimeCap);
        var vm = new RnaTerminalCapsViewModel(common);
        Assert.That(vm.IsFivePrimeCustom, Is.True);
        // the view model round-trips through the formula's canonical text, so compare as formulas
        Assert.That(ChemicalFormula.ParseFormula(vm.FivePrimeTerminus), Is.EqualTo(ChemicalFormula.ParseFormula(FivePrimeCap)));
        Assert.That(vm.IsThreePrimeCustom, Is.False); // H2O4P matches the 3'-phosphate preset
        Assert.That(ChemicalFormula.ParseFormula(vm.ThreePrimeTerminus), Is.EqualTo(ChemicalFormula.ParseFormula(ThreePrimeCap)));
    }

    [Test]
    public void ViewModelRejectsBadCustomFormula()
    {
        var vm = new RnaTerminalCapsViewModel
        {
            SelectedThreePrimeOption = RnaTerminalCapsViewModel.CustomOption,
            ThreePrimeCustomFormula = "zzz!"
        };
        Assert.That(vm.Validate(out string error), Is.False);
        Assert.That(error, Does.Contain("3'"));

        vm.ThreePrimeCustomFormula = "  ";
        Assert.That(vm.ThreePrimeTerminus, Is.Null); // blank custom falls back to the default
        Assert.That(vm.Validate(out _), Is.True);
    }

    [Test]
    public void ViewModelTreatsAnExplicitHydroxylFormulaAsTheDefault()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams(),
            fivePrimeTerminus: NucleicAcid.DefaultFivePrimeTerminus.Formula, threePrimeTerminus: NucleicAcid.DefaultThreePrimeTerminus.Formula);
        var vm = new RnaTerminalCapsViewModel(common);
        Assert.That(vm.SelectedFivePrimeOption.IsDefault, Is.True);
        Assert.That(vm.SelectedThreePrimeOption.IsDefault, Is.True);
        Assert.That(vm.FivePrimeTerminus, Is.Null);
        Assert.That(vm.ThreePrimeTerminus, Is.Null);
    }

    [Test]
    public void ViewModelRestoresACustomThreePrimeCap()
    {
        var common = new CommonParameters(digestionParams: new RnaDigestionParams(), threePrimeTerminus: "C2H3O");
        var vm = new RnaTerminalCapsViewModel(common);
        Assert.That(vm.IsThreePrimeCustom, Is.True);
        Assert.That(vm.IsFivePrimeCustom, Is.False);
        Assert.That(ChemicalFormula.ParseFormula(vm.ThreePrimeCustomFormula), Is.EqualTo(ChemicalFormula.ParseFormula("C2H3O")));
        Assert.That(vm.ThreePrimeTerminus, Is.EqualTo(vm.ThreePrimeCustomFormula));
    }

    [Test]
    public void ViewModelResolvesNoSelectionAndTrimsCustomText()
    {
        var vm = new RnaTerminalCapsViewModel
        {
            SelectedFivePrimeOption = null,
            SelectedThreePrimeOption = RnaTerminalCapsViewModel.CustomOption,
            ThreePrimeCustomFormula = "  H2O4P "
        };
        Assert.That(vm.FivePrimeTerminus, Is.Null);
        Assert.That(vm.IsFivePrimeCustom, Is.False);
        Assert.That(vm.ThreePrimeTerminus, Is.EqualTo("H2O4P"));
        Assert.That(vm.Validate(out string error), Is.True);
        Assert.That(error, Is.Null);
    }

    [Test]
    public void ViewModelRejectsBadCustomFivePrimeFormula()
    {
        var vm = new RnaTerminalCapsViewModel
        {
            SelectedFivePrimeOption = RnaTerminalCapsViewModel.CustomOption,
            FivePrimeCustomFormula = "zzz!"
        };
        Assert.That(vm.Validate(out string error), Is.False);
        Assert.That(error, Does.Contain("5'").And.Contain("zzz!"));
    }

    [Test]
    public void ViewModelRaisesPropertyChangedForEverySetter()
    {
        var vm = new RnaTerminalCapsViewModel();
        var raised = new List<string>();
        vm.PropertyChanged += (_, e) => raised.Add(e.PropertyName);

        vm.SelectedFivePrimeOption = RnaTerminalCapsViewModel.CustomOption;
        vm.SelectedThreePrimeOption = RnaTerminalCapsViewModel.CustomOption;
        vm.FivePrimeCustomFormula = "H";
        vm.ThreePrimeCustomFormula = "OH";

        Assert.That(raised, Is.SupersetOf(new[]
        {
            nameof(RnaTerminalCapsViewModel.SelectedFivePrimeOption), nameof(RnaTerminalCapsViewModel.IsFivePrimeCustom),
            nameof(RnaTerminalCapsViewModel.SelectedThreePrimeOption), nameof(RnaTerminalCapsViewModel.IsThreePrimeCustom),
            nameof(RnaTerminalCapsViewModel.FivePrimeCustomFormula), nameof(RnaTerminalCapsViewModel.ThreePrimeCustomFormula)
        }));
    }

    [Test]
    public void CapOptionsDescribeThemselves()
    {
        var vm = new RnaTerminalCapsViewModel();

        var hydroxyl = vm.FivePrimeOptions.Single(o => o.IsDefault);
        Assert.That(hydroxyl.Name, Is.EqualTo("Hydroxyl (default)"));
        Assert.That(hydroxyl.ToString(), Is.EqualTo(hydroxyl.Name));
        Assert.That(hydroxyl.ToolTip, Does.Contain("Formula:").And.Contain("0.0000 Da"));

        var triphosphate = vm.FivePrimeOptions.Single(o => o.Name == "Triphosphate");
        Assert.That(triphosphate.IsDefault, Is.False);
        Assert.That(triphosphate.Formula, Is.EqualTo(triphosphate.CapFormula.Formula));
        Assert.That(triphosphate.ToolTip, Does.Contain("+239.8990 Da"));

        var custom = RnaTerminalCapsViewModel.CustomOption;
        Assert.That(custom.IsCustom, Is.True);
        Assert.That(custom.Formula, Is.Null);
        Assert.That(custom.MassShift, Is.EqualTo(0));
        Assert.That(custom.ToolTip, Does.Contain("chemical formula"));
    }

    [Test]
    public void ViewModelWithoutParametersOrSelectionFallsBackToDefaults()
    {
        var vm = new RnaTerminalCapsViewModel(null);
        Assert.That(vm.SelectedFivePrimeOption.IsDefault, Is.True);
        Assert.That(vm.SelectedThreePrimeOption.IsDefault, Is.True);

        vm.SelectedThreePrimeOption = null;
        Assert.That(vm.IsThreePrimeCustom, Is.False);
        Assert.That(vm.ThreePrimeTerminus, Is.Null);
    }
}
