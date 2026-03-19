using System.IO;
using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using System.Linq;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Transcriptomics;
using Transcriptomics.Digestion;
using GuiFunctions.Util;
using Chemistry;
using GuiFunctions.Models;
using TaskLayer;

namespace Test.GuiTests;

[TestFixture]
public class FragmentParamsVM
{
    [SetUp]
    public void Setup()
    {
        // Clear cache before each test
        CustomMIonLossManager.ClearCache();

        var path = CustomMIonLossManager.GetCustomMIonLossFilePath();
        if (File.Exists(path))
            File.Delete(path);
    }

    [TearDown]
    public void TearDown()
    {
        // Clear cache after each test
        CustomMIonLossManager.ClearCache();
    }

    [Test]
    [NonParallelizable]
    public void LossOptions_ProteinMode()
    {
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var common = new CommonParameters(fragmentationParams: new FragmentationParams(), digestionParams: new DigestionParams());
        var vm = new FragmentationParamsViewModel(common, new());
        Assert.That(vm.AvailableMIonLosses, Is.Not.Empty);
        Assert.That(vm.AvailableMIonLosses.Count, Is.GreaterThanOrEqualTo(1));
    }

    [Test]
    [NonParallelizable]
    public void LossOptions_RNAMode()
    {
        GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
        var common = new CommonParameters(fragmentationParams: new RnaFragmentationParams(), digestionParams: new RnaDigestionParams());
        var vm = new FragmentationParamsViewModel(common, new());
        Assert.That(vm.AvailableMIonLosses, Is.Not.Empty);
        Assert.That(vm.AvailableMIonLosses.Count, Is.GreaterThan(10));
    }

    [Test]
    [NonParallelizable]
    public void AddAndRemoveAllLosses()
    {
        GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
        var vm = new FragmentationParamsViewModel(new(), new());

        vm.AddAllLossesCommand.Execute(null);
        Assert.That(vm.AvailableMIonLosses.Select(p => p.IsSelected), Is.All.EqualTo(true));

        vm.RemoveAllLossesCommand.Execute(null);
        Assert.That(vm.AvailableMIonLosses.Select(p => p.IsSelected), Is.All.EqualTo(false));
    }

    [Test]
    [NonParallelizable]
    public void Constructor_InitializesProperties()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var common = new CommonParameters(
            fragmentationParams: new FragmentationParams { GenerateMIon = true },
            digestionParams: new DigestionParams());
        var searchParams = new SearchParameters
        {
            MaxFragmentSize = 5000,
            MinAllowedInternalFragmentLength = 4
        };

        // Act
        var vm = new FragmentationParamsViewModel(common, searchParams);

        // Assert
        Assert.That(vm.GenerateMIon, Is.True);
        Assert.That(vm.MaxFragmentMassDa, Is.EqualTo(5000));
        Assert.That(vm.MinInternalIonLength, Is.EqualTo(4));
        Assert.That(vm.GenerateInternalIons, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void GenerateInternalIons_SetsMinLength()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());

        // Act
        vm.GenerateInternalIons = true;

        // Assert
        Assert.That(vm.MinInternalIonLength, Is.EqualTo(4)); // Default value
    }

    [Test]
    [NonParallelizable]
    public void GenerateInternalIons_False_ResetsMinLength()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());
        vm.GenerateInternalIons = true;
        vm.MinInternalIonLength = 6;

        // Act
        vm.GenerateInternalIons = false;

        // Assert
        Assert.That(vm.MinInternalIonLength, Is.EqualTo(0));
    }

    [Test]
    [NonParallelizable]
    public void GetSelectedMIonLosses_ReturnsOnlySelected()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());
        
        // Select only the first loss
        vm.AvailableMIonLosses.First().IsSelected = true;

        // Act
        var selectedLosses = vm.GetSelectedMIonLosses();

        // Assert
        Assert.That(selectedLosses.Count, Is.EqualTo(1));
    }

    [Test]
    [NonParallelizable]
    public void ToFragmentationParams_ProteinMode_ReturnsCorrectType()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());
        vm.GenerateMIon = true;

        // Act
        var fragmentationParams = vm.ToFragmentationParams();

        // Assert
        Assert.That(fragmentationParams, Is.TypeOf<FragmentationParams>());
        Assert.That(fragmentationParams.GenerateMIon, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void ToFragmentationParams_RNAMode_ReturnsCorrectType()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
        var vm = new FragmentationParamsViewModel(new(), new());
        vm.GenerateMIon = true;
        vm.ModsCanSuppressBaseLossFragments = true;

        // Act
        var fragmentationParams = vm.ToFragmentationParams();

        // Assert
        Assert.That(fragmentationParams, Is.TypeOf<RnaFragmentationParams>());
        Assert.That(fragmentationParams.GenerateMIon, Is.True);
        Assert.That(((RnaFragmentationParams)fragmentationParams).ModificationsCanSuppressBaseLossIons, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void ToFragmentationParams_IncludesSelectedLosses()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());
        vm.AddAllLossesCommand.Execute(null);

        // Act
        var fragmentationParams = vm.ToFragmentationParams();

        // Assert
        Assert.That(fragmentationParams.MIonLosses, Is.Not.Empty);
        Assert.That(fragmentationParams.MIonLosses.Count, Is.EqualTo(vm.AvailableMIonLosses.Count));
    }

    [Test]
    [NonParallelizable]
    public void ReloadMIonLosses_LoadsCustomLosses()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var formula = ChemicalFormula.ParseFormula("CO2");
        var customLoss = new CustomMIonLoss("Test Loss", "-CO2", formula, AnalyteType.Peptide);
        CustomMIonLossManager.AddCustomMIonLoss(customLoss);
        
        var vm = new FragmentationParamsViewModel(new(), new());
        var initialCount = vm.AvailableMIonLosses.Count;

        // Act
        vm.ReloadMIonLosses();

        // Assert
        Assert.That(vm.AvailableMIonLosses.Count, Is.GreaterThanOrEqualTo(initialCount));
        Assert.That(vm.AvailableMIonLosses.Any(l => l.Name == "Test Loss"), Is.True);
    }

    [Test]
    [NonParallelizable]
    public void ReloadMIonLosses_PreservesSelections()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());
        
        // Select first loss
        var firstLoss = vm.AvailableMIonLosses.First();
        firstLoss.IsSelected = true;
        var firstLossAnnotation = firstLoss.Annotation;

        // Act
        vm.ReloadMIonLosses();

        // Assert
        var reloadedFirstLoss = vm.AvailableMIonLosses.FirstOrDefault(l => l.Annotation == firstLossAnnotation);
        Assert.That(reloadedFirstLoss, Is.Not.Null);
        Assert.That(reloadedFirstLoss.IsSelected, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void LoadsWithPreviouslySelectedLosses()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var preselectedLoss = new MIonLoss("Water Loss", "-H2O", formula);
        var fragmentationParams = new FragmentationParams
        {
            GenerateMIon = true,
            MIonLosses = new System.Collections.Generic.List<MIonLoss> { preselectedLoss }
        };
        var common = new CommonParameters(fragmentationParams: fragmentationParams);

        // Act
        var vm = new FragmentationParamsViewModel(common, new());

        // Assert
        var waterLoss = vm.AvailableMIonLosses.FirstOrDefault(l => l.Annotation == "-H2O");
        Assert.That(waterLoss, Is.Not.Null);
        Assert.That(waterLoss.IsSelected, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void CustomLosses_LoadedBasedOnMode()
    {
        // Arrange
        var peptideFormula = ChemicalFormula.ParseFormula("SO3");
        var oligoFormula = ChemicalFormula.ParseFormula("H1P1O3");
        
        var peptideLoss = new CustomMIonLoss("Peptide Only", "-PEP", peptideFormula, AnalyteType.Peptide);
        var oligoLoss = new CustomMIonLoss("Oligo Only", "-OLI", oligoFormula, AnalyteType.Oligo);
        
        CustomMIonLossManager.AddCustomMIonLoss(peptideLoss);
        CustomMIonLossManager.AddCustomMIonLoss(oligoLoss);

        // Act - Protein Mode
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var proteinVm = new FragmentationParamsViewModel(new(), new());

        // Act - RNA Mode
        GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
        var common = new CommonParameters(fragmentationParams: new RnaFragmentationParams(), digestionParams: new RnaDigestionParams());
        var rnaVm = new FragmentationParamsViewModel(common, new());

        // Assert
        Assert.That(proteinVm.AvailableMIonLosses.Any(l => l.Name == "Peptide Only"), Is.True);
        Assert.That(proteinVm.AvailableMIonLosses.Any(l => l.Name == "Oligo Only"), Is.False);

        Assert.That(rnaVm.AvailableMIonLosses.Any(l => l.Name == "Oligo Only"), Is.True);
        Assert.That(rnaVm.AvailableMIonLosses.Any(l => l.Name == "Peptide Only"), Is.False);
    }

    [Test]
    [NonParallelizable]
    public void PropertyChanged_RaisesEvent()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var vm = new FragmentationParamsViewModel(new(), new());
        bool eventRaised = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.GenerateMIon)) eventRaised = true; };

        // Act
        vm.GenerateMIon = true;

        // Assert
        Assert.That(eventRaised, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void MIonLossViewModel_Properties_ReflectUnderlyingMIonLoss()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H3N1");
        var mIonLoss = new MIonLoss("Ammonia Loss", "-NH3", formula);

        // Act
        var viewModel = new MIonLossViewModel(mIonLoss, isSelected: true, isCustom: true);

        // Assert
        Assert.That(viewModel.Name, Is.EqualTo("Ammonia Loss"));
        Assert.That(viewModel.Annotation, Is.EqualTo("-NH3"));
        Assert.That(viewModel.MassShift, Is.EqualTo(mIonLoss.MonoisotopicMass).Within(0.0001));
        Assert.That(viewModel.IsSelected, Is.True);
        Assert.That(viewModel.IsCustom, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void MIonLossViewModel_IsSelected_CanBeToggled()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var mIonLoss = new MIonLoss("Water Loss", "-H2O", formula);
        var viewModel = new MIonLossViewModel(mIonLoss);

        // Act
        viewModel.IsSelected = true;
        var firstState = viewModel.IsSelected;
        
        viewModel.IsSelected = false;
        var secondState = viewModel.IsSelected;

        // Assert
        Assert.That(firstState, Is.True);
        Assert.That(secondState, Is.False);
    }

    [Test]
    [NonParallelizable]
    public void MIonLossViewModel_PropertyChanged_RaisesEvent()
    {
        // Arrange
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var mIonLoss = new MIonLoss("Water Loss", "-H2O", formula);
        var viewModel = new MIonLossViewModel(mIonLoss);
        bool eventRaised = false;
        viewModel.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(viewModel.IsSelected)) eventRaised = true; };

        // Act
        viewModel.IsSelected = true;

        // Assert
        Assert.That(eventRaised, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void LeftAndRightSideFragmentIons_InitializedFromTerminus()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var digestionParams = new DigestionParams(fragmentationTerminus: FragmentationTerminus.Both);
        var common = new CommonParameters(digestionParams: digestionParams);

        // Act
        var vm = new FragmentationParamsViewModel(common, new());

        // Assert
        Assert.That(vm.LeftSideFragmentIons, Is.True);
        Assert.That(vm.RightSideFragmentIons, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void OnlyNTerminus_SetsCorrectFlags()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var digestionParams = new DigestionParams(fragmentationTerminus: FragmentationTerminus.N);
        var common = new CommonParameters(digestionParams: digestionParams);

        // Act
        var vm = new FragmentationParamsViewModel(common, new());

        // Assert
        Assert.That(vm.LeftSideFragmentIons, Is.True);
        Assert.That(vm.RightSideFragmentIons, Is.False);
    }

    [Test]
    [NonParallelizable]
    public void OnlyCTerminus_SetsCorrectFlags()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var digestionParams = new DigestionParams(fragmentationTerminus: FragmentationTerminus.C);
        var common = new CommonParameters(digestionParams: digestionParams);

        // Act
        var vm = new FragmentationParamsViewModel(common, new());

        // Assert
        Assert.That(vm.LeftSideFragmentIons, Is.False);
        Assert.That(vm.RightSideFragmentIons, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void OnlyFivePrimeTerminus_SetsCorrectFlags()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
        var digestionParams = new DigestionParams(fragmentationTerminus: FragmentationTerminus.FivePrime);
        var common = new CommonParameters(digestionParams: digestionParams);

        // Act
        var vm = new FragmentationParamsViewModel(common, new());

        // Assert
        Assert.That(vm.LeftSideFragmentIons, Is.True);
        Assert.That(vm.RightSideFragmentIons, Is.False);
    }

    [Test]
    [NonParallelizable]
    public void OnlyThreePrimeTerminus_SetsCorrectFlags()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;
        var digestionParams = new DigestionParams(fragmentationTerminus: FragmentationTerminus.ThreePrime);
        var common = new CommonParameters(digestionParams: digestionParams);

        // Act
        var vm = new FragmentationParamsViewModel(common, new());

        // Assert
        Assert.That(vm.LeftSideFragmentIons, Is.False);
        Assert.That(vm.RightSideFragmentIons, Is.True);
    }

    [Test]
    [NonParallelizable]
    public void SettersForTestCoverage()
    {
        // Arrange
        GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

        // Act
        var vm = new FragmentationParamsViewModel(new(), new());
        vm.MaxFragmentMassDa = 6000;
        vm.MinInternalIonLength = 5;
        vm.RightSideFragmentIons = false;
        vm.LeftSideFragmentIons = true;
        vm.GenerateComplementaryIons = true;

        // Assert
        Assert.That(vm.MaxFragmentMassDa, Is.EqualTo(6000));
        Assert.That(vm.MinInternalIonLength, Is.EqualTo(5));
        Assert.That(vm.RightSideFragmentIons, Is.False);
        Assert.That(vm.LeftSideFragmentIons, Is.True);
        Assert.That(vm.GenerateComplementaryIons, Is.True);
    }


}
