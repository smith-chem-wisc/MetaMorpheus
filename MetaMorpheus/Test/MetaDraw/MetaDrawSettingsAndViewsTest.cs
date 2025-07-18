using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Media;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using OxyPlot;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using System.Collections.ObjectModel;
using System.Threading.Tasks;
using GuiFunctions.MetaDraw;

namespace Test.MetaDraw
{
    public static class MetaDrawSettingsAndViewsTest
    {
        [SetUp]
        public static void SetUp()
        {
            MetaDrawSettings.ResetSettings();
        }

        [Test]
        public static void TestMetaDrawSettingsSnapshot()
        {
            MetaDrawSettingsSnapshot snapshot = new();
            Assert.That(snapshot.DisplayIonAnnotations.Equals(MetaDrawSettings.DisplayIonAnnotations));
            Assert.That(snapshot.AnnotateMzValues.Equals(MetaDrawSettings.AnnotateMzValues));
            Assert.That(snapshot.AnnotateCharges.Equals(MetaDrawSettings.AnnotateCharges));
            Assert.That(snapshot.AnnotationBold.Equals(MetaDrawSettings.AnnotationBold));
            Assert.That(snapshot.ShowDecoys.Equals(MetaDrawSettings.ShowDecoys));
            Assert.That(snapshot.ShowContaminants.Equals(MetaDrawSettings.ShowContaminants));
            Assert.That(snapshot.QValueFilter.Equals(MetaDrawSettings.QValueFilter));
            Assert.That(snapshot.LocalizationLevelStart.Equals(MetaDrawSettings.LocalizationLevelStart));
            Assert.That(snapshot.LocalizationLevelEnd.Equals(MetaDrawSettings.LocalizationLevelEnd));
            Assert.That(snapshot.DisplayInternalIons, Is.EqualTo(MetaDrawSettings.DisplayInternalIons));
            Assert.That(snapshot.SubAndSuperScriptIons, Is.EqualTo(MetaDrawSettings.SubAndSuperScriptIons));
            Assert.That(snapshot.StrokeThicknessAnnotated, Is.EqualTo(MetaDrawSettings.StrokeThicknessAnnotated));
            Assert.That(snapshot.StrokeThicknessUnannotated, Is.EqualTo(MetaDrawSettings.StrokeThicknessUnannotated));
            Assert.That(snapshot.AnnotatedFontSize, Is.EqualTo(MetaDrawSettings.AnnotatedFontSize));
            Assert.That(snapshot.AxisLabelTextSize, Is.EqualTo(MetaDrawSettings.AxisLabelTextSize));
            Assert.That(snapshot.AxisTitleTextSize, Is.EqualTo(MetaDrawSettings.AxisTitleTextSize));
            Assert.That(snapshot.DisplayChimeraLegend, Is.EqualTo(MetaDrawSettings.DisplayChimeraLegend));
            Assert.That(snapshot.ChimeraLegendMainTextType, Is.EqualTo(MetaDrawSettings.ChimeraLegendMainTextType));
            Assert.That(snapshot.ChimeraLegendSubTextType, Is.EqualTo(MetaDrawSettings.ChimeraLegendSubTextType));
            Assert.That(snapshot.SuppressMessageBoxes, Is.EqualTo(MetaDrawSettings.SuppressMessageBoxes));
            Assert.That(snapshot.ChimeraLegendTakeFirstIfAmbiguous, Is.EqualTo(MetaDrawSettings.ChimeraLegendTakeFirstIfAmbiguous));


            MetaDrawSettings.ShowContaminants = true;
            MetaDrawSettings.AnnotateMzValues = false;
            snapshot = MetaDrawSettings.MakeSnapShot();
            Assert.That(snapshot.ShowContaminants.Equals(MetaDrawSettings.ShowContaminants));
            Assert.That(snapshot.AnnotateMzValues.Equals(MetaDrawSettings.AnnotateMzValues));
            Assert.That(snapshot.QValueFilter.Equals(MetaDrawSettings.QValueFilter));
            Assert.That(snapshot.LocalizationLevelStart.Equals(MetaDrawSettings.LocalizationLevelStart));
            Assert.That(snapshot.ExportType.Equals(MetaDrawSettings.ExportType));
            Assert.That(snapshot.DisplayInternalIons, Is.EqualTo(MetaDrawSettings.DisplayInternalIons));
            Assert.That(snapshot.SubAndSuperScriptIons, Is.EqualTo(MetaDrawSettings.SubAndSuperScriptIons));
            Assert.That(snapshot.SuppressMessageBoxes.Equals(MetaDrawSettings.SuppressMessageBoxes));

            var colorValues = MetaDrawSettings.ProductTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            var betaColorValues = MetaDrawSettings.BetaProductTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            var modificationColorValues = MetaDrawSettings.ModificationTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            var coverageColorValues = MetaDrawSettings.CoverageTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            var spectrumDescriptionValues = MetaDrawSettings.SpectrumDescription
                .Select(p => $"{p.Key},{p.Value}").ToList();


            Assert.That(!snapshot.ProductTypeToColorValues.Except(colorValues).Any());
            Assert.That(!snapshot.BetaProductTypeToColorValues.Except(betaColorValues).Any());
            Assert.That(!snapshot.ModificationTypeToColorValues.Except(modificationColorValues).Any());
            Assert.That(!snapshot.CoverageTypeToColorValues.Except(coverageColorValues).Any());
            Assert.That(!snapshot.SpectrumDescriptionValues.Except(spectrumDescriptionValues).Any());

            snapshot.QValueFilter = 0.5;
            snapshot.AnnotateCharges = true;
            snapshot.DisplayIonAnnotations = true;
            snapshot.AnnotateMzValues = true;
            snapshot.AnnotationBold = true;
            snapshot.ShowContaminants = false;
            snapshot.ShowDecoys = true;
            snapshot.ShowLegend = false;
            snapshot.DrawNumbersUnderStationary = false;
            snapshot.SubAndSuperScriptIons = false;
            snapshot.DisplayChimeraLegend = false;
            snapshot.ChimeraLegendMainTextType = LegendDisplayProperty.ProteinAccession;
            snapshot.ChimeraLegendSubTextType = LegendDisplayProperty.FullSequence;
            snapshot.SuppressMessageBoxes = true;
            snapshot.ChimeraLegendTakeFirstIfAmbiguous = true;

            MetaDrawSettings.LoadSettings(snapshot, out bool flaggedError);
            Assert.That(!flaggedError);
            Assert.That(snapshot.DisplayIonAnnotations.Equals(MetaDrawSettings.DisplayIonAnnotations));
            Assert.That(snapshot.AnnotateMzValues.Equals(MetaDrawSettings.AnnotateMzValues));
            Assert.That(snapshot.AnnotateCharges.Equals(MetaDrawSettings.AnnotateCharges));
            Assert.That(snapshot.AnnotationBold.Equals(MetaDrawSettings.AnnotationBold));
            Assert.That(snapshot.ShowDecoys.Equals(MetaDrawSettings.ShowDecoys));
            Assert.That(snapshot.ShowContaminants.Equals(MetaDrawSettings.ShowContaminants));
            Assert.That(snapshot.QValueFilter.Equals(MetaDrawSettings.QValueFilter));
            Assert.That(snapshot.LocalizationLevelStart.Equals(MetaDrawSettings.LocalizationLevelStart));
            Assert.That(snapshot.LocalizationLevelEnd.Equals(MetaDrawSettings.LocalizationLevelEnd));
            Assert.That(snapshot.DisplayInternalIons, Is.EqualTo(MetaDrawSettings.DisplayInternalIons));
            Assert.That(snapshot.SubAndSuperScriptIons, Is.EqualTo(MetaDrawSettings.SubAndSuperScriptIons));
            Assert.That(snapshot.DisplayChimeraLegend, Is.EqualTo(MetaDrawSettings.DisplayChimeraLegend));
            Assert.That(snapshot.ChimeraLegendMainTextType, Is.EqualTo(MetaDrawSettings.ChimeraLegendMainTextType));
            Assert.That(snapshot.ChimeraLegendSubTextType, Is.EqualTo(MetaDrawSettings.ChimeraLegendSubTextType));
            Assert.That(snapshot.SuppressMessageBoxes, Is.EqualTo(MetaDrawSettings.SuppressMessageBoxes));
            Assert.That(snapshot.ChimeraLegendTakeFirstIfAmbiguous, Is.EqualTo(MetaDrawSettings.ChimeraLegendTakeFirstIfAmbiguous));
            colorValues = MetaDrawSettings.ProductTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            betaColorValues = MetaDrawSettings.BetaProductTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            modificationColorValues = MetaDrawSettings.ModificationTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            coverageColorValues = MetaDrawSettings.CoverageTypeToColor
                .Select(p => $"{p.Key},{p.Value.GetColorName()}").ToList();
            spectrumDescriptionValues = MetaDrawSettings.SpectrumDescription
                .Select(p => $"{p.Key},{p.Value}").ToList();
            Assert.That(!snapshot.ProductTypeToColorValues.Except(colorValues).Any());
            Assert.That(!snapshot.BetaProductTypeToColorValues.Except(betaColorValues).Any());
            Assert.That(!snapshot.ModificationTypeToColorValues.Except(modificationColorValues).Any());
            Assert.That(!snapshot.CoverageTypeToColorValues.Except(coverageColorValues).Any());
            Assert.That(!snapshot.SpectrumDescriptionValues.Except(spectrumDescriptionValues).Any());

            snapshot.AnnotatedFontSize = 0;
            snapshot.AxisLabelTextSize = 0;
            snapshot.AxisTitleTextSize = 0;
            snapshot.StrokeThicknessAnnotated = 0;
            snapshot.StrokeThicknessUnannotated = 0;
            MetaDrawSettings.LoadSettings(snapshot, out flaggedError);
            Assert.That(!flaggedError);
            Assert.That(MetaDrawSettings.AnnotatedFontSize, Is.EqualTo(14));
            Assert.That(MetaDrawSettings.AxisLabelTextSize, Is.EqualTo(12));
            Assert.That(MetaDrawSettings.AxisTitleTextSize, Is.EqualTo(14));
            Assert.That(MetaDrawSettings.StrokeThicknessAnnotated, Is.EqualTo(1));
            Assert.That(MetaDrawSettings.StrokeThicknessUnannotated, Is.EqualTo(0.7));
        }

        [Test]
        public static void TestSaveAndLoadDefaultSettings()
        {
            MetaDrawSettingsViewModel model = new MetaDrawSettingsViewModel(false);

            string outputFolder =
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawWithSpectraLibrary");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Assert.That(!Directory.Exists(outputFolder));

            MetaDrawSettingsViewModel.SettingsPath = Path.Combine(outputFolder, @"MetaDrawSettingsDefault.xml");
            Assert.That(model.HasDefaultSaved == false);
            try
            {
                model.LoadSettings();
                Assert.Fail();
            }
            catch (Exception)
            {
                // ignored
            }

            var toTest = model.Modifications.SelectMany(p => p.Children)
                .First(p => p.ModName == "Carbamidomethyl on C");
            Assert.That(toTest.SelectedColor, Is.Not.EqualTo("Blue"));
            toTest.SelectionChanged("Blue");
            Assert.That(toTest.SelectedColor, Is.EqualTo("Blue"));
            model.SaveAsDefault();
            Assert.That(model.HasDefaultSaved == true);
            model.LoadSettings();

            MetaDrawSettingsViewModel model2 = new();
            Assert.That(model2.Modifications.First().Children.First().SelectedColor, Is.EqualTo("Blue"));

            MetaDrawSettingsViewModel model3 = new(false);
            Assert.That(model3.Modifications.First().Children.First().SelectedColor, Is.EqualTo("Blue"));

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestSettingsViewLoading()
        {
            MetaDrawSettingsViewModel BlankSettingsView = MetaDrawSettingsViewModel.Instance;
            BlankSettingsView.Modifications = new();
            BlankSettingsView.IonGroups = new();
            BlankSettingsView.CoverageColors = new();

            Assert.That(BlankSettingsView.IonGroups.Count == 0);
            Assert.That(BlankSettingsView.CoverageColors.Count == 0);
            Assert.That(BlankSettingsView.Modifications.Count == 0);
            Assert.That(!BlankSettingsView.CanOpen);

            BlankSettingsView.LoadIonTypes();
            Assert.That(BlankSettingsView.IonGroups.Count > 0);
            Assert.That(!BlankSettingsView.CanOpen);
            BlankSettingsView.LoadSequenceCoverage();
            Assert.That(BlankSettingsView.CoverageColors.Count > 0);
            Assert.That(!BlankSettingsView.CanOpen);
            BlankSettingsView.LoadPTMs();
            Assert.That(BlankSettingsView.Modifications.Count > 0);
            Assert.That(BlankSettingsView.CanOpen);
        }


        [Test]
        public static void TestLoadMetaDrawSettings()
        {
            // Set up testing environment
            var metaDrawTestingDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw");
            var testingDir = Path.Combine(metaDrawTestingDirectory, "TempTestData");
            if (Directory.Exists(testingDir))
                Directory.Delete(testingDir, true);
            Directory.CreateDirectory(testingDir);

            var masterSettings = Path.Combine(metaDrawTestingDirectory, @"105MetaDrawSettingsSaved.xml");
            string outdatedSettingsPath = Path.Combine(testingDir, @"105MetaDrawSettingsSaved_COPY.xml");
            File.Copy(masterSettings, outdatedSettingsPath);

            // Save default values currently stored in MetaDrawSettings
            var defaultCoverageColors = MetaDrawSettings.CoverageTypeToColor.Values.ToList();
            var defaultColorValues = MetaDrawSettings.ProductTypeToColor.Values.ToList();

            // CASE: user does not have a settings config file stored
            // Desired outcome: All default values are used
            MetaDrawSettingsViewModel viewModel = new MetaDrawSettingsViewModel(false);
            var modelSettingsPath = MetaDrawSettingsViewModel.SettingsPath;
            Assert.That(!File.Exists(modelSettingsPath));
            Assert.That(!viewModel.HasDefaultSaved);

            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(0), Is.EqualTo(defaultCoverageColors[0]));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(1), Is.EqualTo(defaultCoverageColors[1]));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(2), Is.EqualTo(defaultCoverageColors[2]));

            // Make one change and save this viewModel as a modern settings config for next test case
            viewModel.CoverageColors.First().SelectionChanged("Yellow");
            string modernSettingsPath = Path.Combine(testingDir, "temporaryCurrentSettingsConfig.xml");
            MetaDrawSettingsViewModel.SettingsPath = modernSettingsPath;
            viewModel.SaveAsDefault();

            // CASE: user has modern settings config file stored
            // Desired outcome: Read in correctly 
            modelSettingsPath = MetaDrawSettingsViewModel.SettingsPath;
            Assert.That(File.Exists(modelSettingsPath));
            var creationTime = File.GetLastWriteTime(modernSettingsPath);
            viewModel = new MetaDrawSettingsViewModel(false);
            // check that no new file was generated
            var creationTimeAfterLoad = File.GetLastWriteTime(modernSettingsPath);
            Assert.That(creationTimeAfterLoad, Is.EqualTo(creationTime));

            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(0), !Is.EqualTo(defaultCoverageColors[0]));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(0), Is.EqualTo(OxyColors.Yellow));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(1), Is.EqualTo(defaultCoverageColors[1]));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(2), Is.EqualTo(defaultCoverageColors[2]));



            // CASE: user has old settings config file stored
            // Desired outcome: We use what we can salvage, default the rest, delete old file, replace with modern
            MetaDrawSettingsViewModel.SettingsPath = outdatedSettingsPath;

            // Ensure settings are in outdated format
            int commaCount = File.ReadLines(outdatedSettingsPath).Sum(p => p.Count(m => m.Equals(',')));
            Assert.That(commaCount, Is.EqualTo(0));

            creationTime = File.GetLastWriteTime(outdatedSettingsPath);
            viewModel = new MetaDrawSettingsViewModel(false);
            creationTimeAfterLoad = File.GetLastWriteTime(outdatedSettingsPath);
            Assert.That(creationTimeAfterLoad, !Is.EqualTo(creationTime));

            Assert.That(MetaDrawSettings.ProductTypeToColor.Values.ElementAt(0), Is.EqualTo(defaultColorValues[0]));
            Assert.That(MetaDrawSettings.ProductTypeToColor.Values.ElementAt(1), Is.EqualTo(defaultColorValues[1]));
            Assert.That(MetaDrawSettings.ProductTypeToColor.Values.ElementAt(2), Is.EqualTo(defaultColorValues[2]));

            Directory.Delete(testingDir, true);
        }

        [Test]
        public static void TestOldMetaDrawSettingsFileDoesNotCrash()
        {
            // Save default values currently stored in MetaDrawSettings
            var defaultColorValues = MetaDrawSettings.ProductTypeToColor.Values.ToList();
            var defaultBetaColorValues = MetaDrawSettings.BetaProductTypeToColor.Values.ToList();
            var defaultModificationColorValues = MetaDrawSettings.ModificationTypeToColor.Values.ToList();
            var defaultCoverageColors = MetaDrawSettings.CoverageTypeToColor.Values.ToList();
            var defaultSpectrumDescriptionValues = MetaDrawSettings.SpectrumDescription.Values.ToList();

            // Load in an outdated settings file and ensure no crashes occur
            string metaDrawSettingsPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", @"105MetaDrawSettingsSaved.xml");
            MetaDrawSettingsViewModel model = new MetaDrawSettingsViewModel(false);
            Assert.That(!model.HasDefaultSaved);

            MetaDrawSettingsViewModel.SettingsPath = metaDrawSettingsPath;
            Assert.That(model.HasDefaultSaved);

            model.LoadSettings();

            // In this case (6/27/24), the product type, beta product type, and spectrum descriptors will fail
            // As of (7/11/24), new modifications were added and now the modification type will fail when loading into MetaDrawSettings
                // This is okay and working as intended. 

            // check that failed settings loaded to default
            CollectionAssert.AreEqual(defaultColorValues, MetaDrawSettings.ProductTypeToColor.Values);
            CollectionAssert.AreEqual(defaultBetaColorValues, MetaDrawSettings.BetaProductTypeToColor.Values);
            CollectionAssert.AreEqual(defaultSpectrumDescriptionValues, MetaDrawSettings.SpectrumDescription.Values);
            CollectionAssert.AreEqual(defaultModificationColorValues, MetaDrawSettings.ModificationTypeToColor.Values);

            // check successful settings loaded correctly, in this case they were set to aqua
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(0), Is.EqualTo(OxyColors.Aqua));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(0), Is.Not.EqualTo(defaultCoverageColors[0]));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(1), Is.EqualTo(OxyColors.Aqua));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(1), Is.Not.EqualTo(defaultCoverageColors[1]));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(2), Is.EqualTo(OxyColors.Aqua));
            Assert.That(MetaDrawSettings.CoverageTypeToColor.Values.ElementAt(2), Is.Not.EqualTo(defaultCoverageColors[2]));
        }

        [Test] // This test passes by not crashing
        public static void TestMetaDrawSettingsLoadSettingsCases()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", @"105MetaDrawSettingsSavedEditedForTestCoverageFailures.xml");
            var snapShot = XmlReaderWriter.ReadFromXmlFile<MetaDrawSettingsSnapshot>(path);
            MetaDrawSettings.LoadSettings(snapShot, out bool flaggedError);
            Assert.That(flaggedError);

            path = Path.Combine(TestContext.CurrentContext.TestDirectory, "MetaDraw", @"105MetaDrawSettingsSavedEditedForTestCoverageSuccess.xml");
            snapShot = XmlReaderWriter.ReadFromXmlFile<MetaDrawSettingsSnapshot>(path);
            MetaDrawSettings.LoadSettings(snapShot, out flaggedError);
            Assert.That(flaggedError);
        }

        [Test]
        public static void TestSettingsViewSaveAndChildSelectionChanged()
        {
            MetaDrawSettingsViewModel view = new MetaDrawSettingsViewModel(false);
            Assert.That(!view.IonGroups.First().Ions.First().HasChanged);
            view.IonGroups.First().Ions.First().SelectionChanged("Blue");
            Assert.That(view.IonGroups.First().Ions.First().HasChanged);
            Assert.That(view.IonGroups.First().Ions.First().SelectedColor == "Blue");
            Assert.That(view.IonGroups.First().Ions.First().ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromName("Blue").Color);

            Assert.That(!view.Modifications.First().Children.First().HasChanged);
            view.Modifications.First().Children.First().SelectionChanged("Blue");
            Assert.That(view.Modifications.First().Children.First().HasChanged);
            Assert.That(view.Modifications.First().Children.First().SelectedColor == "Blue");
            Assert.That(view.Modifications.First().Children.First().ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromName("Blue").Color);

            Assert.That(!view.CoverageColors.First().HasChanged);
            view.CoverageColors.First().SelectionChanged("Blue");
            Assert.That(view.CoverageColors.First().HasChanged);
            Assert.That(view.CoverageColors.First().SelectedColor == "Blue");
            Assert.That(view.CoverageColors.First().ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromName("Blue").Color);

            var internalIonIonTypeForTreeView = view.IonGroups.First().Ions.First(p => p.IonName == "Internal Ion");
            Assert.That(!internalIonIonTypeForTreeView.HasChanged);
            internalIonIonTypeForTreeView.SelectionChanged("Blue");
            Assert.That(internalIonIonTypeForTreeView.HasChanged);
            Assert.That(internalIonIonTypeForTreeView.SelectedColor == "Blue");
            Assert.That(internalIonIonTypeForTreeView.ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromName("Blue").Color);

            internalIonIonTypeForTreeView = view.IonGroups.First().Ions.First(p => p.IonName == "Unannotated Peak");
            Assert.That(!internalIonIonTypeForTreeView.HasChanged);
            internalIonIonTypeForTreeView.SelectionChanged("Blue");
            Assert.That(internalIonIonTypeForTreeView.HasChanged);
            Assert.That(internalIonIonTypeForTreeView.SelectedColor == "Blue");
            Assert.That(internalIonIonTypeForTreeView.ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromName("Blue").Color);

            view.Save();
            Assert.That(MetaDrawSettings.ProductTypeToColor[view.IonGroups.First().Ions.First().IonType] ==
                        OxyColors.Blue);
            Assert.That(MetaDrawSettings.ModificationTypeToColor[view.Modifications.First().Children.First().ModName] ==
                        OxyColors.Blue);
            Assert.That(MetaDrawSettings.CoverageTypeToColor[view.CoverageColors.First().Name] == OxyColors.Blue);
            Assert.That(MetaDrawSettings.InternalIonColor == OxyColors.Blue);
            Assert.That(MetaDrawSettings.UnannotatedPeakColor == OxyColors.Blue);
        }

        [Test]
        public static void TestCoverageTypeForTreeView()
        {
            CoverageTypeForTreeViewModel coverageTypeForTreeView = new("N-Terminal Color");
            Assert.That(coverageTypeForTreeView.Name == "N-Terminal Color");
            var color = MetaDrawSettings.CoverageTypeToColor["N-Terminal Color"];
            Assert.That(coverageTypeForTreeView.SelectedColor == color.GetColorName());
            Assert.That(coverageTypeForTreeView.ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromOxyColor(color).Color);
        }

        [Test]
        public static void TestMetaDrawSettingsViewModelPublicProperties()
        {
            // Arrange
            var viewModel = new MetaDrawSettingsViewModel(false);

            // Test Modifications property
            Assert.That(viewModel.Modifications, Is.Not.Null);
            Assert.That(viewModel.Modifications, Is.InstanceOf<ObservableCollection<ModTypeForTreeViewModel>>());

            // Test IonGroups property
            Assert.That(viewModel.IonGroups, Is.Not.Null);
            Assert.That(viewModel.IonGroups, Is.InstanceOf<ObservableCollection<IonTypeForTreeViewModel>>());

            // Test CoverageColors property
            Assert.That(viewModel.CoverageColors, Is.Not.Null);
            Assert.That(viewModel.CoverageColors, Is.InstanceOf<ObservableCollection<CoverageTypeForTreeViewModel>>());

            // Test SpectrumDescriptors property
            Assert.That(viewModel.SpectrumDescriptors, Is.Not.Null);
            Assert.That(viewModel.SpectrumDescriptors, Is.InstanceOf<ObservableCollection<SpectrumDescriptorViewModel>>());

            // Test PossibleColors property
            Assert.That(viewModel.PossibleColors, Is.Not.Null);
            Assert.That(viewModel.PossibleColors, Is.InstanceOf<ObservableCollection<string>>());

            // Test HasDefaultSaved property
            Assert.That(viewModel.HasDefaultSaved, Is.TypeOf<bool>());

            // Test CanOpen property
            Assert.That(viewModel.CanOpen, Is.TypeOf<bool>());

            // Test Initialization property
            Assert.That(viewModel.Initialization, Is.Not.Null);
            Assert.That(viewModel.Initialization, Is.InstanceOf<Task>());

            // Test static SettingsPath property
            Assert.That(MetaDrawSettingsViewModel.SettingsPath, Is.Not.Null.And.Not.Empty);

            // Test DisplayIonAnnotations property
            bool originalDisplayIonAnnotations = viewModel.DisplayIonAnnotations;
            viewModel.DisplayIonAnnotations = !originalDisplayIonAnnotations;
            Assert.That(viewModel.DisplayIonAnnotations, Is.EqualTo(!originalDisplayIonAnnotations));
            viewModel.DisplayIonAnnotations = originalDisplayIonAnnotations;

            // Test AnnotateMzValues property
            bool originalAnnotateMzValues = viewModel.AnnotateMzValues;
            viewModel.AnnotateMzValues = !originalAnnotateMzValues;
            Assert.That(viewModel.AnnotateMzValues, Is.EqualTo(!originalAnnotateMzValues));
            viewModel.AnnotateMzValues = originalAnnotateMzValues;

            // Test AnnotateCharges property
            bool originalAnnotateCharges = viewModel.AnnotateCharges;
            viewModel.AnnotateCharges = !originalAnnotateCharges;
            Assert.That(viewModel.AnnotateCharges, Is.EqualTo(!originalAnnotateCharges));
            viewModel.AnnotateCharges = originalAnnotateCharges;

            // Test DisplayInternalIonAnnotations property
            bool originalDisplayInternalIonAnnotations = viewModel.DisplayInternalIonAnnotations;
            viewModel.DisplayInternalIonAnnotations = !originalDisplayInternalIonAnnotations;
            Assert.That(viewModel.DisplayInternalIonAnnotations, Is.EqualTo(!originalDisplayInternalIonAnnotations));
            viewModel.DisplayInternalIonAnnotations = originalDisplayInternalIonAnnotations;

            // Test DisplayInternalIons property
            bool originalDisplayInternalIons = viewModel.DisplayInternalIons;
            viewModel.DisplayInternalIons = !originalDisplayInternalIons;
            Assert.That(viewModel.DisplayInternalIons, Is.EqualTo(!originalDisplayInternalIons));
            viewModel.DisplayInternalIons = originalDisplayInternalIons;

            // Test AnnotationBold property
            bool originalAnnotationBold = viewModel.AnnotationBold;
            viewModel.AnnotationBold = !originalAnnotationBold;
            Assert.That(viewModel.AnnotationBold, Is.EqualTo(!originalAnnotationBold));
            viewModel.AnnotationBold = originalAnnotationBold;

            // Test SubAndSuperScriptIons property
            bool originalSubAndSuperScriptIons = viewModel.SubAndSuperScriptIons;
            viewModel.SubAndSuperScriptIons = !originalSubAndSuperScriptIons;
            Assert.That(viewModel.SubAndSuperScriptIons, Is.EqualTo(!originalSubAndSuperScriptIons));
            viewModel.SubAndSuperScriptIons = originalSubAndSuperScriptIons;

            // Test AnnotatedFontSize property
            int originalAnnotatedFontSize = viewModel.AnnotatedFontSize;
            viewModel.AnnotatedFontSize = originalAnnotatedFontSize + 1;
            Assert.That(viewModel.AnnotatedFontSize, Is.EqualTo(originalAnnotatedFontSize + 1));
            viewModel.AnnotatedFontSize = originalAnnotatedFontSize;

            // Test AxisLabelTextSize property
            int originalAxisLabelTextSize = viewModel.AxisLabelTextSize;
            viewModel.AxisLabelTextSize = originalAxisLabelTextSize + 1;
            Assert.That(viewModel.AxisLabelTextSize, Is.EqualTo(originalAxisLabelTextSize + 1));
            viewModel.AxisLabelTextSize = originalAxisLabelTextSize;

            // Test AxisTitleTextSize property
            int originalAxisTitleTextSize = viewModel.AxisTitleTextSize;
            viewModel.AxisTitleTextSize = originalAxisTitleTextSize + 1;
            Assert.That(viewModel.AxisTitleTextSize, Is.EqualTo(originalAxisTitleTextSize + 1));
            viewModel.AxisTitleTextSize = originalAxisTitleTextSize;

            // Test StrokeThicknessAnnotated property
            double originalStrokeThicknessAnnotated = viewModel.StrokeThicknessAnnotated;
            viewModel.StrokeThicknessAnnotated = originalStrokeThicknessAnnotated + 0.1;
            Assert.That(viewModel.StrokeThicknessAnnotated, Is.EqualTo(originalStrokeThicknessAnnotated + 0.1));
            viewModel.StrokeThicknessAnnotated = originalStrokeThicknessAnnotated;

            // Test StrokeThicknessUnannotated property
            double originalStrokeThicknessUnannotated = viewModel.StrokeThicknessUnannotated;
            viewModel.StrokeThicknessUnannotated = originalStrokeThicknessUnannotated + 0.1;
            Assert.That(viewModel.StrokeThicknessUnannotated, Is.EqualTo(originalStrokeThicknessUnannotated + 0.1));
            viewModel.StrokeThicknessUnannotated = originalStrokeThicknessUnannotated;
        }

        [Test]
        public static void TestModTypeForTreeView()
        {
            var modGroups = GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType).ToList();
            var key = modGroups.First().Key;
            ModTypeForTreeViewModel modTypeForTreeView = new(key, false);
            Assert.That(!modTypeForTreeView.Expanded);
            Assert.That(modTypeForTreeView.DisplayName == key);
            Assert.That(((SolidColorBrush)modTypeForTreeView.Background).Color ==
                        new SolidColorBrush(Colors.Transparent).Color);
            Assert.That(modTypeForTreeView.Children != null);

            modTypeForTreeView = new(modGroups.First().Key, true);
            Assert.That(((SolidColorBrush)modTypeForTreeView.Background).Color ==
                        new SolidColorBrush(Colors.Red).Color);

            modGroups.First().Select(p =>
                    new ModForTreeViewModel(p.ToString(), false, p.IdWithMotif, false, modTypeForTreeView))
                .ForEach(mod => modTypeForTreeView.Children.Add(mod));
            Assert.That(modTypeForTreeView.Children.Count == modGroups.First().Count());
            Assert.That(modTypeForTreeView.Children.All(p => p.Parent == modTypeForTreeView));
            Assert.That(modTypeForTreeView.Children.All(p => p.Use == false));
            modTypeForTreeView.VerifyCheckState();
            Assert.That(modTypeForTreeView.Use == false);

            modTypeForTreeView.Children.First().Use = true;
            modTypeForTreeView.VerifyCheckState();
            Assert.That(modTypeForTreeView.Use == null);

            modTypeForTreeView.Children.ForEach(mod => mod.Use = true);
            modTypeForTreeView.VerifyCheckState();
            Assert.That(modTypeForTreeView.Use == true);
        }

        [Test]
        public static void TestModForTreeView()
        {
            var modGroup = GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType).First();
            var mod = modGroup.First();
            ModTypeForTreeViewModel modTypeForTreeView = new(modGroup.Key, false);
            ModForTreeViewModel modForTreeView = new(mod.ToString(), false, mod.IdWithMotif, false, modTypeForTreeView);
            Assert.That(modForTreeView.ModName == mod.IdWithMotif);
            Assert.That(modForTreeView.DisplayName == mod.IdWithMotif);
            Assert.That(!modForTreeView.Use);
            Assert.That(modForTreeView.ToolTipStuff == mod.ToString());
            var color = MetaDrawSettings.ModificationTypeToColor[mod.IdWithMotif];
            Assert.That(modForTreeView.SelectedColor == color.GetColorName());
            Assert.That(modForTreeView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(color).Color);
        }

        [Test]
        public static void TestIonTypeForTreeView()
        {
            var ions = (ProductType[])Enum.GetValues(typeof(ProductType));
            IonTypeForTreeViewModel ionForTreeViews = new("Common Ions", ions, false);
            Assert.That(ionForTreeViews.GroupName == "Common Ions");
            Assert.That(ionForTreeViews.Ions.Count ==
                        ions.Length + 2); // magic number +2 is for the internal ion color and background peak color
            Assert.That(!ionForTreeViews.Ions.Any(p => p.IsBeta));
            ionForTreeViews = new("Common Ions", ions, true);
            Assert.That(ionForTreeViews.Ions.Any(p => p.IsBeta));
        }

        [Test]
        public static void TestIonForTreeView()
        {
            foreach (var ion in Enum.GetValues<ProductType>())
            {
                IonForTreeViewModel ionForTreeView = new(ion, false);
                Assert.That(!ionForTreeView.IsBeta);
                Assert.That(ionForTreeView.IonType == ion);

                var color = MetaDrawSettings.ProductTypeToColor[ion];
                var name = color.GetColorName();
                Assert.That(ionForTreeView.SelectedColor.Replace(" ", "") == name);
                Assert.That(ionForTreeView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(color).Color);

                ionForTreeView = new(ion, true);
                Assert.That(ionForTreeView.IsBeta);
                Assert.That(ionForTreeView.IonType == ion);

                color = MetaDrawSettings.BetaProductTypeToColor[ion];
                name = color.GetColorName();
                Assert.That(ionForTreeView.SelectedColor.Replace(" ", "") == name);
                Assert.That(ionForTreeView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(color).Color);
            }
        }

        [Test]
        public static void TestPtmLegendViewModels()
        {
            string psmsPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TopDownTestData\TDGPTMDSearchResults.psmtsv");
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmsPath, out List<string> warnings)
                .Where(p => p.AmbiguityLevel == "1").ToList();
            PsmFromTsv psm = psms.First(p =>
                new PeptideWithSetModifications(p.FullSequence, GlobalVariables.AllModsKnownDictionary)
                    .AllModsOneIsNterminus.Values.Distinct().Count() == 2);
            PeptideWithSetModifications pepWithSetMods = new(psm.FullSequence, GlobalVariables.AllModsKnownDictionary);
            var twoMods = pepWithSetMods.AllModsOneIsNterminus.Values.ToList();

            PtmLegendViewModel PtmLegendView = new PtmLegendViewModel(psm, 100);
            PtmLegendView.Visibility = false;
            Assert.That(PtmLegendView.Header == "Legend");
            Assert.That(PtmLegendView.HeaderSize == 12);
            Assert.That(PtmLegendView.TopOffset == 100);
            Assert.That(PtmLegendView.LegendItemViewModels.Count == 2);
            Assert.That(PtmLegendView.LegendItemViewModels.First().Name == twoMods.First().IdWithMotif);
            Assert.That(PtmLegendView.LegendItemViewModels.First().ColorBrush.Color == DrawnSequence
                .ParseColorBrushFromOxyColor(MetaDrawSettings.ModificationTypeToColor[twoMods.First().IdWithMotif])
                .Color);
            Assert.That(PtmLegendView.LegendItemViewModels.First().Name == twoMods.First().IdWithMotif);
            var mod = twoMods.First();
            PtmLegendItemViewModel ptmLegendItemView = new(mod.IdWithMotif);
            Assert.That(ptmLegendItemView.Name == mod.IdWithMotif);
            Assert.That(ptmLegendItemView.ColorBrush.Color == DrawnSequence
                .ParseColorBrushFromOxyColor(MetaDrawSettings.ModificationTypeToColor[mod.IdWithMotif]).Color);

            // test that residue per segment incrementation works and cannot be less than 1
            int residuesPerSegment = PtmLegendView.ResiduesPerSegment;
            PtmLegendView.IncreaseResiduesPerSegment();
            Assert.That(residuesPerSegment + 1 == PtmLegendView.ResiduesPerSegment);
            PtmLegendView.DecreaseResiduesPerSegment();
            Assert.That(residuesPerSegment == PtmLegendView.ResiduesPerSegment);
            PtmLegendView.ResiduesPerSegment = 1;
            try
            {
                PtmLegendView.DecreaseResiduesPerSegment();
                Assert.That(false);
            }
            catch (IndexOutOfRangeException e)
            {
                Assert.That(e.Message == ("ResiduesPerSegment cannot be less than one"));
            }
            catch (Exception)
            {
                Assert.That(false);
            }

            // test that segments per row incrementation works and cannot be less than 1
            int segmentsPerRow = PtmLegendView.SegmentsPerRow;
            PtmLegendView.IncreaseSegmentsPerRow();
            Assert.That(segmentsPerRow + 1 == PtmLegendView.SegmentsPerRow);
            PtmLegendView.DecreaseSegmentsPerRow();
            Assert.That(segmentsPerRow == PtmLegendView.SegmentsPerRow);
            PtmLegendView.SegmentsPerRow = 1;
            try
            {
                PtmLegendView.DecreaseSegmentsPerRow();
                Assert.That(false);
            }
            catch (IndexOutOfRangeException e)
            {
                Assert.That(e.Message == ("SegmentsPerRow cannot be less than one"));
            }
            catch (Exception)
            {
                Assert.That(false);
            }
        }

        [Test]
        public static void TestDrawnSequenceColorConversions()
        {
            OxyColor oxyBlue = OxyColors.Blue;
            Color colorBlue = Colors.Blue;
            SolidColorBrush brushBlue = new SolidColorBrush(colorBlue);

            var colorBrushFromOxy = DrawnSequence.ParseColorBrushFromOxyColor(oxyBlue);
            Assert.That(colorBrushFromOxy.Color == brushBlue.Color);

            var colorBrushfromName = DrawnSequence.ParseColorBrushFromName(oxyBlue.GetColorName());
            Assert.That(colorBrushfromName.Color == brushBlue.Color);
            var colorBrushfromNameBad = DrawnSequence.ParseColorBrushFromName("humbug");
            Assert.That(colorBrushfromNameBad.Color == Colors.Aqua);

            var oxyFromName = DrawnSequence.ParseOxyColorFromName(oxyBlue.GetColorName());
            Assert.That(oxyFromName == oxyBlue);
            var oxyFromNameBad = DrawnSequence.ParseOxyColorFromName("gobbledygook");
            Assert.That(oxyFromNameBad == MetaDrawSettings.FallbackColor);

            var colorFromOxy = DrawnSequence.ParseColorFromOxyColor(oxyBlue);
            Assert.That(colorFromOxy == colorBlue);
        }

        [Test]
        public static void TestSpectrumDescriptorViewModelSyncsWithStaticSettings()
        {
            // Arrange: Reset settings and create a new view model
            MetaDrawSettings.ResetSettings();
            var vm = new MetaDrawSettingsViewModel(false);

            // Assert: Each SpectrumDescriptorViewModel reflects the static settings
            foreach (var descVm in vm.SpectrumDescriptors)
            {
                // The initial IsSelected should match the static MetaDrawSettings.SpectrumDescription
                Assert.That(descVm.IsSelected, Is.EqualTo(MetaDrawSettings.SpectrumDescription[descVm.DisplayName + ": "]));
            }

            // Act: Change IsSelected in the view model for each descriptor
            foreach (var descVm in vm.SpectrumDescriptors)
            {
                bool newValue = !descVm.IsSelected;
                descVm.IsSelected = newValue;

                // Assert: The static settings are updated accordingly
                Assert.That(MetaDrawSettings.SpectrumDescription[descVm.DisplayName + ": "], Is.EqualTo(newValue));
                Assert.That(descVm.IsSelected, Is.EqualTo(newValue));
            }
        }
    }
}
  
