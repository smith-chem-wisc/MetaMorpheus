using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Media;
using EngineLayer;
using GuiFunctions;
using GuiFunctions.ViewModels.Legends;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using OxyPlot;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;

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

            MetaDrawSettings.ShowContaminants = true;
            MetaDrawSettings.AnnotateMzValues = false;
            snapshot = MetaDrawSettings.MakeSnapShot();
            Assert.That(snapshot.ShowContaminants.Equals(MetaDrawSettings.ShowContaminants));
            Assert.That(snapshot.AnnotateMzValues.Equals(MetaDrawSettings.AnnotateMzValues));
            Assert.That(snapshot.QValueFilter.Equals(MetaDrawSettings.QValueFilter));
            Assert.That(snapshot.LocalizationLevelStart.Equals(MetaDrawSettings.LocalizationLevelStart));
            Assert.That(snapshot.ExportType.Equals(MetaDrawSettings.ExportType));

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
            MetaDrawSettingsViewModel BlankSettingsView = new MetaDrawSettingsViewModel(false);
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
        public static void TestModTypeForTreeView()
        {
            var modGroups = GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType);
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
        public static void TestChimeraLegendViewModels()
        {
            // object setup
            string psmsPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"TopDownTestData\TDGPTMDSearchResults.psmtsv");
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmsPath, out List<string> warnings);
            Assert.That(warnings.Count, Is.EqualTo(0));
            List<PsmFromTsv> filteredChimeras =
                psms.Where(p => p.QValue <= 0.01 && p.PEP <= 0.5 && p.PrecursorScanNum == 1557).ToList();
            Assert.That(filteredChimeras.Count, Is.EqualTo(3));

            // test chimera legend basic functionality
            ChimeraLegendViewModel chimeraLegend = new ChimeraLegendViewModel(filteredChimeras);
            Assert.That(chimeraLegend.ChimeraLegendItems.Count == 2);
            Assert.That(chimeraLegend.TopOffset == 0);
            Assert.That(chimeraLegend.DisplaySharedIonLabel == true);
            Assert.That(chimeraLegend.Visibility == true);
            Assert.That(chimeraLegend.ChimeraLegendItems.Values.First().Count == 3);
            Assert.That(chimeraLegend.ChimeraLegendItems.Values.ToList()[1].Count == 1);

            // test chimera legend overflow colors
            // more unique proteins than colored
            List<PsmFromTsv> overflowInducingProteins = psms.DistinctBy(p => p.BaseSeq)
                .Take(ChimeraSpectrumMatchPlot.ColorByProteinDictionary.Keys.Count + 1).ToList();
            chimeraLegend = new(overflowInducingProteins);
            Assert.That(chimeraLegend.ChimeraLegendItems.Values.DistinctBy(p =>
                p.Select(m => m.ColorBrush.Color)).Count(), Is.EqualTo(overflowInducingProteins.Count()));
            Assert.That(chimeraLegend.ChimeraLegendItems.First().Value.First().ColorBrush.Color !=
                        chimeraLegend.ChimeraLegendItems[overflowInducingProteins[1].BaseSeq].First().ColorBrush.Color);
            Assert.That(chimeraLegend.ChimeraLegendItems.First().Value.First().ColorBrush.Color ==
                        chimeraLegend.ChimeraLegendItems.Last().Value.First().ColorBrush.Color);

            // more unique proteoforms than colored
            overflowInducingProteins = psms
                .Take(ChimeraSpectrumMatchPlot.ColorByProteinDictionary.First().Value.Count)
                .Select(p => p = new(p, overflowInducingProteins.First().FullSequence, 0,
                    overflowInducingProteins.First().BaseSeq)).ToList();
            Assert.That(overflowInducingProteins.All(p => p.BaseSeq == overflowInducingProteins.First().BaseSeq));
            Assert.That(overflowInducingProteins.All(p =>
                p.FullSequence == overflowInducingProteins.First().FullSequence));
            chimeraLegend = new(overflowInducingProteins);
            Assert.That(chimeraLegend.ChimeraLegendItems.First().Value.DistinctBy(p => p.ColorBrush.Color).Count(),
                Is.EqualTo(overflowInducingProteins.Count() + 1));
            Assert.That(chimeraLegend.ChimeraLegendItems.First().Value.Count() == overflowInducingProteins.Count + 1);
            Assert.That(chimeraLegend.ChimeraLegendItems.First().Value.Last().ColorBrush.Color, Is.EqualTo(DrawnSequence
                .ParseColorBrushFromOxyColor(ChimeraSpectrumMatchPlot.OverflowColors.Dequeue()).Color));

            // test chimera legend item
            ChimeraLegendItemViewModel chimeraLegendItem = new("tacos", OxyColors.Chocolate);
            Assert.That(chimeraLegendItem.Name == "tacos");
            Assert.That(chimeraLegendItem.ColorBrush.Color ==
                        DrawnSequence.ParseColorBrushFromOxyColor(OxyColors.Chocolate).Color);
            chimeraLegendItem = new("", OxyColors.Chocolate);
            Assert.That(chimeraLegendItem.Name == "No Modifications");
            chimeraLegendItem = new(null, OxyColors.Chocolate);
            Assert.That(chimeraLegendItem.Name == "No Modifications");

            chimeraLegend = new ChimeraLegendViewModel(new List<PsmFromTsv>() { psms.First() });
            Assert.That(chimeraLegend.DisplaySharedIonLabel == false);
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
    }
}
  
