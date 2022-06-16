using EngineLayer;
using GuiFunctions;
using IO.MzML;
using NUnit.Framework;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace Test
{
    public static class MetaDrawSettingsAndViewsTest
    {
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
            
            MetaDrawSettings.ShowContaminants = true;
            MetaDrawSettings.AnnotateMzValues = false;
            snapshot = MetaDrawSettings.MakeSnapShot();
            Assert.That(snapshot.ShowContaminants.Equals(MetaDrawSettings.ShowContaminants));
            Assert.That(snapshot.AnnotateMzValues.Equals(MetaDrawSettings.AnnotateMzValues));
            Assert.That(snapshot.QValueFilter.Equals(MetaDrawSettings.QValueFilter));
            Assert.That(snapshot.LocalizationLevelStart.Equals(MetaDrawSettings.LocalizationLevelStart));
            Assert.That(snapshot.ExportType.Equals(MetaDrawSettings.ExportType));
            var colorValues = MetaDrawSettings.ProductTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            var betaColorValues = MetaDrawSettings.BetaProductTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            var modificationColorValues = MetaDrawSettings.ModificationTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            var coverageColorValues = MetaDrawSettings.CoverageTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            var spectrumDescriptionValues = MetaDrawSettings.SpectrumDescription.Values.ToList();
            Assert.That(snapshot.ProductTypeToColorValues.Except(colorValues).Count() == 0);
            Assert.That(snapshot.BetaProductTypeToColorValues.Except(betaColorValues).Count() == 0);
            Assert.That(snapshot.ModificationTypeToColorValues.Except(modificationColorValues).Count() == 0);
            Assert.That(snapshot.CoverageTypeToColorValues.Except(coverageColorValues).Count() == 0);
            Assert.That(snapshot.SpectrumDescriptionValues.Except(spectrumDescriptionValues).Count() == 0);

            snapshot.QValueFilter = 0.5;
            snapshot.AnnotateCharges = true;
            snapshot.DisplayIonAnnotations = true;
            snapshot.AnnotateMzValues = true;
            snapshot.AnnotationBold = true;
            snapshot.ShowContaminants = false;
            snapshot.ShowDecoys = true;
            snapshot.ShowLegend = false;
            snapshot.DrawNumbersUnderStationary = false;
            MetaDrawSettings.LoadSettings(snapshot);
            Assert.That(snapshot.DisplayIonAnnotations.Equals(MetaDrawSettings.DisplayIonAnnotations));
            Assert.That(snapshot.AnnotateMzValues.Equals(MetaDrawSettings.AnnotateMzValues));
            Assert.That(snapshot.AnnotateCharges.Equals(MetaDrawSettings.AnnotateCharges));
            Assert.That(snapshot.AnnotationBold.Equals(MetaDrawSettings.AnnotationBold));
            Assert.That(snapshot.ShowDecoys.Equals(MetaDrawSettings.ShowDecoys));
            Assert.That(snapshot.ShowContaminants.Equals(MetaDrawSettings.ShowContaminants));
            Assert.That(snapshot.QValueFilter.Equals(MetaDrawSettings.QValueFilter));
            Assert.That(snapshot.LocalizationLevelStart.Equals(MetaDrawSettings.LocalizationLevelStart));
            Assert.That(snapshot.LocalizationLevelEnd.Equals(MetaDrawSettings.LocalizationLevelEnd));
            colorValues = MetaDrawSettings.ProductTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            betaColorValues = MetaDrawSettings.BetaProductTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            modificationColorValues = MetaDrawSettings.ModificationTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            coverageColorValues = MetaDrawSettings.CoverageTypeToColor.Values.Select(p => p.GetColorName()).ToList();
            spectrumDescriptionValues = MetaDrawSettings.SpectrumDescription.Values.ToList();
            Assert.That(snapshot.ProductTypeToColorValues.Except(colorValues).Count() == 0);
            Assert.That(snapshot.BetaProductTypeToColorValues.Except(betaColorValues).Count() == 0);
            Assert.That(snapshot.ModificationTypeToColorValues.Except(modificationColorValues).Count() == 0);
            Assert.That(snapshot.CoverageTypeToColorValues.Except(coverageColorValues).Count() == 0);
            Assert.That(snapshot.SpectrumDescriptionValues.Except(spectrumDescriptionValues).Count() == 0);
        }

        [Test]
        public static void TestSaveAndLoadDefaultSettings()
        {
            SettingsViewModel model = new SettingsViewModel(false);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawWithSpectraLibrary");
            Directory.CreateDirectory(outputFolder);

            SettingsViewModel.SettingsPath = Path.Combine(outputFolder, @"MetaDrawSettingsDefault.xml");
            Assert.That(model.HasDefaultSaved == false);
            try
            {
                model.LoadSettings();
                Assert.Fail();
            }
            catch (Exception e) { }

            Assert.That(model.Modifications.First().Children.First().SelectedColor == "Green");
            model.Modifications.First().Children.First().SelectionChanged("Blue");
            Assert.That(model.Modifications.First().Children.First().SelectedColor == "Blue");
            model.SaveAsDefault();
            Assert.That(model.HasDefaultSaved == true);
            model.LoadSettings();

            SettingsViewModel model2 = new();
            Assert.That(model2.Modifications.First().Children.First().SelectedColor == "Blue");

            SettingsViewModel model3 = new(false);
            Assert.That(model3.Modifications.First().Children.First().SelectedColor == "Blue");

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestSettingsViewLoading()
        {
            SettingsViewModel BlankSettingsView = new SettingsViewModel(false);
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
        public static void TestSettingsViewSaveAndChildSelectionChanged()
        {
            SettingsViewModel view = new SettingsViewModel(false);
            Assert.That(!view.IonGroups.First().Ions.First().HasChanged);
            view.IonGroups.First().Ions.First().SelectionChanged("Blue");
            Assert.That(view.IonGroups.First().Ions.First().HasChanged);
            Assert.That(view.IonGroups.First().Ions.First().SelectedColor == "Blue");
            Assert.That(view.IonGroups.First().Ions.First().ColorBrush.Color == DrawnSequence.ParseColorBrushFromName("Blue").Color);

            Assert.That(!view.Modifications.First().Children.First().HasChanged);
            view.Modifications.First().Children.First().SelectionChanged("Blue");
            Assert.That(view.Modifications.First().Children.First().HasChanged);
            Assert.That(view.Modifications.First().Children.First().SelectedColor == "Blue");
            Assert.That(view.Modifications.First().Children.First().ColorBrush.Color == DrawnSequence.ParseColorBrushFromName("Blue").Color);

            Assert.That(!view.CoverageColors.First().HasChanged);
            view.CoverageColors.First().SelectionChanged("Blue");
            Assert.That(view.CoverageColors.First().HasChanged);
            Assert.That(view.CoverageColors.First().SelectedColor == "Blue");
            Assert.That(view.CoverageColors.First().ColorBrush.Color == DrawnSequence.ParseColorBrushFromName("Blue").Color);

            view.Save();
            Assert.That(MetaDrawSettings.ProductTypeToColor[view.IonGroups.First().Ions.First().IonType] == OxyColors.Blue);
            Assert.That(MetaDrawSettings.ModificationTypeToColor[view.Modifications.First().Children.First().ModName] == OxyColors.Blue);
            Assert.That(MetaDrawSettings.CoverageTypeToColor[view.CoverageColors.First().Name] == OxyColors.Blue);
        }

        [Test]
        public static void TestCoverageTypeForTreeView()
        {
            CoverageTypeForTreeViewModel coverageTypeForTreeView = new("N-Terminal Color");
            Assert.That(coverageTypeForTreeView.Name == "N-Terminal Color");
            var color = MetaDrawSettings.CoverageTypeToColor["N-Terminal Color"];
            Assert.That(coverageTypeForTreeView.SelectedColor == color.GetColorName());
            Assert.That(coverageTypeForTreeView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(color).Color);
        }

        [Test]
        public static void TestModTypeForTreeView()
        {
            var modGroups = GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType);
            var key = modGroups.First().Key;
            ModTypeForTreeViewModel modTypeForTreeView = new(key, false);
            Assert.That(!modTypeForTreeView.Expanded);
            Assert.That(modTypeForTreeView.DisplayName == key);
            Assert.That(((SolidColorBrush)modTypeForTreeView.Background).Color == new SolidColorBrush(Colors.Transparent).Color);
            Assert.That(modTypeForTreeView.Children != null);

            modTypeForTreeView = new(modGroups.First().Key, true);
            Assert.That(((SolidColorBrush)modTypeForTreeView.Background).Color == new SolidColorBrush(Colors.Red).Color);
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
            Assert.That(ionForTreeViews.Ions.Count == ions.Length);
            Assert.That(!ionForTreeViews.Ions.Any(p => p.IsBeta));
            ionForTreeViews = new("Common Ions", ions, true);
            Assert.That(ionForTreeViews.Ions.Any(p => p.IsBeta));
        }

        [Test]
        public static void TestIonForTreeView()
        {
            var ion = ((ProductType[])Enum.GetValues(typeof(ProductType))).First();
            IonForTreeViewModel ionForTreeView = new(ion, false);
            Assert.That(!ionForTreeView.IsBeta);
            Assert.That(ionForTreeView.IonType == ion);
            var color = MetaDrawSettings.ProductTypeToColor[ion];
            Assert.That(ionForTreeView.SelectedColor == color.GetColorName());
            Assert.That(ionForTreeView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(color).Color);

            ionForTreeView = new(ion, true);
            Assert.That(ionForTreeView.IsBeta);
            color = MetaDrawSettings.BetaProductTypeToColor[ion];
            Assert.That(ionForTreeView.SelectedColor == color.GetColorName());
            Assert.That(ionForTreeView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(color).Color);
        }

        [Test]
        public static void TestPtmLegendViews()
        {
            var modGroup = GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType).First();
            var twoMods = modGroup.Take(2).ToList();
            PtmLegendViewModel PtmLegendView = new PtmLegendViewModel(twoMods);
            Assert.That(PtmLegendView.Header == "Legend");
            Assert.That(PtmLegendView.HeaderSize == 12);
            Assert.That(PtmLegendView.LegendItems.Count == 2);
            Assert.That(PtmLegendView.LegendItems.First().ModName == twoMods.First().IdWithMotif);
            Assert.That(PtmLegendView.LegendItems.First().ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(MetaDrawSettings.ModificationTypeToColor[twoMods.First().IdWithMotif]).Color);
            Assert.That(PtmLegendView.LegendItems.First().ModName == twoMods.First().IdWithMotif);
            var mod = twoMods.First();
            PtmLegendItemViewModel ptmLegendItemView = new(mod.IdWithMotif);
            Assert.That(ptmLegendItemView.ModName == mod.IdWithMotif);
            Assert.That(ptmLegendItemView.ColorBrush.Color == DrawnSequence.ParseColorBrushFromOxyColor(MetaDrawSettings.ModificationTypeToColor[mod.IdWithMotif]).Color);

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
            catch (Exception e)
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
            catch (Exception e)
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

            var oxyFromName = DrawnSequence.ParseOxyColorFromName(oxyBlue.GetColorName());
            Assert.That(oxyFromName == oxyBlue);

            var colorFromOxy = DrawnSequence.ParseColorFromOxyColor(oxyBlue);
            Assert.That(colorFromOxy == colorBlue);
        }
    }
}
