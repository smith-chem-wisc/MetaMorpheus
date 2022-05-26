using EngineLayer;
using EngineLayer.GlycoSearch;
using GuiFunctions;
using Nett;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawSettingsWindow.xaml
    /// </summary>
    public partial class MetaDrawSettingsWindow : Window
    {
        private readonly ObservableCollection<ModTypeForTreeView> Modifications = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<IonTypeForTreeViewModel> IonGroups = new ObservableCollection<IonTypeForTreeViewModel>();

        public MetaDrawSettingsWindow()
        {
            InitializeComponent();
            PopulateChoices();
        }

        private void PopulateChoices()
        {
            foreach (string level in System.Enum.GetNames(typeof(LocalizationLevel)))
            {
                CmbGlycanLocalizationLevelStart.Items.Add(level);
                CmbGlycanLocalizationLevelEnd.Items.Add(level);
            }

            DisplayAnnotationsCheckBox.IsChecked = MetaDrawSettings.DisplayIonAnnotations;
            MZCheckBox.IsChecked = MetaDrawSettings.AnnotateMzValues;
            ChargesCheckBox.IsChecked = MetaDrawSettings.AnnotateCharges;
            BoldTextCheckBox.IsChecked = MetaDrawSettings.AnnotationBold;
            DecoysCheckBox.IsChecked = MetaDrawSettings.ShowDecoys;
            ContaminantsCheckBox.IsChecked = MetaDrawSettings.ShowContaminants;
            PrecursorChargeCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Charge: "];
            PrecursorMassCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Mass: "];
            TheoreticalMassCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Theoretical Mass: "];
            ScoreCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Score: "];
            ProteinAccessionCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Protein Accession: "];
            ProteinCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Protein: "];
            DecoyContaminantTargetCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "];
            QValueCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Q-Value: "];
            SequenceLengthCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Sequence Length: "];
            ProFormaLevelCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["ProForma Level: "];
            PEPCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["PEP: "];
            PEPQValueCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["PEP Q-Value: "];
            qValueBox.Text = MetaDrawSettings.QValueFilter.ToString();
            TextSizeBox.Text = MetaDrawSettings.AnnotatedFontSize.ToString();
            CmbGlycanLocalizationLevelStart.SelectedItem = MetaDrawSettings.LocalizationLevelStart.ToString();
            CmbGlycanLocalizationLevelEnd.SelectedItem = MetaDrawSettings.LocalizationLevelEnd.ToString();


            var ions = ((ProductType[])Enum.GetValues(typeof(ProductType)));
            var common = ions.Where(p => p.ToString().Equals("a") || p.ToString().Equals("b") || p.ToString().Equals("c")
                                          || p.ToString().Equals("x") || p.ToString().Equals("y") || p.ToString().Equals("zDot"));
            var lessCommon = ions.Where(p => !p.ToString().Equals("a") || !p.ToString().Equals("b") || !p.ToString().Equals("c")
                                          || !p.ToString().Equals("x") || !p.ToString().Equals("y") || !p.ToString().Equals("zDot"));
            IonGroups.Add(new IonTypeForTreeViewModel("Common Ions", common, false));
            IonGroups.Add(new IonTypeForTreeViewModel("Less Common Ions", lessCommon, false));
            IonGroups.Add(new IonTypeForTreeViewModel("Cross Linked Beta Peptide", ions, true));

            TestIonColorTreeView.DataContext = IonGroups;

            foreach (var modGroup in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(modGroup.Key, false);
                Modifications.Add(theModType);
                foreach (var mod in modGroup)
                {
                    theModType.Children.Add(new ModForTreeView(mod.ToString(), false, mod.IdWithMotif, false, theModType));
                }
            }
            gptmdModsTreeView.DataContext = Modifications;




            #region Colors

            aIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            bIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            cIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            xIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            yIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            zIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            aStarIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            aDegreeIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            bStarIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            bDegreeIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            yStarIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            yDegreeIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            zPlusOneIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            mIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            dIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            yCoreIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            YIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            aIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.a]];
            bIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.b]];
            cIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.c]];
            xIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.x]];
            yIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.y]];
            zIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.zDot]];
            aStarIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.aStar]];
            aDegreeIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.aDegree]];
            bStarIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.bStar]];
            bDegreeIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.bDegree]];
            yStarIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.yStar]];
            yDegreeIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.yDegree]];
            zPlusOneIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.zPlusOne]];
            mIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.M]];
            dIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.D]];
            yCoreIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.Ycore]];
            YIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.ProductTypeToColor[ProductType.Y]];

            BetaaIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetabIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetacIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetaxIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetayIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetazIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetaaStarIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetaaDegreeIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetabStarIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetabDegreeIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetayStarIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetayDegreeIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetazPlusOneIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetamIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetadIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetayCoreIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetaYIonComboBox.ItemsSource = MetaDrawSettings.PossibleColors.Values;
            BetaaIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.a]];
            BetabIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.b]];
            BetacIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.c]];
            BetaxIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.x]];
            BetayIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.y]];
            BetazIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.zDot]];
            BetaaStarIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.aStar]];
            BetaaDegreeIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.aDegree]];
            BetabStarIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.bStar]];
            BetabDegreeIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.bDegree]];
            BetayStarIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.yStar]];
            BetayDegreeIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.yDegree]];
            BetazPlusOneIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.zPlusOne]];
            BetamIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.M]];
            BetadIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.D]];
            BetayCoreIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.Ycore]];
            BetaYIonComboBox.SelectedItem = MetaDrawSettings.PossibleColors[MetaDrawSettings.BetaProductTypeToColor[ProductType.Y]];
            #endregion
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawSettings.DisplayIonAnnotations = DisplayAnnotationsCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotateMzValues = MZCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotateCharges = ChargesCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotationBold = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowDecoys = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowContaminants = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Precursor Charge: "] = PrecursorChargeCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Precursor Mass: "] = PrecursorMassCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Theoretical Mass: "] = TheoreticalMassCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Score: "] = ScoreCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Protein Accession: "] = ProteinAccessionCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Protein: "] = ProteinCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "] = DecoyContaminantTargetCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Q-Value: "] = QValueCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Sequence Length: "] = SequenceLengthCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["ProForma Level: "] = ProFormaLevelCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["PEP: "] = PEPCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["PEP Q-Value: "] = PEPQValueCheckBox.IsChecked.Value;
            MetaDrawSettings.LocalizationLevelStart = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelStart.SelectedItem.ToString());
            MetaDrawSettings.LocalizationLevelEnd = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelEnd.SelectedItem.ToString());

            #region Colors

            MetaDrawSettings.ProductTypeToColor[ProductType.a] = MetaDrawSettings.NameToOxyColorConverter(aIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.b] = MetaDrawSettings.NameToOxyColorConverter(bIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.c] = MetaDrawSettings.NameToOxyColorConverter(cIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.x] = MetaDrawSettings.NameToOxyColorConverter(xIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.y] = MetaDrawSettings.NameToOxyColorConverter(yIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.zDot] = MetaDrawSettings.NameToOxyColorConverter(zIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.aStar] = MetaDrawSettings.NameToOxyColorConverter(aStarIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.aDegree] = MetaDrawSettings.NameToOxyColorConverter(aDegreeIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.bStar] = MetaDrawSettings.NameToOxyColorConverter(bStarIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.bDegree] = MetaDrawSettings.NameToOxyColorConverter(bDegreeIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.yStar] = MetaDrawSettings.NameToOxyColorConverter(yStarIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.yDegree] = MetaDrawSettings.NameToOxyColorConverter(yDegreeIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.zPlusOne] = MetaDrawSettings.NameToOxyColorConverter(zPlusOneIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.M] = MetaDrawSettings.NameToOxyColorConverter(mIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.D] = MetaDrawSettings.NameToOxyColorConverter(dIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.Ycore] = MetaDrawSettings.NameToOxyColorConverter(yCoreIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.ProductTypeToColor[ProductType.Y] = MetaDrawSettings.NameToOxyColorConverter(yIonComboBox.SelectedItem.ToString());

            MetaDrawSettings.BetaProductTypeToColor[ProductType.a] = MetaDrawSettings.NameToOxyColorConverter(BetaaIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.b] = MetaDrawSettings.NameToOxyColorConverter(BetabIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.c] = MetaDrawSettings.NameToOxyColorConverter(BetacIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.x] = MetaDrawSettings.NameToOxyColorConverter(BetaxIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.y] = MetaDrawSettings.NameToOxyColorConverter(BetayIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.zDot] = MetaDrawSettings.NameToOxyColorConverter(BetazIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.aStar] = MetaDrawSettings.NameToOxyColorConverter(BetaaStarIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.aDegree] = MetaDrawSettings.NameToOxyColorConverter(BetaaDegreeIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.bStar] = MetaDrawSettings.NameToOxyColorConverter(BetabStarIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.bDegree] = MetaDrawSettings.NameToOxyColorConverter(BetabDegreeIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.yStar] = MetaDrawSettings.NameToOxyColorConverter(BetayStarIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.yDegree] = MetaDrawSettings.NameToOxyColorConverter(BetayDegreeIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.zPlusOne] = MetaDrawSettings.NameToOxyColorConverter(BetazPlusOneIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.M] = MetaDrawSettings.NameToOxyColorConverter(BetamIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.D] = MetaDrawSettings.NameToOxyColorConverter(BetadIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.Ycore] = MetaDrawSettings.NameToOxyColorConverter(BetayCoreIonComboBox.SelectedItem.ToString());
            MetaDrawSettings.BetaProductTypeToColor[ProductType.Y] = MetaDrawSettings.NameToOxyColorConverter(BetayIonComboBox.SelectedItem.ToString());

            foreach (var group in IonGroups)
            {
                foreach (var ion in group)
                {
                    if (ion.HasChanged)
                    {
                        if (ion.IsBeta)
                            MetaDrawSettings.BetaProductTypeToColor[ion.IonType] = MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(ion.SelectedColor.Replace(" ", ""))).First();
                        else
                            MetaDrawSettings.ProductTypeToColor[ion.IonType] = MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(ion.SelectedColor.Replace(" ", ""))).First();
                    }
                }
            }





            #endregion

            if (!string.IsNullOrWhiteSpace(qValueBox.Text))
            {
                if (double.TryParse(qValueBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out double qValueFilter) && qValueFilter >= 0 && qValueFilter <= 1)
                {
                    MetaDrawSettings.QValueFilter = qValueFilter;
                }
                else
                {
                    MessageBox.Show("Could not parse q-value filter; must be number between 0 and 1 inclusive");
                    return;
                }
            }
            else
            {
                MetaDrawSettings.QValueFilter = 1;
            }

            if (!string.IsNullOrWhiteSpace(TextSizeBox.Text))
            {
                if (int.TryParse(TextSizeBox.Text, out int fontSize))
                {
                    if (fontSize > 15)
                    {
                        MessageBox.Show("Font size must be <= 15");
                        return;
                    }

                    MetaDrawSettings.AnnotatedFontSize = fontSize;
                }
                else
                {
                    MessageBox.Show("Could not parse font size; must be a positive integer");
                    return;
                }
            }
            else
            {
                MetaDrawSettings.AnnotatedFontSize = 12;
            }

            DialogResult = true;
        }

        private void setDefaultbutton_Click(object sender, RoutedEventArgs e)
        {
            Save_Click(sender, e);
            MetaDrawSettingsSnapshot settings = MetaDrawSettings.MakeSnapShot();
            Toml.WriteFile<MetaDrawSettingsSnapshot>(settings, Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"MetaDrawSettingsDefault.toml"));
        }

        private void ComboBox_SelectionChanged(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            ((IonForTreeViewModel)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
        }
    }
}
