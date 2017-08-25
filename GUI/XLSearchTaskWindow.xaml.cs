using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class XLSearchTaskWindow : Window
    {
        #region Private Fields

        private readonly DataContextForSearchTaskWindow dataContextForSearchTaskWindow;

        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();

        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();

        #endregion Private Fields

        #region Public Constructors

        public XLSearchTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new XLSearchTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the XLSearch Task";

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                //ModExpanderTitle =
                //"fixed: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
                //+ " variable: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
                //+ " localize: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = dataContextForSearchTaskWindow;
        }

        public XLSearchTaskWindow(XLSearchTask task)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = task;
            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                //ModExpanderTitle =
                //"fixed: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
                //+ " variable: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
                //+ " localize: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = dataContextForSearchTaskWindow;
        }

        #endregion Public Constructors

        #region Internal Properties

        internal XLSearchTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void PopulateChoices()
        {
            foreach (string crosslinkerName in Enum.GetNames(typeof(CrosslinkerType)))
                cbCrosslinker.Items.Add(crosslinkerName);

            cbbXLprecusorMsTl.Items.Add("Absolute");
            cbbXLprecusorMsTl.Items.Add("Ppm");

            cbbXLBetaprecusorMsTl.Items.Add("Absolute");
            cbbXLBetaprecusorMsTl.Items.Add("Ppm");

            foreach (Protease protease in GlobalEngineLevelSettings.ProteaseDictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            productMassToleranceComboBox.Items.Add("Absolute");
            productMassToleranceComboBox.Items.Add("Ppm");

            //foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
            //    productMassToleranceComboBox.Items.Add(toleranceUnit);

            foreach (var hm in GlobalEngineLevelSettings.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                fixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            fixedModsTreeView.DataContext = fixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalEngineLevelSettings.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                variableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            variableModsTreeView.DataContext = variableModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalEngineLevelSettings.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                localizeModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            localizeModsTreeView.DataContext = localizeModTypeForTreeViewObservableCollection;
        }

        private void UpdateFieldsFromTask(XLSearchTask task)
        {
            //Crosslink search para
            cbCrosslinker.SelectedIndex = (int)task.XlSearchParameters.CrosslinkerType;
            txtXLTopNum.Text = task.XlSearchParameters.CrosslinkSearchTopNum.ToString(CultureInfo.InvariantCulture);
            ckbSearchWithXLAllBeta.IsChecked = task.XlSearchParameters.CrosslinkSearchWithAllBeta;
            txtUdXLKerName.Text = task.XlSearchParameters.UdXLkerName;
            ckbUdXLkerCleavable.IsChecked = task.XlSearchParameters.UdXLkerCleavable;
            txtUdXLkerTotalMs.Text = task.XlSearchParameters.UdXLkerTotalMass.HasValue ? task.XlSearchParameters.UdXLkerTotalMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerShortMass.Text = task.XlSearchParameters.UdXLkerShortMass.HasValue ? task.XlSearchParameters.UdXLkerShortMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerLongMass.Text = task.XlSearchParameters.UdXLkerLongMass.HasValue ? task.XlSearchParameters.UdXLkerLongMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerAminoAcid.Text = task.XlSearchParameters.UdXLkerResidue.ToString();
            cbbXLprecusorMsTl.SelectedIndex = task.XlSearchParameters.XlPrecusorMsTl is AbsoluteTolerance ? 0 : 1;
            txtXLPrecusorMsTl.Text = task.XlSearchParameters.XlPrecusorMsTl.Value.ToString(CultureInfo.InvariantCulture);
            cbbXLBetaprecusorMsTl.SelectedIndex = task.XlSearchParameters.XlPrecusorMsTl is AbsoluteTolerance ? 0 : 1;
            txtXLBetaPrecusorMsTl.Text = task.XlSearchParameters.XlPrecusorMsTl.Value.ToString(CultureInfo.InvariantCulture);

            checkBoxDecoy.IsChecked = task.XlSearchParameters.SearchDecoy;
            int temp = task.CommonParameters.MaxMissedCleavages.Value;
            missedCleavagesTextBox.Text = temp.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParameters.MinPeptideLength.HasValue ? task.CommonParameters.MinPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtMaxPeptideLength.Text = task.CommonParameters.MaxPeptideLength.HasValue ? task.CommonParameters.MaxPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            proteaseComboBox.SelectedItem = task.CommonParameters.Protease;
            temp = task.CommonParameters.MaxModificationIsoforms.Value;
            maxModificationIsoformsTextBox.Text = temp.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = fixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    fixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = variableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    variableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }

            localizeAllCheckBox.IsChecked = task.CommonParameters.LocalizeAll;
            if (task.CommonParameters.LocalizeAll.Value)
            {
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                    heh.Use = true;
            }
            else
            {
                foreach (var mod in task.CommonParameters.ListOfModsLocalize)
                {
                    var theModType = localizeModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                    if (theModType != null)
                    {
                        var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                        if (theMod != null)
                            theMod.Use = true;
                        else
                            theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                    else
                    {
                        theModType = new ModTypeForTreeView(mod.Item1, true);
                        localizeModTypeForTreeViewObservableCollection.Add(theModType);
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
            }

            foreach (var ye in variableModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in localizeModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.XlSearchParameters.CrosslinkSearchTopNum = int.Parse(txtXLTopNum.Text, CultureInfo.InvariantCulture);
            TheTask.XlSearchParameters.CrosslinkSearchWithAllBeta = ckbSearchWithXLAllBeta.IsChecked.Value;
            TheTask.XlSearchParameters.CrosslinkerType = (CrosslinkerType)cbCrosslinker.SelectedIndex;
            if (cbbXLprecusorMsTl.SelectedIndex == 0)
                TheTask.XlSearchParameters.XlPrecusorMsTl = new AbsoluteTolerance(double.Parse(txtXLPrecusorMsTl.Text, CultureInfo.InvariantCulture));
            else
                TheTask.XlSearchParameters.XlPrecusorMsTl = new PpmTolerance(double.Parse(txtXLPrecusorMsTl.Text, CultureInfo.InvariantCulture));

            if (cbbXLBetaprecusorMsTl.SelectedIndex == 0)
                TheTask.XlSearchParameters.XlBetaPrecusorMsTl = new AbsoluteTolerance(double.Parse(txtXLBetaPrecusorMsTl.Text, CultureInfo.InvariantCulture));
            else
                TheTask.XlSearchParameters.XlBetaPrecusorMsTl = new PpmTolerance(double.Parse(txtXLBetaPrecusorMsTl.Text, CultureInfo.InvariantCulture));

            if (TheTask.XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                TheTask.XlSearchParameters.UdXLkerName = txtUdXLKerName.Text;
                TheTask.XlSearchParameters.UdXLkerCleavable = ckbUdXLkerCleavable.IsChecked.Value;
                TheTask.XlSearchParameters.UdXLkerLongMass = double.Parse(txtUdXLkerLongMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerShortMass = double.Parse(txtUdXLkerShortMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerTotalMass = double.Parse(txtUdXLkerTotalMs.Text, CultureInfo.InvariantCulture);
            }

            TheTask.XlSearchParameters.SearchDecoy = checkBoxDecoy.IsChecked.Value;
            TheTask.CommonParameters.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.MinPeptideLength = int.TryParse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out int temp) ? (int?)temp : null;
            TheTask.CommonParameters.MaxPeptideLength = int.TryParse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out temp) ? (int?)temp : null;
            TheTask.CommonParameters.Protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.CommonParameters.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            if (productMassToleranceComboBox.SelectedIndex == 0)
                TheTask.CommonParameters.ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            else
                TheTask.CommonParameters.ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            TheTask.CommonParameters.BIons = bCheckBox.IsChecked.Value;
            TheTask.CommonParameters.YIons = yCheckBox.IsChecked.Value;
            TheTask.CommonParameters.CIons = cCheckBox.IsChecked.Value;
            TheTask.CommonParameters.ZdotIons = zdotCheckBox.IsChecked.Value;

            //TheTask.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;

            TheTask.CommonParameters.ListOfModsVariable = new List<Tuple<string, string>>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
                TheTask.CommonParameters.ListOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => new Tuple<string, string>(b.Parent.DisplayName, b.DisplayName)));
            TheTask.CommonParameters.ListOfModsFixed = new List<Tuple<string, string>>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
                TheTask.CommonParameters.ListOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => new Tuple<string, string>(b.Parent.DisplayName, b.DisplayName)));
            if (localizeAllCheckBox.IsChecked.Value)
            {
                TheTask.CommonParameters.ListOfModsLocalize = null;
                TheTask.CommonParameters.LocalizeAll = true;
            }
            else
            {
                TheTask.CommonParameters.LocalizeAll = false;
                TheTask.CommonParameters.ListOfModsLocalize = new List<Tuple<string, string>>();
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                    TheTask.CommonParameters.ListOfModsLocalize.AddRange(heh.Children.Where(b => b.Use).Select(b => new Tuple<string, string>(b.Parent.DisplayName, b.DisplayName)));
            }

            DialogResult = true;
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            dataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
            //dataContextForSearchTaskWindow.ModExpanderTitle =
            //    "fixed: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
            //    + " variable: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
            //    + " localize: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName));
            dataContextForSearchTaskWindow.AnalysisExpanderTitle = "Some analysis properties...";
            dataContextForSearchTaskWindow.SearchModeExpanderTitle = "Some search properties...";
        }

        private void ModExpander_Expanded(object sender, RoutedEventArgs e)
        {
        }

        private void ModificationsDataGrid_Loaded(object sender, RoutedEventArgs e)
        {
        }

        private void ModificationsDataGrid_DataContextChanged(object sender, DependencyPropertyChangedEventArgs e)
        {
        }

        private void ModificationsDataGrid_AutoGeneratedColumns(object sender, EventArgs e)
        {
            //if (!TheTask.WritePrunedDatabase)
            //    modificationsDataGrid.Columns[3].Visibility = Visibility.Collapsed;
        }

        #endregion Private Methods
    }
}