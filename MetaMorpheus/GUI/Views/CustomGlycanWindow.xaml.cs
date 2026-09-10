using EngineLayer;
using System;
using System.IO;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomGlycanWindow.xaml -- adds one glycan to the user's own O-glycan or
    /// N-glycan database.
    /// </summary>
    /// <remarks>
    /// The window collects and sanity-checks; <see cref="GlycanDatabase.PersistCustomGlycan"/> is what
    /// decides whether the entry is a glycan and writes it, so the format the file is read in and the format
    /// it is written in cannot drift apart. Same split as CustomMonosaccharideWindow.
    /// </remarks>
    public partial class CustomGlycanWindow : Window
    {
        public CustomGlycanWindow()
        {
            InitializeComponent();
        }

        /// <summary>Whether the O-glycan database is the one selected. Index 0 is O, index 1 is N.</summary>
        private bool IsOGlycanSelected => databaseComboBox.SelectedIndex == 0;

        private string SelectedDatabasePath => IsOGlycanSelected
            ? GlobalVariables.CustomOGlycanDatabasePath
            : GlobalVariables.CustomNGlycanDatabasePath;

        private void SaveCustomGlycan_Click(object sender, RoutedEventArgs e)
        {
            string glycanText = glycanTextBox.Text.Trim();

            if (ErrorsDetected(glycanText))
            {
                return;
            }

            string databasePath = SelectedDatabasePath;
            try
            {
                GlycanDatabase.PersistCustomGlycan(glycanText, databasePath, IsOGlycanSelected);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error saving glycan: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            MessageBox.Show(
                $"The glycan \"{glycanText}\" was added to {Path.GetFileName(databasePath)}. Select that database in a "
                + "GlycoSearch task to search for it.",
                "Success", MessageBoxButton.OK, MessageBoxImage.Information);

            DialogResult = true;
        }

        private void CancelCustomGlycan_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        /// <summary>
        /// The checks worth answering before the user has left the window. Everything else -- whether the
        /// glycan parses, whether its format matches the file, whether it is already there -- is the
        /// engine's call, and its message is shown as-is.
        /// </summary>
        private bool ErrorsDetected(string glycanText)
        {
            if (string.IsNullOrWhiteSpace(glycanText))
            {
                MessageBox.Show("Please enter a glycan.", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (glycanText.StartsWith("#", StringComparison.Ordinal))
            {
                MessageBox.Show(
                    "A line beginning with '#' is a comment and would be ignored when the database is read. Enter the glycan itself.",
                    "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            return false;
        }
    }
}
