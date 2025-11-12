using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Data;
using GuiFunctions;

namespace MetaMorpheusGUI
{
    public class ProteinRnaTerminologyConverter : IValueConverter
    {
        // Add all relevant mappings here (case-insensitive, will match whole words)
        private static readonly Dictionary<string, string> ProteinToRna = new()
        {
            { "protease", "rnase" },
            { "Protease", "Rnase" },
            { "protein", "transcript" },
            { "proteins", "transcripts" },
            { "peptide", "oligo" },
            { "peptides", "oligos" },
            { "Protein", "Transcript" },
            { "Proteins", "Transcripts" },
            { "Peptide", "Oligo" },
            { "Peptides", "Oligos" },
            { "proteome", "transcriptome" },
            { "proteomic", "transcriptomic" },
            { "Proteoform", "Oligo" },
            { "post-translational", "post-transcriptional" },
            { "PTM", "PTrM" },
            { "PSM", "OSM" },
            { "PSMs", "OSMs" },
            { "amino acid", "nucleic acid" },
            { "ProteolyticDigestion", "Digestion" },
            { "proteases.tsv", "rnases.tsv" },
            { "N-Terminal", "5'-Terminal" },
            { "N-Terminus", "5'-Terminus" },
            { "C-Terminal", "3'-Terminal" },
            { "C-Terminus", "3'-Terminus" },
            // Add more as needed
        };

        public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            string input = parameter as string ?? value?.ToString() ?? "";

            // Try to get the ViewModel, but handle if it's not initialized yet
            try
            {
                if (GuiGlobalParamsViewModel.Instance?.IsRnaMode ?? false)
                {
                    // Replace all protein terms with their RNA counterparts
                    foreach (var kvp in ProteinToRna)
                    {
                        // Use word boundaries to avoid partial replacements
                        input = Regex.Replace(input, $@"\b{Regex.Escape(kvp.Key)}\b", kvp.Value);
                    }
                }
            }
            catch
            {
                // If instance creation fails, just return the input unchanged
                // This can happen during design time before everything is initialized
            }

            return input;
        }

        public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
