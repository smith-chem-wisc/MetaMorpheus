using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibSpectralAveraging;
using SpectralAveraging;

namespace MetaMorpheusGUI
{
    public class DataStructureToStringConverter : BaseValueConverter<DataStructureToStringConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            string val = value.ToString();

            // spectra file processing specific
            if (value.GetType() == typeof(SpectraFileProcessingType))
            {
                if (val == "AverageEverynScans")
                    return "Average Every n Scans";
                else if (val == "AverageEverynScansWithOverlap")
                    return "Average Every n Scans with Overlap";
                else if (val == "AverageDDAScans")
                    return "Average DDA Scans";
                else if (val == "AverageDDAScansWithOverlap")
                    return "Average DDA Scans with Overlap";
            }

            // output type specific
            if (value.GetType() == typeof(OutputType))
            {
                if (val == "txt")
                    return "Text File";
                else if (val == "mzML")
                    return "MzML file";
            }


            // empty
            if (string.IsNullOrWhiteSpace(val))
                    return "";

            // split all others at the capital letter
            StringBuilder newText = new StringBuilder(val.Length * 2);
            newText.Append(val[0]);
            for (int i = 1; i < val.Length; i++)
            {
                if (char.IsUpper(val[i]) && val[i - 1] != ' ')
                    newText.Append(' ');
                newText.Append(val[i]);
            }
            return new String(newText.ToString());

        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
