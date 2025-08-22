using System;
using System.Globalization;

namespace MetaMorpheusGUI;

public class FilePathToNameConverter : BaseValueConverter<FilePathToNameConverter>
{
    public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
    {
        if (value is string path)
        {
            return System.IO.Path.GetFileName(path);
        }
        return value;
    }

    public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
    {
        // Not supported: cannot convert a file name back to a full path
        throw new NotSupportedException("ConvertBack is not supported for FilePathToNameConverter.");
    }
}