using System.Windows;

namespace GuiFunctions;

public static class MessageBoxHelper
{
    public static bool SuppressMessageBoxes { get; set; } = false;
    public static void Show(string message)
    {
        if (!SuppressMessageBoxes)
            MessageBox.Show(message);
    }
    
    public static void Warn(string message)
    {
        if (!SuppressMessageBoxes)
            MessageBox.Show(message, "Warning", MessageBoxButton.OK, MessageBoxImage.Warning);
    }
}
