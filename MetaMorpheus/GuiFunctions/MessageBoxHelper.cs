using System.Windows;

namespace GuiFunctions;

/// <summary>
/// Provides helper methods for displaying message boxes in a WPF application,
/// with support for suppressing message boxes globally. Suppression is typically for testing methods. 
/// </summary>
public static class MessageBoxHelper
{
    /// <summary>
    /// Gets or sets a value indicating whether message boxes should be suppressed.
    /// When set to <c>true</c>, no message boxes will be shown.
    /// </summary>
    public static bool SuppressMessageBoxes { get; set; } = false;

    /// <summary>
    /// Displays a standard message box with the specified message, unless message boxes are suppressed.
    /// </summary>
    /// <param name="message">The message to display in the message box.</param>
    public static void Show(string message)
    {
        if (!SuppressMessageBoxes)
            MessageBox.Show(message);
    }

    /// <summary>
    /// Displays a warning message box with the specified message, unless message boxes are suppressed.
    /// </summary>
    /// <param name="message">The warning message to display in the message box.</param>
    public static void Warn(string message)
    {
        if (!SuppressMessageBoxes)
            MessageBox.Show(message, "Warning", MessageBoxButton.OK, MessageBoxImage.Warning);
    }
}
