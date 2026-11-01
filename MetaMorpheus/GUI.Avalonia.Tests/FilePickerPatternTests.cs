using System.Linq;
using Avalonia.Headless.NUnit;
using EngineLayer;
using MetaMorpheusAvalonia.ViewModels;
using NUnit.Framework;

namespace Test.AvaloniaGui;

/// <summary>
/// The file-picker patterns carry two platform invariants that a comment alone cannot enforce.
/// They live on the view model rather than on MainWindow precisely so they can be tested - the
/// window itself is [ExcludeFromCodeCoverage], for the same reason Test.csproj does not reference
/// the WPF GUI project.
/// </summary>
public class FilePickerPatternTests
{
    /// <summary>
    /// Lowercase is not cosmetic. Avalonia's FreeDesktop backend hands these to the XDG portal as
    /// GlobStyle globs, and portal matching is case-sensitive, so "*.mzML" hides sample.mzml on
    /// Linux - the platform this GUI exists for, and the one platform where running the app on a
    /// developer's machine would never reveal it, because Windows and macOS match case-insensitively.
    /// </summary>
    [AvaloniaTest]
    public void SpectraPatternsAreLowercase()
    {
        string[] patterns = MainWindowViewModel.SpectraPatterns;

        Assert.That(patterns, Is.Not.Empty);
        foreach (string pattern in patterns)
        {
            Assert.That(pattern, Is.EqualTo(pattern.ToLowerInvariant()),
                "an uppercase glob silently hides matching files under the XDG portal");
        }
    }

    /// <summary>
    /// Bruker data is a directory and OpenFilePickerAsync cannot select one, so offering "*.d" would
    /// give the user a filter that matches nothing. They pick the .tdf inside it instead, and
    /// AddSpectraFiles maps that back to the folder.
    /// </summary>
    [AvaloniaTest]
    public void SpectraPatternsExcludeBrukerDirectoriesButNothingElse()
    {
        Assert.That(GlobalVariables.AcceptedSpectraFormats, Contains.Item(".d"),
            "if .d stops being an accepted format, this exclusion needs revisiting");

        Assert.That(MainWindowViewModel.SpectraPatterns, Does.Not.Contain("*.d"));

        // every other accepted format must still be offered, or the picker quietly hides valid input
        foreach (string extension in GlobalVariables.AcceptedSpectraFormats.Where(e => e != ".d"))
        {
            Assert.That(MainWindowViewModel.SpectraPatterns, Contains.Item("*" + extension.ToLowerInvariant()));
        }
    }

    /// <summary>Databases are routinely distributed gzipped, so both forms have to be offered.</summary>
    [AvaloniaTest]
    public void DatabasePatternsOfferPlainAndGzippedForms()
    {
        string[] patterns = MainWindowViewModel.DatabasePatterns;

        foreach (string extension in GlobalVariables.AcceptedDatabaseFormats)
        {
            Assert.That(patterns, Contains.Item("*" + extension));
            Assert.That(patterns, Contains.Item("*" + extension + ".gz"));
        }
    }
}
