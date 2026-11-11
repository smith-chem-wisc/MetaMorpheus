#nullable enable
using Readers;

namespace TaskLayer.ParallelSearch.IO;

public abstract class ParallelSearchResultFile<TResult> : ResultFile<TResult>, IResultFile
{
    public override SupportedFileType FileType { get; }
    public override Software Software { get; set; } = Software.MetaMorpheus;

    protected ParallelSearchResultFile(string filePath) : base(filePath) { }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    protected ParallelSearchResultFile() : base() { }

    protected static string EscapeCsv(string? value)
    {
        if (string.IsNullOrEmpty(value))
            return string.Empty;

        // If the value contains comma, quote, or newline, wrap in quotes and escape internal quotes
        if (value.Contains(',') || value.Contains('"') || value.Contains('\n'))
        {
            return $"\"{value.Replace("\"", "\"\"")}\"";
        }

        return value;
    }
}