#nullable enable
using System;
using System.Buffers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;

namespace EngineLayer.DatabaseLoading;

public class DbForTask
{
    public DbForTask(string filePath, bool isContaminant, string? decoyIdentifier = null)
    {
        FilePath = filePath;
        IsContaminant = isContaminant;
        FileName = System.IO.Path.GetFileName(filePath);
        var ext = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant();
        IsSpectralLibrary = ext == ".msp" || ext == ".msl";
        DecoyIdentifier = decoyIdentifier ?? GlobalVariables.DecoyIdentifier;
    }

    public bool IsSpectralLibrary { get; }
    public string FilePath { get; }
    public bool IsContaminant { get; }
    public string FileName { get; }
    public string DecoyIdentifier { get; }
    public int? BioPolymerCount { get; internal set; } = null;
    public int? TargetCount { get; internal set; } = null;
    public int? DecoyCount { get; internal set; } = null;

    /// <summary>
    /// Lower-case hex SHA-256 of the file as given (the compressed bytes for a .gz), taken when the database was
    /// loaded or written. Null until then, or when the file could not be read (see <see cref="IdentityError"/>).
    /// A file name is not an identity: the same name is reused across UniProt releases, review states and
    /// GPTMD runs, so this is what says exactly which database a task searched.
    /// </summary>
    public string? Sha256 { get; private set; }

    /// <summary>The file's size in bytes when <see cref="Sha256"/> was taken.</summary>
    public long? SizeBytes { get; private set; }

    /// <summary>Why <see cref="Sha256"/> could not be taken, or null.</summary>
    public string? IdentityError { get; private set; }

    /// <summary>
    /// The databases this one was written from (GPTMD's output names its inputs), each with the identity it had
    /// when it was read. Empty for a database the user supplied.
    /// </summary>
    public IReadOnlyList<DbForTask> DerivedFrom { get; init; } = Array.Empty<DbForTask>();

    /// <summary>
    /// Takes the file's size and SHA-256 now. Called every time a task loads the database, because the same path
    /// can hold different content from one task (or run) to the next. A file that cannot be read is recorded as
    /// such rather than failing the task: loading it reports that problem in its own terms.
    /// The file is read in 1 MiB sequential blocks, which on a 1 GB UniProt XML in the OS cache hashes in about
    /// 0.75 s against 1.4 s with File.OpenRead's 4 KiB buffer. The share mode is the database loaders' own
    /// (read, others may read), so the parse can read the file at the same time.
    /// The block is borrowed from the shared array pool and the stream does no buffering of its own: a new 1 MiB
    /// buffer per call lands on the large object heap, and ParallelSearch loads (so hashes) one database per search,
    /// tens of thousands in a run.
    /// </summary>
    public void RecordIdentity()
    {
        byte[] block = ArrayPool<byte>.Shared.Rent(IdentityBlockSize);
        try
        {
            using var stream = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.Read, bufferSize: 0, FileOptions.SequentialScan);
            SizeBytes = stream.Length;
            using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
            int read;
            while ((read = stream.Read(block, 0, IdentityBlockSize)) > 0)
                hash.AppendData(block, 0, read);
            Sha256 = Convert.ToHexString(hash.GetHashAndReset()).ToLowerInvariant();
            IdentityError = null;
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            Sha256 = null;
            SizeBytes = null;
            IdentityError = e.Message;
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(block);
        }
    }

    private const int IdentityBlockSize = 1 << 20;

    /// <summary>The identity as one line of text: "SHA-256 &lt;hex&gt;, &lt;n&gt; bytes", or why it is unavailable.</summary>
    public string IdentityText() => Sha256 is null
        ? "SHA-256 unavailable" + (IdentityError is null ? "" : " (" + IdentityError + ")")
        : $"SHA-256 {Sha256}, {SizeBytes} bytes";

    /// <summary>"Derived from: a.xml (SHA-256 …); b.xml (SHA-256 …)", or null for a database the user supplied.</summary>
    public string? DerivedFromText() => DerivedFrom.Count == 0
        ? null
        : "Derived from: " + string.Join("; ", DerivedFrom.Select(p => p.FileName + " (" + p.IdentityText() + ")"));
}
