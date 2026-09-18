using System;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;

namespace EngineLayer
{
    /// <summary>
    /// The helper for seeding a data file that the USER is expected to edit.
    ///
    /// <para>
    /// Every such file -- custom proteases, rnases, crosslinkers, modifications, RNA modifications,
    /// monosaccharides, amino acids -- follows the same contract, and it exists because breaking it
    /// produced a data-loss bug (#2752: a user's edited custom monosaccharides were overwritten on
    /// install, repair and upgrade). The rules are:
    /// </para>
    ///
    /// <list type="number">
    ///   <item><description>
    ///     <b>Never touch a file that exists.</b> Seeding is guarded by <c>!File.Exists</c>. An
    ///     existing custom file is never rewritten, reformatted, migrated or validated on startup.
    ///     This is the rule that actually protects the user's work; everything else supports it.
    ///   </description></item>
    ///   <item><description>
    ///     <b>Seed a documented template, not an empty file and not a copy of the data.</b> Whatever
    ///     shape it takes it carries no data rows, so the user opens a file that explains its own format
    ///     and has nothing to delete first. Where the file has a shipped sibling that opens with a banner
    ///     and a header row, <see cref="BannerAndHeaderFrom(Stream, string)"/> derives the template from
    ///     it, which keeps the two from drifting. Where it does not -- a glycan database is a bare list,
    ///     with no header row to stop at and no bannered sibling to derive from -- the template is a
    ///     hand-written banner, read whole with <see cref="EmbeddedText(Assembly, string)"/>. It is still
    ///     a file in the repository rather than a string literal, so it is reviewed like one.
    ///   </description></item>
    ///   <item><description>
    ///     <b>The custom file is never installer- or build-managed.</b> It must not appear in
    ///     <c>Product.wxs</c> as a <c>&lt;File&gt;</c> and must not carry a
    ///     <c>&lt;None Update ... CopyToOutputDirectory&gt;</c> rule. Those two are what let an
    ///     installer or a build overwrite the user's copy, and both were removed by the #2752 fix.
    ///     The shipped file a template is read from may be either; only the custom file may not.
    ///   </description></item>
    ///   <item><description>
    ///     <b>Failing to seed is reported, not swallowed.</b> The user gets a
    ///     <see cref="MetaMorpheusException"/> naming the file, rather than a silently absent one.
    ///   </description></item>
    /// </list>
    ///
    /// <para>
    /// Not every file is seeded through this class, or from a shipped sibling:
    /// </para>
    /// <list type="bullet">
    ///   <item><description>
    ///     <c>CustomAminoAcids.txt</c> is seeded by <c>GlobalVariables.WriteAminoAcidsFile</c> with a
    ///     full A-Z dump of the existing residues rather than a bare header, because its whole purpose is
    ///     letting a user adjust the mass of a residue that already exists.
    ///   </description></item>
    ///   <item><description>
    ///     <c>MonosaccharidesCustom.tsv</c> is seeded by
    ///     <c>GlycanDatabase.EnsureCustomMonosaccharideFileExists</c>, which also carries over a legacy
    ///     copy, from its own embedded, documented template.
    ///   </description></item>
    ///   <item><description>
    ///     <c>CustomModifications.txt</c> and <c>RnaCustomModifications.txt</c> use a hand-written banner
    ///     (<c>GlobalVariables.CustomModificationsTemplate</c>), since <c>Mods.txt</c> has no header row to
    ///     cut at. <c>CustomCrosslinkers.tsv</c> gets a header alone, because the shipped
    ///     <c>Crosslinkers.tsv</c> has no banner.
    ///   </description></item>
    /// </list>
    /// </summary>
    public static class CustomDataFile
    {
        /// <summary>
        /// Creates <paramref name="path"/> from <paramref name="buildTemplate"/> if, and only if, it does
        /// not already exist. An existing file -- which is to say, one the user may have edited -- is left
        /// exactly as it is.
        /// </summary>
        /// <param name="path">Full path to the custom file.</param>
        /// <param name="buildTemplate">
        /// Produces the template contents. Deferred rather than passed as a string so that the cost of
        /// reading an embedded resource is not paid on every startup, which is the common case where the
        /// file is already there.
        /// </param>
        /// <param name="description">
        /// What the file is for, in words, used in the exception message when seeding fails.
        /// </param>
        public static void EnsureExists(string path, Func<string> buildTemplate, string description)
        {
            if (File.Exists(path))
            {
                return;
            }

            try
            {
                string directory = Path.GetDirectoryName(path);
                if (!string.IsNullOrEmpty(directory))
                {
                    Directory.CreateDirectory(directory);
                }

                // A blank file would be protected by rule 1 on every later startup and never repaired.
                string template = buildTemplate();
                if (string.IsNullOrWhiteSpace(template))
                {
                    throw new MetaMorpheusException("the template was empty");
                }

                // Written beside the destination and moved into place. File.WriteAllText truncates its
                // destination before it streams, so writing straight to `path` means an interrupted write
                // -- a full disk is the realistic one -- leaves a partial file that rule 1 then protects on
                // every later startup, and nothing ever repairs it. The move is the only step that creates
                // the real path, and it is a rename within one folder.
                string partialPath = path + ".tmp";
                try
                {
                    File.WriteAllText(partialPath, template);
                    File.Move(partialPath, path);
                }
                catch
                {
                    TryDelete(partialPath);
                    throw;
                }
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException(
                    $"Error creating the default {description} file at {path}: {e.Message}", e);
            }
        }

        /// <summary>
        /// Best-effort removal of the half-written file left by a failed seed. It is already the unhappy
        /// path, so a failure to clean up must not replace the error that got us here.
        /// </summary>
        private static void TryDelete(string path)
        {
            try
            {
                if (File.Exists(path))
                {
                    File.Delete(path);
                }
            }
            catch
            {
                // nothing useful to do, and the caller's exception is the one worth reporting
            }
        }

        /// <summary>
        /// The comment banner and header row of a shipped tab-separated file, with every data row
        /// dropped -- the template a user should be handed for the custom counterpart.
        /// </summary>
        /// <param name="shipped">The shipped file. Disposed by this method.</param>
        /// <param name="headerPrefix">
        /// How the header row starts, e.g. <c>"Name\t"</c>. Everything up to and including the first line
        /// that starts with this (ignoring leading whitespace) is kept, with the comment and blank lines
        /// above it.
        /// </param>
        /// <exception cref="MetaMorpheusException">A data row comes before any header row.</exception>
        public static string BannerAndHeaderFrom(Stream shipped, string headerPrefix)
        {
            var template = new StringBuilder();

            using (var reader = new StreamReader(shipped))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    string trimmed = line.TrimStart();

                    // a blank line inside the banner is part of the banner, not the end of it
                    if (trimmed.Length == 0 || trimmed.StartsWith("#", StringComparison.Ordinal))
                    {
                        template.AppendLine(line);
                        continue;
                    }

                    if (trimmed.StartsWith(headerPrefix, StringComparison.Ordinal))
                    {
                        template.AppendLine(line);
                        return template.ToString();
                    }

                    // a data row before any header: a template built from here would be headerless
                    break;
                }
            }

            throw new MetaMorpheusException(
                $"No header row starting with '{headerPrefix.Trim()}' was found before the data rows of the shipped file.");
        }

        /// <summary>
        /// <see cref="BannerAndHeaderFrom(Stream, string)"/> over an embedded resource.
        /// </summary>
        public static string BannerAndHeaderFrom(Assembly assembly, string resourceName, string headerPrefix)
        {
            Stream stream = assembly.GetManifestResourceStream(resourceName)
                ?? throw new MetaMorpheusException(
                    $"Embedded resource '{resourceName}' was not found in {assembly.GetName().Name}.");

            return BannerAndHeaderFrom(stream, headerPrefix);
        }

        /// <summary>
        /// The whole of an embedded template, verbatim -- for a file whose format has no header row to
        /// stop at, so <see cref="BannerAndHeaderFrom(Assembly, string, string)"/> has nothing to derive.
        /// A glycan database is one: it is a bare list of glycans, so its template is a hand-written
        /// banner of comment lines and the rules of rule 2 are met by the banner alone.
        /// </summary>
        /// <remarks>
        /// The template is still an embedded resource rather than a string literal here, so that what the
        /// user is handed is reviewable as a file and cannot drift from the format it documents.
        /// </remarks>
        public static string EmbeddedText(Assembly assembly, string resourceName)
        {
            Stream stream = assembly.GetManifestResourceStream(resourceName)
                ?? throw new MetaMorpheusException(
                    $"Embedded resource '{resourceName}' was not found in {assembly.GetName().Name}.");

            using (var reader = new StreamReader(stream))
            {
                return reader.ReadToEnd();
            }
        }

        /// <summary>
        /// <see cref="BannerAndHeaderFrom(Stream, string)"/> over a shipped file on disk.
        /// </summary>
        public static string BannerAndHeaderFromFile(string shippedPath, string headerPrefix)
        {
            return BannerAndHeaderFrom(File.OpenRead(shippedPath), headerPrefix);
        }
    }
}
