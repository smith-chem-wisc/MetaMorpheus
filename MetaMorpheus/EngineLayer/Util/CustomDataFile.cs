using System;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;

namespace EngineLayer
{
    /// <summary>
    /// The one place MetaMorpheus creates a data file that the USER is expected to edit.
    ///
    /// <para>
    /// Every such file -- custom proteases, rnases, crosslinkers, modifications, RNA modifications,
    /// monosaccharides -- follows the same recipe, and it exists because breaking it produced a
    /// data-loss bug (#2752: a user's edited custom monosaccharides were overwritten on install,
    /// repair and upgrade). The rules are:
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
    ///     <b>The template ships inside an assembly, never as a file on disk.</b> It must not appear
    ///     in <c>Product.wxs</c> as a <c>&lt;File&gt;</c> and must not carry a
    ///     <c>&lt;None Update ... CopyToOutputDirectory&gt;</c> rule. Those two are what let an
    ///     installer or a build overwrite the user's copy, and both were removed by the #2752 fix.
    ///   </description></item>
    ///   <item><description>
    ///     <b>Failing to seed is reported, not swallowed.</b> The user gets a
    ///     <see cref="MetaMorpheusException"/> naming the file, rather than a silently absent one.
    ///   </description></item>
    /// </list>
    ///
    /// <para>
    /// One deliberate exception: <c>CustomAminoAcids.txt</c> is seeded with a full A-Z dump of the
    /// existing residues rather than a bare header, because its whole purpose is letting a user
    /// adjust the mass of a residue that already exists. A header-only template would make the
    /// common case harder, not easier. See <c>GlobalVariables.WriteAminoAcidsFile</c>.
    /// </para>
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

                File.WriteAllText(path, buildTemplate());
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException(
                    $"Error creating the default {description} file at {path}: {e.Message}", e);
            }
        }

        /// <summary>
        /// The comment banner and header row of a shipped tab-separated file, with every data row
        /// dropped -- the template a user should be handed for the custom counterpart.
        /// </summary>
        /// <param name="shipped">The shipped file. Disposed by this method.</param>
        /// <param name="headerPrefix">
        /// How the header row starts, e.g. <c>"Name\t"</c>. Everything up to and including the first line
        /// that starts with this (ignoring leading whitespace, and skipping comment lines) is kept.
        /// </param>
        public static string BannerAndHeaderFrom(Stream shipped, string headerPrefix)
        {
            var template = new StringBuilder();

            using (var reader = new StreamReader(shipped))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    bool isComment = line.StartsWith("#", StringComparison.Ordinal);
                    bool isHeader = !isComment && line.TrimStart().StartsWith(headerPrefix, StringComparison.Ordinal);

                    if (!isComment && !isHeader)
                    {
                        // a data row: the banner is over and the header was already written, or there
                        // was no header to find
                        break;
                    }

                    template.AppendLine(line);

                    if (isHeader)
                    {
                        break;
                    }
                }
            }

            return template.ToString();
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
