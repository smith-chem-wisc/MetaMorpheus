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
    ///     <b>Seed a documented template, not an empty file and not a copy of the data.</b> The
    ///     template is the shipped sibling's comment banner plus its header row, and no data rows --
    ///     so the user opens a file that explains its own format and has nothing to delete first.
    ///     Where there is a shipped sibling, <see cref="BannerAndHeaderFrom(Stream, string)"/> derives
    ///     exactly that from it, which keeps the two from drifting.
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

                File.WriteAllText(path, template);
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
        /// <see cref="BannerAndHeaderFrom(Stream, string)"/> over a shipped file on disk.
        /// </summary>
        public static string BannerAndHeaderFromFile(string shippedPath, string headerPrefix)
        {
            return BannerAndHeaderFrom(File.OpenRead(shippedPath), headerPrefix);
        }
    }
}
