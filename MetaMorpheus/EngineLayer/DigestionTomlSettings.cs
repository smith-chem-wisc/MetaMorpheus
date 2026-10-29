using Nett;
using Proteomics.ProteolyticDigestion;
using Transcriptomics.Digestion;

namespace EngineLayer
{
    /// <summary>
    /// How MetaMorpheus writes and reads digestion parameters in TOML: the protease and RNase conversions, and the
    /// derived properties a saved task leaves out. This is the only definition of it.
    /// </summary>
    /// <remarks>
    /// Two places need the same text. <c>MetaMorpheusTask.tomlConfig</c> saves and reloads task settings with it, and
    /// the index cache key (<c>IndexingEngine.ToString</c>) writes the digestion parameters with it, so every setting a
    /// saved task carries is in the key. TaskLayer references EngineLayer and not the reverse, so the definition lives
    /// here and <c>tomlConfig</c> is built on it. A new digestion parameters type, digestion agent or derived property
    /// is configured here, once, and both get it.
    ///
    /// <c>tomlConfig</c> keeps one digestion entry of its own: the read for <see cref="Omics.Digestion.IDigestionParams"/>,
    /// which picks the concrete type from the table and translates removed protease names, warning through the task.
    /// It only reads, and the key only writes.
    /// </remarks>
    public static class DigestionTomlSettings
    {
        /// <summary>
        /// Adds the digestion types to a TOML settings builder.
        /// </summary>
        public static TomlSettings.ITomlSettingsBuilder ConfigureDigestionTypes(this TomlSettings.ITomlSettingsBuilder cfg) => cfg
            .ConfigureType<Protease>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.ToString())
                    // GetProtease, not the raw indexer: it maps the names older releases wrote, such as
                    // "chymotrypsin (don't cleave before proline)", onto the enzyme they meant.
                    .FromToml(tmlString => ProteaseDictionary.GetProtease(tmlString.Value))))
            .ConfigureType<DigestionParams>(type => type
                .IgnoreProperty(p => p.DigestionAgent)
                .IgnoreProperty(p => p.SpecificDigestionAgent)
                .IgnoreProperty(p => p.MaxMods)
                .IgnoreProperty(p => p.MaxLength)
                .IgnoreProperty(p => p.MinLength))
            .ConfigureType<RnaDigestionParams>(type => type
                .IgnoreProperty(p => p.DigestionAgent)
                .IgnoreProperty(p => p.SpecificDigestionAgent))
            .ConfigureType<Rnase>(type => type
                .WithConversionFor<TomlString>(convert => convert
                    .ToToml(custom => custom.Name)
                    .FromToml(tmlString => RnaseDictionary.Dictionary[tmlString.Value])));

        /// <summary>
        /// The digestion types and nothing else, for writing digestion parameters on their own.
        /// </summary>
        public static readonly TomlSettings Settings = TomlSettings.Create(cfg => cfg.ConfigureDigestionTypes());
    }
}
