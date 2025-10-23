#nullable enable
using MzLibUtil;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Omics.Digestion;
using Transcriptomics;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

namespace EngineLayer.DatabaseLoading;
public class DatabaseLoadingEngine(
    CommonParameters commonParameters,
    List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
    List<string> nestedIds,
    List<DbForTask> dbFilenameList,
    string taskId,
    DecoyType decoyType,
    bool generateTargets = true,
    List<string> localizableMods = null,
    TargetContaminantAmbiguity tcAmbiguity = TargetContaminantAmbiguity.RemoveContaminant)
    : MetaMorpheusEngine(commonParameters, fileSpecificParameters, nestedIds)
{
    private static readonly ListPool<string> _sequencePool = new();

    public string TaskId { get; } = taskId;
    public bool GenerateTargets { get; } = generateTargets;
    public DecoyType DecoyType { get; } = decoyType;
    public TargetContaminantAmbiguity TcAmbiguity { get; } = tcAmbiguity;
    public List<string> LocalizableMods { get; } = localizableMods ?? [];
    public List<DbForTask> DbForTask { get; } = dbFilenameList;

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        int scrambled = 0;
        int removed;
        int renamed;

        Status($"Loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s...");

        var bioPolymers = LoadBioPolymers(TaskId, DbForTask, GenerateTargets, DecoyType, LocalizableMods, CommonParameters, out var errors);
        foreach (var error in errors)
            Warn(error);

        if (bioPolymers.Any(p => p.IsDecoy))
            scrambled = ScrambleHomologousDecoys(bioPolymers, CommonParameters.DigestionParams);

        (removed, renamed) = SanitizeBioPolymerDatabase(bioPolymers, TcAmbiguity, out errors);
        foreach (var error in errors)
            Warn(error);

        Status($"Done loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s");
        var results = new DatabaseLoadingEngineResults(this, DbForTask, bioPolymers, removed, renamed, scrambled);
        return results;
    }

    public static List<IBioPolymer> LoadBioPolymers(string taskId, List<DbForTask> dbFilenameList, bool searchTarget, DecoyType decoyType, List<string> localizeableModificationTypes, CommonParameters commonParameters, out List<string> errors)
    {
        int emptyEntries = 0;
        errors = new();
        List<IBioPolymer> bioPolymerList = new();
        foreach (var db in dbFilenameList.Where(p => !p.IsSpectralLibrary))
        {
            IEnumerable<IBioPolymer> dbBioPolymers;
            int targetCount = 0;
            int decoyCount = 0;
            if (GlobalVariables.AnalyteType == AnalyteType.Oligo)
            {
                dbBioPolymers = LoadOligoDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out int emptyOligoEntriesForThisDb, commonParameters, db.DecoyIdentifier);
                emptyEntries += emptyOligoEntriesForThisDb;
            }
            else
            {
                dbBioPolymers = LoadProteinDb(db.FilePath, searchTarget, decoyType, localizeableModificationTypes, db.IsContaminant, out Dictionary<string, Modification> unknownModifications, out int emptyProteinEntriesForThisDb, commonParameters, db.DecoyIdentifier);
                emptyEntries += emptyProteinEntriesForThisDb;
            }

            foreach (var bioPol in dbBioPolymers)
            {
                bioPolymerList.Add(bioPol);
                if (bioPol.IsDecoy)
                    decoyCount++;
                else
                    targetCount++;
            }

            db.BioPolymerCount = targetCount + decoyCount;
            db.TargetCount = targetCount;
            db.DecoyCount = decoyCount;
        }
        if (!bioPolymerList.Any())
        {
            errors.Add($"Warning: No {GlobalVariables.AnalyteType.GetBioPolymerLabel()} entries were found in the database");
        }
        else if (emptyEntries > 0)
        {
            errors.Add("Warning: " + emptyEntries + $" empty {GlobalVariables.AnalyteType.GetBioPolymerLabel()} entries ignored");
        }

        return bioPolymerList;
    }

    public static IEnumerable<RNA> LoadOligoDb(string fileName, bool generateTargets, DecoyType decoyType,
        List<string> localizeableModificationTypes, bool isContaminant,
        out Dictionary<string, Modification> unknownMods, out int emptyEntriesCount,
        CommonParameters commonParameters, string? decoyIdentifier = null)
    {
        decoyIdentifier ??= GlobalVariables.DecoyIdentifier;
        List<RNA> rnaList;

        string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
        bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
        theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

        if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
        {
            unknownMods = null;
            rnaList = RnaDbLoader.LoadRnaFasta(fileName, generateTargets, decoyType, isContaminant, out var dbErrors);
        }
        else
        {
            GlobalVariables.AddMods(ProteinDbLoader.GetPtmListFromProteinXml(fileName), true, true);
            // TODO: Add in variant params when fixed in MzLib. 
            List<string> modTypesToExclude = GlobalVariables.AllRnaModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
            rnaList = RnaDbLoader.LoadRnaXML(fileName, generateTargets, decoyType, isContaminant, GlobalVariables.AllRnaModsKnown, modTypesToExclude, out unknownMods, commonParameters.MaxThreadsToUsePerFile, decoyIdentifier: decoyIdentifier);
        }

        emptyEntriesCount = rnaList.Count(p => p.BaseSequence.Length == 0);
        return rnaList.Where(p => p.BaseSequence.Length > 0);
    }

    public static IEnumerable<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um,
            out int emptyEntriesCount, CommonParameters commonParameters, string? decoyIdentifier = null)
    {
        decoyIdentifier ??= GlobalVariables.DecoyIdentifier;
        List<Protein> proteinList;

        string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
        bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
        theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

        if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
        {
            um = null;
            proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, out var dbErrors,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, commonParameters.MaxThreadsToUsePerFile, addTruncations: commonParameters.AddTruncations);
        }
        else
        {
            List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
            GlobalVariables.AddMods(ProteinDbLoader.GetPtmListFromProteinXml(fileName), true, false);
            //proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, commonParameters.MaxHeterozygousVariants, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations, decoyIdentifier: decoyIdentifier);
            proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, 0, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations, decoyIdentifier: decoyIdentifier);
        }
        emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
        return proteinList.Where(p => p.BaseSequence.Length > 0);
    }

    /// <summary>
    /// handles Accession Collisions
    /// </summary>
    /// <exception cref="ArgumentException"></exception>
    public static (int Removed, int Renamed) SanitizeBioPolymerDatabase<TBioPolymer>(List<TBioPolymer> bioPolymers, TargetContaminantAmbiguity tcAmbiguity, out List<string> errors) where TBioPolymer : IBioPolymer
    {
        int removed = 0;
        int rename = 0;
        errors = new();
        List<TBioPolymer> toRemove = new();
        bioPolymers.RemoveAll(p => p == null);
        foreach (var accessionGroup in bioPolymers.GroupBy(p => p.Accession)
                     .Where(group => group.Count() > 1) // only keep the ones with multiple entries sharing an accession
                     .Select(group => group.OrderBy(p => p.OneBasedPossibleLocalizedModifications.Count) // order by mods then truncation products (this is what was here before)
                         .ThenBy(p => p.TruncationProducts.Count)
                         .ToList()) // Individual ordered accession group to list from IEnumerable
                     .ToList()) // Collapse entire group and sort enumerable to a list so we can modify the bioPolymers collection. 
        {
            toRemove.Clear();
            string accession = accessionGroup.First().Accession;

            if (tcAmbiguity == TargetContaminantAmbiguity.RenameProtein)
            {
                int bioPolymerNumber = 1;
                errors.Add("The protein '" + accession + "' has multiple entries. Protein accessions must be unique. Protein " + accession + " was renamed.");
                foreach (var originalBioPolymer in accessionGroup)
                {
                    //accession is private and there's no clone method, so we need to make a whole new bioPolymer... TODO: put this in mzlib
                    //use PROTEIN_D1 instead of PROTEIN_1 so it doesn't look like an isoform (D for Duplicate)
                    IBioPolymer renamed;
                    if (originalBioPolymer is RNA r)
                    {
                        renamed = new RNA(originalBioPolymer.BaseSequence, originalBioPolymer.Accession + "_D" + bioPolymerNumber,
                            r.OneBasedPossibleLocalizedModifications, r.FivePrimeTerminus, r.ThreePrimeTerminus, r.Name, r.Organism,
                            r.DatabaseFilePath, r.IsContaminant, r.IsDecoy, r.GeneNames, r.AdditionalDatabaseFields, r.TruncationProducts,
                            r.SequenceVariations, r.AppliedSequenceVariations, r.SampleNameForVariants, r.FullName);
                    }
                    else
                    {
                        Protein p = originalBioPolymer as Protein ?? throw new ArgumentException($"Database sanitization assumed BioPolymer was a protein when it was {originalBioPolymer.GetType()}");
                        renamed = new Protein(originalBioPolymer.BaseSequence, originalBioPolymer.Accession + "_D" + bioPolymerNumber, originalBioPolymer.Organism,
                            originalBioPolymer.GeneNames, originalBioPolymer.OneBasedPossibleLocalizedModifications, p.TruncationProducts, originalBioPolymer.Name, originalBioPolymer.FullName,
                            originalBioPolymer.IsDecoy, originalBioPolymer.IsContaminant, p.DatabaseReferences, p.SequenceVariations, p.AppliedSequenceVariations,
                            p.SampleNameForVariants, p.DisulfideBonds, p.SpliceSites, originalBioPolymer.DatabaseFilePath);
                    }

                    bioPolymers.Add((TBioPolymer)renamed);
                    bioPolymers.RemoveAll(m => ReferenceEquals(m, originalBioPolymer));
                    bioPolymerNumber++;
                    rename++;
                }

                continue;
            }

            // if we are not renaming, we need to remove the duplicates
            if (tcAmbiguity == TargetContaminantAmbiguity.RemoveContaminant)
            {
                // remove contaminants
                toRemove.AddRange(accessionGroup.Where(p => p.IsContaminant));
                if (toRemove.Any())
                {
                    errors.Add($"The {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s '" + accession + $"' has multiple entries. {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s accessions must be unique. Contaminant {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s " + accession + " was ignored.");
                }
            }
            else if (tcAmbiguity == TargetContaminantAmbiguity.RemoveTarget)
            {
                // remove targets
                toRemove.AddRange(accessionGroup.Where(p => !p.IsDecoy && !p.IsContaminant));
                if (toRemove.Any())
                {
                    errors.Add($"The {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s '" + accession + $"' has multiple entries. {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s accessions must be unique. Target {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s " + accession + " was ignored.");
                }
            }

            //remove the bioPolymers specified above
            foreach (var bioPolymer in toRemove.Where(_ => accessionGroup.Count > 1))
            {
                bioPolymers.RemoveAll(p => ReferenceEquals(p, bioPolymer));
                accessionGroup.RemoveAll(p => ReferenceEquals(p, bioPolymer));
            }
            removed += toRemove.Count;

            // most ambiguity should be handled by now, but for edge cases and decoys:
            // remove bioPolymers so that only 1 bioPolymer with this accession remains
            for (int i = 0; i < accessionGroup.Count - 1; i++) //-1 to keep the last one (most mods)
            {
                bioPolymers.RemoveAll(p => ReferenceEquals(p, accessionGroup[i]));
                removed++;
            }
        }
        return (removed, rename);
    }

    /// <summary>
    /// Handles Decoy-Target sequence collisions. If a decoy produces any peptides that are identical to target peptides, the decoy sequence is scrambled until no collisions remain.   
    /// </summary>
    public static int ScrambleHomologousDecoys<TBioPolymer>(List<TBioPolymer> bioPolymers, IDigestionParams digestionParams) where TBioPolymer : IBioPolymer
    {
        // Sanitize the decoys
        // TODO: Fix this so that it accounts for multi-protease searches. Currently, we only consider the first protease
        // when looking for target/decoy collisions
        HashSet<string> targetPeptideSequences = new();
        foreach (var bioPolymer in bioPolymers.Where(p => !p.IsDecoy))
        {
            // When thinking about decoy collisions, we can ignore modifications
            foreach (var peptide in bioPolymer.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
            {
                targetPeptideSequences.Add(peptide.BaseSequence);
            }
        }

        var variable = new List<Modification>();
        var fixedM = new List<Modification>();
        int scrambled = 0;

        // Now, we iterate through the decoys and scramble the sequences that correspond to target peptides
        for (int i = 0; i < bioPolymers.Count; i++)
        {
            var bioPolymer = bioPolymers[i];
            if (!bioPolymer.IsDecoy)
                continue;

            // Use a pooled list to avoid repeated allocations in tight loops
            List<string> peptidesToReplace = null;
            foreach (var peptide in bioPolymer.Digest(digestionParams, variable, fixedM))
            {
                if (targetPeptideSequences.Contains(peptide.BaseSequence))
                {
                    peptidesToReplace ??= _sequencePool.Get();
                    peptidesToReplace.Add(peptide.BaseSequence);
                }
            }

            if (peptidesToReplace is { Count: > 0 })
            {
                bioPolymers[i] = DecoySequenceValidator.ScrambleDecoyBioPolymer(
                    bioPolymer,
                    digestionParams,
                    forbiddenSequences: targetPeptideSequences,
                    peptidesToReplace);
                _sequencePool.Return(peptidesToReplace);
                scrambled++;
            }
        }
        return scrambled;
    }
}

public class DatabaseLoadingEngineResults(DatabaseLoadingEngine engine, List<DbForTask> loadedDbs, List<IBioPolymer> bioPolymers, int removedBySanitization, int renamedBySanitization, int scrambled)
    : MetaMorpheusEngineResults(engine)
{
    public List<IBioPolymer> BioPolymers { get; init; } = bioPolymers;
    public int BioPolymersRemovedBySanitization { get; init; } = removedBySanitization;
    public int BioPolymersRenamedBySanitization { get; init; } = renamedBySanitization;
    public int DecoysScrambledDueToHomology { get; init; } = scrambled;

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        sb.AppendLine("Databases loaded:");
        foreach (var db in loadedDbs)
        {
            sb.AppendLine(
                $"\t{Path.GetFileName(db.FilePath)}: {db.BioPolymerCount} {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s");
            sb.AppendLine($"\t\t{db.TargetCount} Target {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s \t {db.DecoyCount} Decoy {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s");
        }
        sb.AppendLine($"Total {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s loaded: {BioPolymers.Count + BioPolymersRemovedBySanitization}");
        sb.AppendLine($"Total {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s removed by sanitization: {BioPolymersRemovedBySanitization}");
        sb.AppendLine($"Total {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s renamed by sanitization: {BioPolymersRenamedBySanitization}");
        sb.AppendLine($"Total decoy {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s scrambled due to homology: {DecoysScrambledDueToHomology}");
        sb.AppendLine($"Total {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s remaining: {BioPolymers.Count}");
        return sb.ToString();
    }
}

