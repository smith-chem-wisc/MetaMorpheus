using Omics;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
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
    public string TaskId { get; } = taskId;
    public bool GenerateTargets { get; } = generateTargets;
    public DecoyType DecoyType { get; } = decoyType;
    public TargetContaminantAmbiguity TcAmbiguity { get; } = tcAmbiguity;
    public List<string> LocalizableMods { get; } = localizableMods ?? [];
    public List<DbForTask> DbForTask { get; } = dbFilenameList;

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        Status("Loading Databases...");

        var bioPolymers = LoadBioPolymers(TaskId, DbForTask, GenerateTargets, DecoyType, LocalizableMods, CommonParameters);
        SanitizeBioPolymerDatabase(bioPolymers, TcAmbiguity);

        Status("Done.");
        var results = new DatabaseLoadingEngineResults(this, DbForTask, bioPolymers);
        return results;
    }

    public List<IBioPolymer> LoadBioPolymers(string taskId, List<DbForTask> dbFilenameList, bool searchTarget, DecoyType decoyType, List<string> localizeableModificationTypes, CommonParameters commonParameters)
    {
        Status($"Loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s...");
        int emptyEntries = 0;
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
            Warn($"Warning: No {GlobalVariables.AnalyteType.GetBioPolymerLabel()} entries were found in the database");
        }
        else if (emptyEntries > 0)
        {
            Warn("Warning: " + emptyEntries + $" empty {GlobalVariables.AnalyteType.GetBioPolymerLabel()} entries ignored");
        }

        // We are not generating decoys, so just return the read in database
        if (!bioPolymerList.Any(p => p.IsDecoy))
        {
            Status($"Done loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s");
            return bioPolymerList;
        }

        // Sanitize the decoys
        // TODO: Fix this so that it accounts for multi-protease searches. Currently, we only consider the first protease
        // when looking for target/decoy collisions
        HashSet<string> targetPeptideSequences = new();
        foreach (var bioPolymer in bioPolymerList.Where(p => !p.IsDecoy))
        {
            // When thinking about decoy collisions, we can ignore modifications
            foreach (var peptide in bioPolymer.Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>()))
            {
                targetPeptideSequences.Add(peptide.BaseSequence);
            }
        }
        // Now, we iterate through the decoys and scramble the sequences that correspond to target peptides
        for (int i = 0; i < bioPolymerList.Count; i++)
        {
            if (bioPolymerList[i].IsDecoy)
            {
                var peptidesToReplace = bioPolymerList[i]
                    .Digest(commonParameters.DigestionParams, new List<Modification>(), new List<Modification>())
                    .Select(p => p.BaseSequence)
                    .Where(targetPeptideSequences.Contains)
                    .ToList();
                if (peptidesToReplace.Any())
                {
                    bioPolymerList[i] = DecoySequenceValidator.ScrambleDecoyBioPolymer(bioPolymerList[i], commonParameters.DigestionParams, forbiddenSequences: targetPeptideSequences, peptidesToReplace);
                }
            }
        }

        Status($"Done loading {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s");
        return bioPolymerList;
    }

    public static IEnumerable<RNA> LoadOligoDb(string fileName, bool generateTargets, DecoyType decoyType,
        List<string> localizeableModificationTypes, bool isContaminant,
        out Dictionary<string, Modification> unknownMods, out int emptyEntriesCount,
        CommonParameters commonParameters, string decoyIdentifier = "DECOY")
    {
        List<string> dbErrors = new List<string>();
        List<RNA> rnaList = new List<RNA>();

        string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
        bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
        theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

        if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
        {
            unknownMods = null;
            rnaList = RnaDbLoader.LoadRnaFasta(fileName, generateTargets, decoyType, isContaminant, out dbErrors);
        }
        else
        {
            // TODO: Add in variant params when fixed in MzLib. 
            List<string> modTypesToExclude = GlobalVariables.AllRnaModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
            rnaList = RnaDbLoader.LoadRnaXML(fileName, generateTargets, decoyType, isContaminant, GlobalVariables.AllRnaModsKnown, modTypesToExclude, out unknownMods, commonParameters.MaxThreadsToUsePerFile, decoyIdentifier: decoyIdentifier);
        }

        emptyEntriesCount = rnaList.Count(p => p.BaseSequence.Length == 0);
        return rnaList.Where(p => p.BaseSequence.Length > 0);
    }

    public static IEnumerable<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, out Dictionary<string, Modification> um,
            out int emptyEntriesCount, CommonParameters commonParameters, string decoyIdentifier = "DECOY")
    {
        List<string> dbErrors = new List<string>();
        List<Protein> proteinList = new List<Protein>();

        string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
        bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
        theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

        if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
        {
            um = null;
            proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, out dbErrors,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, commonParameters.MaxThreadsToUsePerFile, addTruncations: commonParameters.AddTruncations);
        }
        else
        {
            List<string> modTypesToExclude = GlobalVariables.AllModTypesKnown.Where(b => !localizeableModificationTypes.Contains(b)).ToList();
            //proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, commonParameters.MaxHeterozygousVariants, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations, decoyIdentifier: decoyIdentifier);
            proteinList = ProteinDbLoader.LoadProteinXML(fileName, generateTargets, decoyType, GlobalVariables.AllModsKnown, isContaminant, modTypesToExclude, out um, commonParameters.MaxThreadsToUsePerFile, 0, commonParameters.MinVariantDepth, addTruncations: commonParameters.AddTruncations, decoyIdentifier: decoyIdentifier);
        }

        emptyEntriesCount = proteinList.Count(p => p.BaseSequence.Length == 0);
        return proteinList.Where(p => p.BaseSequence.Length > 0);
    }

    public void SanitizeBioPolymerDatabase<TBioPolymer>(List<TBioPolymer> bioPolymers, TargetContaminantAmbiguity tcAmbiguity)
            where TBioPolymer : IBioPolymer
    {
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
                Warn("The protein '" + accession + "' has multiple entries. Protein accessions must be unique. Protein " + accession + " was renamed.");
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
                    Warn($"The {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s '" + accession + $"' has multiple entries. {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s accessions must be unique. Contaminant {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s " + accession + " was ignored.");
                }
            }
            else if (tcAmbiguity == TargetContaminantAmbiguity.RemoveTarget)
            {
                // remove targets
                toRemove.AddRange(accessionGroup.Where(p => !p.IsDecoy && !p.IsContaminant));
                if (toRemove.Any())
                {
                    Warn($"The {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s '" + accession + $"' has multiple entries. {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s accessions must be unique. Target {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s " + accession + " was ignored.");
                }
            }

            //remove the bioPolymers specified above
            foreach (var bioPolymer in toRemove.Where(_ => accessionGroup.Count > 1))
            {
                bioPolymers.RemoveAll(p => ReferenceEquals(p, bioPolymer));
                accessionGroup.RemoveAll(p => ReferenceEquals(p, bioPolymer));
            }

            // most ambiguity should be handled by now, but for edge cases and decoys:
            // remove bioPolymers so that only 1 bioPolymer with this accession remains
            for (int i = 0; i < accessionGroup.Count - 1; i++) //-1 to keep the last one (most mods)
            {
                bioPolymers.RemoveAll(p => ReferenceEquals(p, accessionGroup[i]));
            }
        }
    }
}

public class DatabaseLoadingEngineResults(DatabaseLoadingEngine engine, List<DbForTask> loadedDbs, List<IBioPolymer> bioPolymers)
    : MetaMorpheusEngineResults(engine)
{
    private List<DbForTask> loadedDbs = loadedDbs;
    public List<IBioPolymer> BioPolymers { get; init; } = bioPolymers;

    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine(base.ToString());
        sb.AppendLine("Databases loaded:");
        foreach (var db in loadedDbs)
            sb.AppendLine($"\t{Path.GetFileName(db.FilePath)}: {db.BioPolymerCount} {GlobalVariables.AnalyteType.GetBioPolymerLabel()}s ({db.TargetCount} target, {db.DecoyCount} decoy)");
        return sb.ToString();
    }
}

