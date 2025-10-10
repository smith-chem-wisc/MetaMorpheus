using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Proteomics;
using UsefulProteomicsDatabases;
using Omics.BioPolymer;
using Omics.Modifications;
using EngineLayer;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using TaskLayer;

public static class VariantDbDiagnostics
{
    /// <summary>
    /// Loads a protein XML database, simulates variant isoform generation, and logs any SequenceVariation errors.
    /// </summary>
    /// <param name="xmlDbPath">Path to the protein XML database.</param>
    /// <param name="logPath">Path to write the error log.</param>
    public static void DiagnoseProteinXmlDbRealistic(string xmlDbPath, string logPath)
    {
        var errors = new List<string>();
        Dictionary<string, Modification> unknownMods;
        int emptyEntriesCount;
        CommonParameters parameters = new CommonParameters
        {

        };
        try
        {
            // Use the same parameters as MetaMorpheus default
            var proteins = ProteinDbLoader.LoadProteinXML(
                xmlDbPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: GlobalVariables.AllModsKnown, // or as used in your search
                isContaminant: false,
                modTypesToExclude: new List<string>(), // or as used in your search
                unknownModifications: out unknownMods,
                maxThreads: parameters.MaxThreadsToUsePerFile,
                maxSequenceVariantsPerIsoform: parameters.MaxSequenceVariantsPerIsoform,
                minAlleleDepth: parameters.MinAlleleDepth,
                maxSequenceVariantIsoforms: parameters.MaxSequenceVariantIsoforms,
                addTruncations: parameters.AddTruncations
            );

            int proteinIndex = 0;
            foreach (var protein in proteins)
            {
                proteinIndex++;
                // Simulate variant isoform generation as in a real search
                List<Protein> isoforms;
                try
                {
                    isoforms = protein.GetVariantBioPolymers(
                        minAlleleDepth: 0,
                        maxSequenceVariantsPerIsoform: 0,
                        maxSequenceVariantIsoforms: 1
                    ).OfType<Protein>().ToList();
                }
                catch (Exception ex)
                {
                    errors.Add($"Protein Accession: {protein.Accession}, Index: {proteinIndex}, Exception during GetVariantBioPolymers: {ex.Message}");
                    continue;
                }

                int isoformIndex = 0;
                foreach (var isoform in isoforms)
                {
                    isoformIndex++;
                    // Check SequenceVariations
                    if (isoform.SequenceVariations != null)
                    {
                        foreach (var sv in isoform.SequenceVariations)
                        {
                            if (!sv.AreValid())
                            {
                                errors.Add(
                                    $"Protein Accession: {protein.Accession}, Isoform: {isoformIndex}, SequenceVariation: [{sv.OneBasedBeginPosition}-{sv.OneBasedEndPosition}] {sv.OriginalSequence}->{sv.VariantSequence}, Description: {sv.Description}, VCF: {sv.VariantCallFormatData?.Description ?? "null"} (SequenceVariations)");
                            }
                        }
                    }
                    // Check AppliedSequenceVariations
                    if (isoform.AppliedSequenceVariations != null)
                    {
                        foreach (var sv in isoform.AppliedSequenceVariations)
                        {
                            if (!sv.AreValid())
                            {
                                errors.Add(
                                    $"Protein Accession: {protein.Accession}, Isoform: {isoformIndex}, SequenceVariation: [{sv.OneBasedBeginPosition}-{sv.OneBasedEndPosition}] {sv.OriginalSequence}->{sv.VariantSequence}, Description: {sv.Description}, VCF: {sv.VariantCallFormatData?.Description ?? "null"} (AppliedSequenceVariations)");
                            }
                        }
                    }
                }
            }
        }
        catch (Exception ex)
        {
            errors.Add($"Exception during database load or isoform generation: {ex.Message}\n{ex.StackTrace}");
        }

        if (errors.Count == 0)
        {
            File.WriteAllText(logPath, "No SequenceVariation errors detected during realistic isoform generation.");
            Console.WriteLine("No SequenceVariation errors detected during realistic isoform generation.");
        }
        else
        {
            File.WriteAllLines(logPath, errors);
            Console.WriteLine($"Logged {errors.Count} SequenceVariation errors to {logPath}");
        }
    }

    [TestFixture]
    public class VariantDbDiagnosticsTests
    {
        [Test]
        public void DiagnoseProteinXmlDbRealistic_FindsNoInvalidSequenceVariations()
        {
            // Set your test database and log file paths here
            string xmlDbPath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_07.xml";
            string logPath = @"E:\Projects\Mann_11cell_lines\A549\A549_1\log.txt";

            // Run the improved diagnostic
            VariantDbDiagnostics.DiagnoseProteinXmlDbRealistic(xmlDbPath, logPath);

            // Read the log and assert
            string logContents = File.ReadAllText(logPath);
            Assert.That(logContents.Contains("No SequenceVariation errors detected"),
                $"SequenceVariation errors were found. See log at {logPath}:\n{logContents}");
        }
    }
}

