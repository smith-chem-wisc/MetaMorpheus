using EngineLayer;
using GuiFunctions.Models;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GuiFunctions.Util;

/// <summary>
/// Manages loading and saving of custom M-Ion losses with caching
/// </summary>
public static class CustomMIonLossManager
{
    private static readonly string CustomMIonLossFileName = "CustomMIonLosses.txt";
    private static List<CustomMIonLoss> _cachedLosses;
    private static DateTime _lastFileWriteTime;

    /// <summary>
    /// Gets the full path to the custom M-Ion losses file
    /// </summary>
    public static string GetCustomMIonLossFilePath()
    {
        string modsDirectory = Path.Combine(GlobalVariables.DataDir, "Mods");
        if (!Directory.Exists(modsDirectory))
        {
            Directory.CreateDirectory(modsDirectory);
        }
        return Path.Combine(modsDirectory, CustomMIonLossFileName);
    }

    /// <summary>
    /// Clears the cache, forcing a reload on next access
    /// </summary>
    public static void ClearCache()
    {
        _cachedLosses = null;
    }

    /// <summary>
    /// Loads all custom M-Ion losses from the file with caching
    /// </summary>
    public static List<CustomMIonLoss> LoadCustomMIonLosses()
    {
        string filePath = GetCustomMIonLossFilePath();

        // Check if cache is valid
        var currentFileWriteTime = File.GetLastWriteTime(filePath);
        if (_cachedLosses != null && _lastFileWriteTime == currentFileWriteTime)
        {
            return [.._cachedLosses]; // Return a copy to prevent modification
        }

        // Load from file
        var customLosses = new List<CustomMIonLoss>();
        try
        {
            var lines = File.ReadAllLines(filePath);
            foreach (var line in lines)
            {
                // Skip empty lines and comments
                if (string.IsNullOrWhiteSpace(line) || line.TrimStart().StartsWith("#"))
                {
                    continue;
                }

                try
                {
                    var customLoss = CustomMIonLoss.Deserialize(line);
                    customLosses.Add(customLoss);
                }
                catch (Exception ex)
                {
                    // Log the error but continue loading other losses
                    System.Diagnostics.Debug.WriteLine($"Error loading custom M-Ion loss: {ex.Message}");
                }
            }

            // Update cache
            _cachedLosses = customLosses;
            _lastFileWriteTime = currentFileWriteTime;
        }
        catch (Exception ex)
        {
            System.Diagnostics.Debug.WriteLine($"Error reading custom M-Ion losses file: {ex.Message}");
        }

        return [..customLosses];
    }

    /// <summary>
    /// Saves all custom M-Ion losses to the file and updates cache
    /// </summary>
    public static void SaveCustomMIonLosses(IEnumerable<CustomMIonLoss> losses)
    {
        string filePath = GetCustomMIonLossFilePath();

        try
        {
            var lines = new List<string>
            {
                "# Custom M-Ion Losses for MetaMorpheus",
                "# Format: Name|Annotation|ChemicalFormula|AnalyteType",
                "# AnalyteType can be: Peptide, Oligo, or Proteoform",
                "#"
            };

            lines.AddRange(losses.Select(loss => loss.Serialize()));

            File.WriteAllLines(filePath, lines);
            
            // Update cache
            _cachedLosses = new List<CustomMIonLoss>(losses);
            _lastFileWriteTime = File.GetLastWriteTime(filePath);
        }
        catch (Exception ex)
        {
            throw new IOException($"Failed to save custom M-Ion losses: {ex.Message}", ex);
        }
    }

    /// <summary>
    /// Adds a new custom M-Ion loss to the file and updates cache
    /// </summary>
    public static void AddCustomMIonLoss(CustomMIonLoss loss)
    {
        var existingLosses = LoadCustomMIonLosses();
        
        // Check for duplicates
        if (existingLosses.Any(l => l.Name.Equals(loss.Name, StringComparison.OrdinalIgnoreCase) 
                                    && l.ApplicableAnalyteType == loss.ApplicableAnalyteType))
        {
            throw new InvalidOperationException($"A custom M-Ion loss with the name '{loss.Name}' already exists for {loss.ApplicableAnalyteType} mode.");
        }

        existingLosses.Add(loss);
        SaveCustomMIonLosses(existingLosses);
    }

    /// <summary>
    /// Gets all M-Ion losses (built-in + custom) for the specified analyte type
    /// </summary>
    public static List<MIonLoss> GetAllMIonLossesForAnalyteType(AnalyteType analyteType)
    {
        var allLosses = new List<MIonLoss>();

        // Add built-in losses from MIonLoss.AllMIonLosses
        if (MIonLoss.AllMIonLosses != null)
        {
            allLosses.AddRange(MIonLoss.AllMIonLosses.Values);
        }

        // Add applicable custom losses
        var customLosses = LoadCustomMIonLosses()
            .Where(l => l.IsApplicableToCurrentMode(analyteType))
            .Select(l => l.ToMIonLoss());
        
        allLosses.AddRange(customLosses);

        return allLosses;
    }

    /// <summary>
    /// Deletes a custom M-Ion loss from the file and updates cache
    /// </summary>
    public static void DeleteCustomMIonLoss(string name, AnalyteType analyteType)
    {
        var existingLosses = LoadCustomMIonLosses();
        var lossToRemove = existingLosses.FirstOrDefault(l => 
            l.Name.Equals(name, StringComparison.OrdinalIgnoreCase) 
            && l.ApplicableAnalyteType == analyteType);

        if (lossToRemove != null)
        {
            existingLosses.Remove(lossToRemove);
            SaveCustomMIonLosses(existingLosses);
        }
    }
}
