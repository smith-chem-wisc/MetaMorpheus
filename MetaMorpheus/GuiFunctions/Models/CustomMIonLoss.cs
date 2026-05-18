using Chemistry;
using EngineLayer;
using Omics.Fragmentation;
using System;

namespace GuiFunctions.Models;

/// <summary>
/// Represents a custom M-Ion loss that can be defined by the user.
/// Extends MIonLoss with AnalyteType information for filtering between Peptide and Oligo modes.
/// </summary>
public class CustomMIonLoss : MIonLoss
{
    /// <summary>
    /// Indicates which analyte type(s) this loss applies to
    /// </summary>
    public AnalyteType ApplicableAnalyteType { get; set; }

    public CustomMIonLoss(string name, string annotation, ChemicalFormula chemicalFormula, AnalyteType applicableAnalyteType) 
        : base(name, annotation, chemicalFormula)
    {
        ApplicableAnalyteType = applicableAnalyteType;
    }

    /// <summary>
    /// Converts a CustomMIonLoss to a regular MIonLoss for use in fragmentation
    /// </summary>
    public MIonLoss ToMIonLoss()
    {
        return new MIonLoss(Name, Annotation, ThisChemicalFormula);
    }

    /// <summary>
    /// Creates a CustomMIonLoss from a MIonLoss and AnalyteType
    /// </summary>
    public static CustomMIonLoss FromMIonLoss(MIonLoss mIonLoss, AnalyteType analyteType)
    {
        return new CustomMIonLoss(mIonLoss.Name, mIonLoss.Annotation, mIonLoss.ThisChemicalFormula, analyteType);
    }

    /// <summary>
    /// Serializes the custom M-Ion loss to a string format for file storage
    /// Format: Name|Annotation|ChemicalFormula|AnalyteType
    /// </summary>
    public string Serialize()
    {
        return $"{Name}|{Annotation}|{ThisChemicalFormula.Formula}|{ApplicableAnalyteType}";
    }

    /// <summary>
    /// Deserializes a custom M-Ion loss from a string
    /// </summary>
    public static CustomMIonLoss Deserialize(string line)
    {
        var parts = line.Split('|');
        if (parts.Length != 4)
        {
            throw new FormatException($"Invalid custom M-Ion loss format: {line}");
        }

        string name = parts[0];
        string annotation = parts[1];
        string formulaString = parts[2];
        
        if (!Enum.TryParse<AnalyteType>(parts[3], out var analyteType))
        {
            throw new FormatException($"Invalid AnalyteType in custom M-Ion loss: {parts[3]}");
        }

        ChemicalFormula formula;
        try
        {
            formula = ChemicalFormula.ParseFormula(formulaString);
        }
        catch (Exception ex)
        {
            throw new FormatException($"Invalid chemical formula '{formulaString}' in custom M-Ion loss", ex);
        }

        return new CustomMIonLoss(name, annotation, formula, analyteType);
    }

    /// <summary>
    /// Checks if this loss is applicable to the current analyte type
    /// </summary>
    public bool IsApplicableToCurrentMode(AnalyteType? currentAnalyteType = null)
    {
        switch (ApplicableAnalyteType)
        {
            // If the loss is for both types, it's always applicable
            case AnalyteType.Peptide when currentAnalyteType == AnalyteType.Peptide:
            case AnalyteType.Peptide when currentAnalyteType == null && !GuiGlobalParamsViewModel.Instance.IsRnaMode:
            case AnalyteType.Oligo when currentAnalyteType == AnalyteType.Oligo:
            case AnalyteType.Oligo when currentAnalyteType == null && GuiGlobalParamsViewModel.Instance.IsRnaMode:
            case AnalyteType.Proteoform when currentAnalyteType == AnalyteType.Proteoform:
            case AnalyteType.Proteoform when currentAnalyteType == null && !GuiGlobalParamsViewModel.Instance.IsRnaMode:
                return true;
            default:
                return false;
        }
    }
}
