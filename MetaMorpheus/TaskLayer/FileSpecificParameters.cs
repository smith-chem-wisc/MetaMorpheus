using System.Collections.Generic;
using EngineLayer;
using MzLibUtil;
using Nett;
using Proteomics.ProteolyticDigestion;
using MassSpectrometry;
using Omics.Digestion;
using Omics.Fragmentation;
using Transcriptomics.Digestion;

namespace TaskLayer
{
    public class FileSpecificParameters
    {
        public FileSpecificParameters(TomlTable tomlTable)
        {
            foreach (var keyValuePair in tomlTable)
            {
                switch (keyValuePair.Key)
                {
                    // we're using the name of the variable here and not a fixed string
                    // in case the variable name changes at some point
                    case nameof(PrecursorMassTolerance):
                        PrecursorMassTolerance = keyValuePair.Value.Get<Tolerance>(); break;
                    case nameof(ProductMassTolerance):
                        ProductMassTolerance = keyValuePair.Value.Get<Tolerance>(); break;
                    case nameof(DigestionAgent): // Support new tomls that labeled by Digestion Agent Type instead of specific type
                        string valueString = keyValuePair.Value.Get<string>();

                        // If type is top-down we only have the analyte type to go off of. 
                        // This should be updated if there is ever an Rnase and Protease that share a same name. 
                        if (valueString.Contains("top-down"))
                        {
                            if (GlobalVariables.AnalyteType == AnalyteType.Oligo)
                                DigestionAgent = keyValuePair.Value.Get<Rnase>();
                            else
                                DigestionAgent = keyValuePair.Value.Get<Protease>();
                        }

                        // If type is not top-down we can check the digestion agent dictionaries. 
                        if (ProteaseDictionary.Dictionary.ContainsKey(valueString))
                            DigestionAgent = keyValuePair.Value.Get<Protease>();
                        else if (RnaseDictionary.Dictionary.ContainsKey(valueString))
                            DigestionAgent = keyValuePair.Value.Get<Rnase>();
                        else
                            throw new MetaMorpheusException("Unrecognized digestion agent type \"" + valueString + "\" in file-specific parameters toml");
                        break;
                    case "Rnase": // Support old tomls that labeled by Digestion Agent Type instead of DigestionAgent
                        DigestionAgent = keyValuePair.Value.Get<Rnase>();
                        break;
                    case "Protease": // Support old tomls that labeled by Digestion Agent Type instead of DigestionAgent
                        DigestionAgent = keyValuePair.Value.Get<Protease>();
                        break;
                    case nameof(MinPeptideLength):
                        MinPeptideLength = keyValuePair.Value.Get<int>(); break;
                    case nameof(MaxPeptideLength):
                        MaxPeptideLength = keyValuePair.Value.Get<int>(); break;
                    case nameof(MaxMissedCleavages):
                        MaxMissedCleavages = keyValuePair.Value.Get<int>(); break;
                    case nameof(MaxModsForPeptide):
                        MaxModsForPeptide = keyValuePair.Value.Get<int>(); break;
                    case nameof(SeparationType):
                        SeparationType = keyValuePair.Value.Get<string>(); break;


                    case nameof(DissociationType):
                        DissociationType = keyValuePair.Value.Get<DissociationType>(); break;
                    case nameof(CustomIons):
                        CustomIons = keyValuePair.Value.Get<List<ProductType>>(); break;


                    default:
                        throw new MetaMorpheusException("Unrecognized parameter \"" + keyValuePair.Key + "\" in file-specific parameters toml");
                }
            }
        }

        public FileSpecificParameters()
        {
            // everything initialized to null
        }

        public Tolerance PrecursorMassTolerance { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public DigestionAgent DigestionAgent { get; set; }
        public int? MinPeptideLength { get; set; }
        public int? MaxPeptideLength { get; set; }
        public int? MaxMissedCleavages { get; set; }
        public int? MaxModsForPeptide { get; set; }
        public DissociationType? DissociationType { get; set; }
        public string SeparationType { get; set; }

        public List<ProductType> CustomIons { get; set; }

        public FileSpecificParameters Clone()
        {
            return (FileSpecificParameters)this.MemberwiseClone();
        }
    }
}