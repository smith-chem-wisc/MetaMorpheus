using EngineLayer;
using MzLibUtil;
using Nett;
using Proteomics.ProteolyticDigestion;
using MassSpectrometry;

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
                    case nameof(Protease):
                        Protease = keyValuePair.Value.Get<Protease>(); break;
                    case nameof(MinPeptideLength):
                        MinPeptideLength = keyValuePair.Value.Get<int>(); break;
                    case nameof(MaxPeptideLength):
                        MaxPeptideLength = keyValuePair.Value.Get<int>(); break;
                    case nameof(MaxMissedCleavages):
                        MaxMissedCleavages = keyValuePair.Value.Get<int>(); break;
                    case nameof(MaxModsForPeptide):
                        MaxModsForPeptide = keyValuePair.Value.Get<int>(); break;


                    case nameof(DissociationType):
                        DissociationType = keyValuePair.Value.Get<MassSpectrometry.DissociationType>(); break;


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
        public Protease Protease { get; set; }
        public int? MinPeptideLength { get; set; }
        public int? MaxPeptideLength { get; set; }
        public int? MaxMissedCleavages { get; set; }
        public int? MaxModsForPeptide { get; set; }
        public DissociationType? DissociationType { get; set; }

        // This method is to make sure developers keep consistent naming between CommonParameters and FileSpecificParameters.
        // It's supposed to immediately crash MetaMorpheus if you rename a Common Parameter and don't rename it here.
        // The reason this method exists is to make sure toml settings are written and parsed consistently between the tasks
        // and the file-specific settings.
        public static void ValidateFileSpecificVariableNames()
        {
            CommonParameters temp = new CommonParameters();

            if (!nameof(temp.PrecursorMassTolerance).Equals(nameof(PrecursorMassTolerance)))
                throw new MetaMorpheusException("Precursor tol variable name is inconsistent");
            if (!nameof(temp.ProductMassTolerance).Equals(nameof(ProductMassTolerance)))
                throw new MetaMorpheusException("Product tol variable name is inconsistent");
            if (!nameof(temp.DigestionParams.Protease).Equals(nameof(Protease)))
                throw new MetaMorpheusException("Protease variable name is inconsistent");
            if (!nameof(temp.DigestionParams.MinPeptideLength).Equals(nameof(MinPeptideLength)))
                throw new MetaMorpheusException("Min peptide length variable name is inconsistent");
            if (!nameof(temp.DigestionParams.MaxPeptideLength).Equals(nameof(MaxPeptideLength)))
                throw new MetaMorpheusException("Max peptide length variable name is inconsistent");
            if (!nameof(temp.DigestionParams.MaxMissedCleavages).Equals(nameof(MaxMissedCleavages)))
                throw new MetaMorpheusException("Max missed cleavages variable name is inconsistent");
            if (!nameof(temp.DigestionParams.MaxModsForPeptide).Equals(nameof(MaxModsForPeptide)))
                throw new MetaMorpheusException("Max mods per peptide variable name is inconsistent");

        }

        public FileSpecificParameters Clone()
        {
            return (FileSpecificParameters)this.MemberwiseClone();
        }
    }
}