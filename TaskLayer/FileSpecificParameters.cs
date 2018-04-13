﻿using EngineLayer;
using MzLibUtil;
using Nett;

namespace TaskLayer
{
    public class FileSpecificParameters
    {
        #region Public Constructors

        public FileSpecificParameters(TomlTable tomlTable)
        {
            foreach(var keyValuePair in tomlTable)
            {
                switch(keyValuePair.Key)
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
                    case nameof(BIons):
                        BIons = keyValuePair.Value.Get<bool>(); break;
                    case nameof(YIons):
                        YIons = keyValuePair.Value.Get<bool>(); break;
                    case nameof(CIons):
                        CIons = keyValuePair.Value.Get<bool>(); break;
                    case nameof(ZdotIons):
                        ZdotIons = keyValuePair.Value.Get<bool>(); break;
                    default:
                        throw new MetaMorpheusException("Unrecognized parameter \"" + keyValuePair.Key + "\" in file-specific parameters toml");
                }
            }
        }

        public FileSpecificParameters()
        {
            // everything initialized to null
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public Protease Protease { get; set; }
        public int? MinPeptideLength { get; set; }
        public int? MaxPeptideLength { get; set; }
        public int? MaxMissedCleavages { get; set; }
        public int? MaxModsForPeptide { get; set; }
        public bool? BIons { get; set; }
        public bool? YIons { get; set; }
        public bool? CIons { get; set; }
        public bool? ZdotIons { get; set; }

        #endregion Public Properties

        #region Public Methods

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
            if (!nameof(temp.BIons).Equals(nameof(BIons)))
                throw new MetaMorpheusException("B ion variable name is inconsistent");
            if (!nameof(temp.YIons).Equals(nameof(YIons)))
                throw new MetaMorpheusException("Y ion variable name is inconsistent");
            if (!nameof(temp.ZdotIons).Equals(nameof(ZdotIons)))
                throw new MetaMorpheusException("Zdot ion variable name is inconsistent");
            if (!nameof(temp.CIons).Equals(nameof(CIons)))
                throw new MetaMorpheusException("C ion variable name is inconsistent");
        }

        public FileSpecificParameters Clone()
        {
            return (FileSpecificParameters)this.MemberwiseClone();
        }

        #endregion Public Methods
    }
}