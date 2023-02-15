using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using Nett;
using SpectralAveraging;

namespace GuiFunctions
{
    public class SpectralAveragingParametersViewModel : BaseViewModel
    {
        #region Private Members

        private SpectralAveragingParameters spectralAveragingParameters;
        private int previousOverlap;
        private string? savedPath;

        #endregion

        #region Public Properties

        public SpectralAveragingParameters SpectralAveragingParameters
        {
            get { return spectralAveragingParameters; }
            set { spectralAveragingParameters = value; OnPropertyChanged(nameof(SpectralAveragingParameters)); }
        }

        public OutlierRejectionType RejectionType
        {
            get { return spectralAveragingParameters.OutlierRejectionType; }
            set { spectralAveragingParameters.OutlierRejectionType = value; OnPropertyChanged(nameof(RejectionType)); }
        }

        public SpectraWeightingType WeightingType
        {
            get { return spectralAveragingParameters.SpectralWeightingType; }
            set { spectralAveragingParameters.SpectralWeightingType = value; OnPropertyChanged(nameof(WeightingType)); }
        }

        public NormalizationType NormalizationType
        {
            get { return spectralAveragingParameters.NormalizationType; }
            set { spectralAveragingParameters.NormalizationType = value; OnPropertyChanged(nameof(NormalizationType)); }
        }

        public double Percentile
        {
            get { return spectralAveragingParameters.Percentile; }
            set { spectralAveragingParameters.Percentile = value; OnPropertyChanged(nameof(Percentile)); }
        }

        public double MinSigmaVale
        {
            get { return spectralAveragingParameters.MinSigmaValue; }
            set { spectralAveragingParameters.MinSigmaValue = value; OnPropertyChanged(nameof(MinSigmaVale)); }
        }

        public double MaxSigmaValue
        {
            get { return spectralAveragingParameters.MaxSigmaValue; }
            set { spectralAveragingParameters.MaxSigmaValue = value; OnPropertyChanged(nameof(MaxSigmaValue)); }
        }

        public double BinSize
        {
            get { return spectralAveragingParameters.BinSize; }
            set { spectralAveragingParameters.BinSize = value; OnPropertyChanged(nameof(BinSize)); }
        }

        public SpectraFileAveragingType SpectraFileProcessingType
        {
            get { return spectralAveragingParameters.SpectraFileAveragingType; }
            set
            {
                spectralAveragingParameters.SpectraFileAveragingType = value;
                if ((SpectraFileAveragingType)value is SpectraFileAveragingType.AverageDdaScans
                    or SpectraFileAveragingType.AverageEverynScans)
                {
                    if (previousOverlap == 0)
                    {
                        previousOverlap = ScanOverlap;
                    }

                    ScanOverlap = 0;
                }
                else if (previousOverlap != 0)
                {
                    ScanOverlap = previousOverlap;
                }

                OnPropertyChanged(nameof(SpectraFileProcessingType));
            }
        }

        public int NumberOfScansToAverage
        {
            get { return spectralAveragingParameters.NumberOfScansToAverage; }
            set { spectralAveragingParameters.NumberOfScansToAverage = value; OnPropertyChanged(nameof(NumberOfScansToAverage)); }
        }

        public int ScanOverlap
        {
            get { return spectralAveragingParameters.ScanOverlap; }
            set
            {
                if (value >= NumberOfScansToAverage)
                    MessageBox.Show("Overlap cannot be greater than or equal to the number of spectra averaged");
                else
                {
                    spectralAveragingParameters.ScanOverlap = value;
                    OnPropertyChanged(nameof(ScanOverlap));
                }
            }
        }

        public string Name
        {
            get
            {
                if (savedPath != null)
                    return Path.GetFileNameWithoutExtension(savedPath);
                else
                    return "Default Options";
            }
        }

        public string SavedPath
        {
            get => savedPath;
            set { savedPath = value; OnPropertyChanged(nameof(SavedPath)); }
        }

        public OutputType OutputType
        {
            get { return spectralAveragingParameters.OutputType; }
            set { spectralAveragingParameters.OutputType = value; OnPropertyChanged(nameof(OutputType)); }
        }

        public OutlierRejectionType[] RejectionTypes { get; set; }
        public SpectraWeightingType[] WeightingTypes { get; set; }
        public SpectraFileAveragingType[] SpectraFileProcessingTypes { get; set; }
        public OutputType[] OutputTypes { get; set; }
        public NormalizationType[] NormalizationTypes { get; set; }

        #endregion

        #region Commands

        #endregion

        #region Constructor

        public SpectralAveragingParametersViewModel(SpectralAveragingParameters parameters)
        {
            // value initialization
            spectralAveragingParameters = parameters;
            RejectionTypes = (OutlierRejectionType[])Enum.GetValues(typeof(OutlierRejectionType));
            WeightingTypes = (SpectraWeightingType[])Enum.GetValues(typeof(SpectraWeightingType));
            SpectraFileProcessingTypes = (SpectraFileAveragingType[])Enum.GetValues(typeof(SpectraFileAveragingType));
            OutputTypes = (OutputType[])Enum.GetValues(typeof(OutputType));
            NormalizationTypes = (NormalizationType[])Enum.GetValues(typeof(NormalizationType));

            // command assignment
        }

        #endregion

        #region Helpers

        public void ResetDefaults()
        {
            SpectralAveragingParameters.SetDefaultValues();
            UpdateVisualRepresentation();
        }

        public void UpdateVisualRepresentation()
        {
            OnPropertyChanged(nameof(SpectralAveragingParameters));
            OnPropertyChanged(nameof(RejectionType));
            OnPropertyChanged(nameof(WeightingType));
            OnPropertyChanged(nameof(NormalizationType));
            OnPropertyChanged(nameof(Percentile));
            OnPropertyChanged(nameof(MinSigmaVale));
            OnPropertyChanged(nameof(MaxSigmaValue));
            OnPropertyChanged(nameof(BinSize));
            OnPropertyChanged(nameof(SpectraFileProcessingType));
            OnPropertyChanged((nameof(NumberOfScansToAverage)));
            OnPropertyChanged(nameof(ScanOverlap));
            OnPropertyChanged(nameof(OutputType));
            OnPropertyChanged(nameof(Name));
        }

        #endregion

    }

    /// <summary>
    /// Model for design time viewing
    /// </summary>
    [ExcludeFromCodeCoverage]
    public class SpectralAveragingParametersModel : SpectralAveragingParametersViewModel
    {
        public static SpectralAveragingParametersModel Instance => new SpectralAveragingParametersModel();

        public SpectralAveragingParametersModel() : base(new SpectralAveragingParameters())
        {

        }
    }
}
