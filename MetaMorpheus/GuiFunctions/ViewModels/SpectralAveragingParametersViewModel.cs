using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Input;
using MathNet.Numerics.Providers.LinearAlgebra;
using Nett;
using SpectralAveraging;

namespace GuiFunctions
{
    public class SpectralAveragingParametersViewModel : BaseViewModel, IEquatable<SpectralAveragingParametersViewModel>
    {
        #region Private Members

        private SpectralAveragingParameters spectralAveragingParameters;
        private int previousOverlap;
        private string? savedPath;

        #endregion

        #region Public Properties

        enum PresetAveragingParameters
        {
            Dda1,
            Dda2,
            DirectInjection
        }

        public SpectralAveragingParameters SpectralAveragingParameters
        {
            get => spectralAveragingParameters;
            set { spectralAveragingParameters = value; OnPropertyChanged(nameof(SpectralAveragingParameters)); UpdateVisualRepresentation(); }
        }

        public OutlierRejectionType RejectionType
        {
            get => spectralAveragingParameters.OutlierRejectionType;
            set { spectralAveragingParameters.OutlierRejectionType = value; OnPropertyChanged(nameof(RejectionType)); }
        }

        public SpectraWeightingType WeightingType
        {
            get => spectralAveragingParameters.SpectralWeightingType;
            set { spectralAveragingParameters.SpectralWeightingType = value; OnPropertyChanged(nameof(WeightingType)); }
        }

        public SpectraFileAveragingType SpectraFileAveragingType
        {
            get => spectralAveragingParameters.SpectraFileAveragingType;
            set { spectralAveragingParameters.SpectraFileAveragingType = value; OnPropertyChanged(nameof(SpectraFileAveragingType)); }
        }

        public bool PerformNormalization
        {
            get => spectralAveragingParameters.NormalizationType == NormalizationType.RelativeToTics ? true : false;
            set
            {
                spectralAveragingParameters.NormalizationType = value == true ? NormalizationType.RelativeToTics : NormalizationType.NoNormalization;
                OnPropertyChanged(nameof(PerformNormalization));
            }
        }

        public double Percentile
        {
            get => spectralAveragingParameters.Percentile;
            set { spectralAveragingParameters.Percentile = value; OnPropertyChanged(nameof(Percentile)); }
        }

        public double MinSigmaValue
        {
            get => spectralAveragingParameters.MinSigmaValue;
            set { spectralAveragingParameters.MinSigmaValue = value; OnPropertyChanged(nameof(MinSigmaValue)); }
        }

        public double MaxSigmaValue
        {
            get => spectralAveragingParameters.MaxSigmaValue;
            set { spectralAveragingParameters.MaxSigmaValue = value; OnPropertyChanged(nameof(MaxSigmaValue)); }
        }

        public double BinSize
        {
            get => spectralAveragingParameters.BinSize;
            set { spectralAveragingParameters.BinSize = value; OnPropertyChanged(nameof(BinSize)); }
        }

        public int NumberOfScansToAverage
        {
            get => spectralAveragingParameters.NumberOfScansToAverage;
            set 
            { 
                spectralAveragingParameters.NumberOfScansToAverage = value;
                spectralAveragingParameters.ScanOverlap = value - 1;
                OnPropertyChanged(nameof(NumberOfScansToAverage));
            }
        }

        public int MaxThreads
        {
            get => spectralAveragingParameters.MaxThreadsToUsePerFile;
            set
            {
                spectralAveragingParameters.MaxThreadsToUsePerFile = value;
                OnPropertyChanged(nameof(MaxThreads));
            }
        }

        public OutlierRejectionType[] RejectionTypes { get; set; }
        public SpectraWeightingType[] WeightingTypes { get; set; }
        public SpectraFileAveragingType[] SpectraFileAveragingTypes { get; set; }

        #endregion

        #region Constructor

        public SpectralAveragingParametersViewModel(SpectralAveragingParameters parameters)
        {
            // value initialization
            spectralAveragingParameters = parameters;
            RejectionTypes = (OutlierRejectionType[])Enum.GetValues(typeof(OutlierRejectionType));
            WeightingTypes = new [] { SpectraWeightingType.WeightEvenly, SpectraWeightingType.TicValue};
            SpectraFileAveragingTypes = new[] { SpectraFileAveragingType.AverageAll, SpectraFileAveragingType.AverageDdaScansWithOverlap, SpectraFileAveragingType.AverageEverynScansWithOverlap};
            UpdateVisualRepresentation();
        }

        #endregion

        #region Command Methods

        public ICommand SetOtherParametersCommand { get; set; }

        /// <summary>
        /// Used to set a few default/preset parameter types
        /// </summary>
        /// <param name="settingsNameToSet"></param>
        /// <exception cref="ArgumentException"></exception>
        public void SetOtherParameters(object settingsNameToSet)
        {
            var parameters = new SpectralAveragingParameters()
            {
                OutputType = OutputType.MzML,
                NormalizationType = NormalizationType.RelativeToTics,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
                BinSize = 0.01,
                SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScansWithOverlap,
        };

            switch (Enum.Parse<PresetAveragingParameters>(settingsNameToSet.ToString() ?? "Dda1"))
            {
                case (PresetAveragingParameters.Dda1):
                    parameters.NumberOfScansToAverage = 5;
                    parameters.ScanOverlap = 4;
                    parameters.MaxSigmaValue = 3;
                    parameters.MinSigmaValue = 0.5;
                    parameters.OutlierRejectionType = OutlierRejectionType.SigmaClipping;
                    break;

                case (PresetAveragingParameters.Dda2):
                    parameters.NumberOfScansToAverage = 5;
                    parameters.ScanOverlap = 4;
                    parameters.MaxSigmaValue = 3;
                    parameters.MinSigmaValue = 0.5;
                    parameters.OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping;
                    break;

                case (PresetAveragingParameters.DirectInjection):
                    parameters.NumberOfScansToAverage = 15;
                    parameters.ScanOverlap = 14;
                    parameters.OutlierRejectionType = OutlierRejectionType.MinMaxClipping;
                    break;

                default:
                    throw new ArgumentException("This should never be hit!");
            }

            spectralAveragingParameters = parameters;
            UpdateVisualRepresentation();
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
            OnPropertyChanged(nameof(SpectraFileAveragingType));
            OnPropertyChanged(nameof(PerformNormalization));
            OnPropertyChanged(nameof(Percentile));
            OnPropertyChanged(nameof(MinSigmaValue));
            OnPropertyChanged(nameof(MaxSigmaValue));
            OnPropertyChanged(nameof(BinSize));
            OnPropertyChanged((nameof(NumberOfScansToAverage)));
            OnPropertyChanged(nameof(MaxThreads));
        }

        #endregion

        public bool Equals(SpectralAveragingParametersViewModel other)
        {
            if (other is null) return false;

            // check view model
            if (RejectionType != other.RejectionType) return false;
            if (WeightingType != other.WeightingType) return false;
            if (SpectraFileAveragingType != other.SpectraFileAveragingType) return false;
            if (PerformNormalization != other.PerformNormalization) return false;
            if (Math.Abs(Percentile - other.Percentile) > 0.001) return false;
            if (Math.Abs(MinSigmaValue - other.MinSigmaValue) > 0.001) return false;
            if (Math.Abs(MaxSigmaValue - other.MaxSigmaValue) > 0.001) return false;
            if (Math.Abs(BinSize - other.BinSize) > 0.001) return false;
            if (NumberOfScansToAverage != other.NumberOfScansToAverage) return false;
            if (MaxThreads != other.MaxThreads) return false;


            // check internal parameters not attached to view model directly
            if (SpectralAveragingParameters.SpectralAveragingType != other.SpectralAveragingParameters.SpectralAveragingType) return false;
            if (SpectralAveragingParameters.OutputType != other.SpectralAveragingParameters.OutputType) return false;
            if (SpectralAveragingParameters.ScanOverlap != other.SpectralAveragingParameters.ScanOverlap) return false;
             return true;
        }
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
