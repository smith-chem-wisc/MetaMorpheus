using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using MzLibSpectralAveraging;
using Nett;
using SpectralAveraging;

namespace GuiFunctions
{
    public class MzLibAveragingOptionsViewModel : BaseViewModel
    {
        #region Private Members

        private MzLibSpectralAveragingOptions mzLibSpectralAveragingOptions;
        private int previousOverlap;
        private string? savedPath;

        #endregion

        #region Public Properties

        public MzLibSpectralAveragingOptions MzLibSpectralAveragingOptions
        {
            get { return mzLibSpectralAveragingOptions; }
            set { mzLibSpectralAveragingOptions = value; OnPropertyChanged(nameof(MzLibSpectralAveragingOptions)); }
        }

        public RejectionType RejectionType
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.RejectionType; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.RejectionType = value; OnPropertyChanged(nameof(RejectionType)); }
        }

        public WeightingType WeightingType
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.WeightingType; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.WeightingType = value; OnPropertyChanged(nameof(WeightingType)); }
        }

        public bool PerformNormalization
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.PerformNormalization; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.PerformNormalization = value; OnPropertyChanged(nameof(PerformNormalization)); }
        }

        public double Percentile
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.Percentile; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.Percentile = value; OnPropertyChanged(nameof(Percentile)); }
        }

        public double MinSigmaVale
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.MinSigmaValue; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.MinSigmaValue = value; OnPropertyChanged(nameof(MinSigmaVale)); }
        }

        public double MaxSigmaValue
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.MaxSigmaValue; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.MaxSigmaValue = value; OnPropertyChanged(nameof(MaxSigmaValue)); }
        }

        public double BinSize
        {
            get { return mzLibSpectralAveragingOptions.SpectralAveragingOptions.BinSize; }
            set { mzLibSpectralAveragingOptions.SpectralAveragingOptions.BinSize = value; OnPropertyChanged(nameof(BinSize)); }
        }

        public SpectraFileProcessingType SpectraFileProcessingType
        {
            get { return mzLibSpectralAveragingOptions.SpectraFileProcessingType; }
            set
            {
                mzLibSpectralAveragingOptions.SpectraFileProcessingType = value;
                if ((SpectraFileProcessingType)value is SpectraFileProcessingType.AverageDDAScans
                    or SpectraFileProcessingType.AverageEverynScans)
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
            get { return mzLibSpectralAveragingOptions.NumberOfScansToAverage; }
            set { mzLibSpectralAveragingOptions.NumberOfScansToAverage = value; OnPropertyChanged(nameof(NumberOfScansToAverage)); }
        }

        public int ScanOverlap
        {
            get { return mzLibSpectralAveragingOptions.ScanOverlap; }
            set
            {
                if (value >= NumberOfScansToAverage)
                    MessageBox.Show("Overlap cannot be greater than or equal to the number of spectra averaged");
                else
                {
                    mzLibSpectralAveragingOptions.ScanOverlap = value;
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
            get { return mzLibSpectralAveragingOptions.OutputType; }
            set { mzLibSpectralAveragingOptions.OutputType = value; OnPropertyChanged(nameof(OutputType)); }
        }

        public RejectionType[] RejectionTypes { get; set; }
        public WeightingType[] WeightingTypes { get; set; }
        public SpectraFileProcessingType[] SpectraFileProcessingTypes { get; set; }
        public OutputType[] OutputTypes { get; set; }

        #endregion

        #region Commands

        #endregion

        #region Constructor

        public MzLibAveragingOptionsViewModel(MzLibSpectralAveragingOptions options)
        {
            // value initialization
            mzLibSpectralAveragingOptions = options;
            RejectionTypes = ((RejectionType[])Enum.GetValues(typeof(RejectionType))).Where(p => p != RejectionType.Thermo).ToArray();
            WeightingTypes = ((WeightingType[])Enum.GetValues(typeof(WeightingType)));
            SpectraFileProcessingTypes = ((SpectraFileProcessingType[])Enum.GetValues(typeof(SpectraFileProcessingType)));
            OutputTypes = ((OutputType[])Enum.GetValues(typeof(OutputType)));

            // command assignment
        }

        #endregion

        #region Helpers

        public void ResetDefaults()
        {
            MzLibSpectralAveragingOptions.SetDefaultValues();
            UpdateVisualRepresentation();
        }

        public void UpdateVisualRepresentation()
        {
            OnPropertyChanged(nameof(MzLibSpectralAveragingOptions));
            OnPropertyChanged(nameof(RejectionType));
            OnPropertyChanged(nameof(WeightingType));
            OnPropertyChanged(nameof(PerformNormalization));
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
    public class MzLibSpectralAveragingModel : MzLibAveragingOptionsViewModel
    {
        public static MzLibSpectralAveragingModel Instance => new MzLibSpectralAveragingModel();

        public MzLibSpectralAveragingModel() : base(new MzLibSpectralAveragingOptions())
        {

        }
    }
}
