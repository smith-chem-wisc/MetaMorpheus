using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace GuiFunctions
{
    public abstract class DeconParamsViewModel : BaseViewModel
    {
        public DeconvolutionType DeconvolutionType => Parameters.DeconvolutionType;
        public abstract DeconvolutionParameters Parameters { get; protected set; }


        public int MinAssumedChargeState
        {
            get => Parameters.MinAssumedChargeState;
            set
            {
                Parameters.MinAssumedChargeState = value;
                OnPropertyChanged(nameof(MinAssumedChargeState));
            }
        }

        public int MaxAssumedChargeState
        {
            get => Parameters.MaxAssumedChargeState;
            set
            {
                Parameters.MaxAssumedChargeState = value;
                OnPropertyChanged(nameof(MaxAssumedChargeState));
            }
        }

        public Polarity Polarity
        {
            get => Parameters.Polarity;
            set
            {
                Parameters.Polarity = value;
                OnPropertyChanged(nameof(Polarity));
            }
        }
    }
}
