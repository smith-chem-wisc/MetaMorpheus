using OxyPlot;
using Omics.Fragmentation;

namespace GuiFunctions
{
    public class IonForTreeViewModel : ColorForTreeViewModel
    {
        public ProductType IonType { get; set; }    
        public bool IsBeta { get; set; }

        public IonForTreeViewModel(ProductType type, bool beta) : base(type + " - Ion", beta ? MetaDrawSettings.BetaProductTypeToColor[type] : MetaDrawSettings.ProductTypeToColor[type])
        {
            IonType = type;
            IsBeta = beta;
        }

        // should only be used for unannotated peak or internal ion color as it is not an actual product type
        public IonForTreeViewModel(string type, bool beta)
            : base(
                type,
                type.Equals("Unannotated Peak") ? MetaDrawSettings.UnannotatedPeakColor :
                type.Equals("Internal Ion") ? MetaDrawSettings.InternalIonColor :
                OxyColors.Transparent // fallback color if type is not recognized
            )
        {
            IsBeta = beta;
            IonType = default; // Not an actual ProductType
        }
    }
}
