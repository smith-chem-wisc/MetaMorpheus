using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Linq;
using System.Windows;
using Easy.Common.Extensions;
using EngineLayer;
using iText.StyledXmlParser.Jsoup;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;

namespace GuiFunctions
{
    /// <summary>
    /// Class representing all GUI interactions that occur within the Additional Fragment Types section of the MetaDraw GUI
    /// </summary>
    public class FragmentationReanalysisViewModel : BaseViewModel
    {
        private readonly bool _isProtein;

        public FragmentationReanalysisViewModel(bool isProtein = true)
        {
            _isProtein = isProtein;
            UseInternalIons = false;
            MinInternalIonLength = 10;
            ProductIonMassTolerance = 20;
            if (isProtein)
            {
                DissociationTypes = new ObservableCollection<DissociationType>(Enum.GetValues<DissociationType>()
                    .Where(p => p != DissociationType.AnyActivationType && Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType.TryGetValue(p, out var prod) && prod.Count != 0));
                PossibleProducts = new ObservableCollection<FragmentViewModel>(GetPossibleProducts(_isProtein));
                SelectedDissociationType = DissociationType.HCD;
            }
            else
            {
                // TODO: Possible product types for RNA
                PossibleProducts = new ObservableCollection<FragmentViewModel>(GetPossibleProducts(_isProtein));
                SelectedDissociationType = DissociationType.LowCID;
            }

        }

        private ObservableCollection<FragmentViewModel> _possibleProducts;
        public ObservableCollection<FragmentViewModel> PossibleProducts
        {
            get => _possibleProducts;
            set { _possibleProducts = value; OnPropertyChanged(nameof(PossibleProducts)); }
        }

        private IEnumerable<ProductType> _productsToUse => PossibleProducts.Where(p => p.Use).Select(p => p.ProductType);

        private bool _persist;
        public bool Persist
        {
            get => _persist;
            set { _persist = value; OnPropertyChanged(nameof(Persist)); }
        }

        private ObservableCollection<DissociationType> _dissociationTypes;
        public ObservableCollection<DissociationType> DissociationTypes
        {
            get => _dissociationTypes;
            set { _dissociationTypes = value; OnPropertyChanged(nameof(DissociationTypes)); }
        }

        private DissociationType _selectedDissociationType;
        public DissociationType SelectedDissociationType
        {
            get => _selectedDissociationType;
            set
            {
                _selectedDissociationType = value;
                SetUseForFragmentsBasedUponDissociationType(value, _isProtein);
                OnPropertyChanged(nameof(SelectedDissociationType));
            }
        }

        private int _minInternalIonLength;
        public int MinInternalIonLength
        {
            get => _minInternalIonLength;
            set { _minInternalIonLength = value; OnPropertyChanged(nameof(MinInternalIonLength)); }
        }

        private bool _useInternalIons;
        public bool UseInternalIons
        {
            get => _useInternalIons; 
            set { _useInternalIons = value; OnPropertyChanged(nameof(UseInternalIons)); }
        }

        private double _productIonMassTolerance;

        public double ProductIonMassTolerance
        {
            get => _productIonMassTolerance;
            set { _productIonMassTolerance = value; OnPropertyChanged(nameof(ProductIonMassTolerance)); }
        }

        private IEnumerable<FragmentViewModel> GetPossibleProducts(bool isProtein)
        {
            foreach (var product in Enum.GetValues<ProductType>())
            {
                switch (product)
                {
                    // primary ions
                    case ProductType.a:
                    case ProductType.b:
                    case ProductType.bWaterLoss:
                    case ProductType.c:
                        yield return new FragmentViewModel(false, product);
                        break;
                    case ProductType.d:
                        if (!isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    case ProductType.x:
                    case ProductType.y:
                    case ProductType.yWaterLoss:
                    case ProductType.z:
                        yield return new FragmentViewModel(false, product);
                        break;
                    case ProductType.w:
                        if (!isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    case ProductType.M:
                        yield return new FragmentViewModel(false, product);
                        break;

                    // protein specific ions
                    case ProductType.aStar:
                    case ProductType.aDegree:
                    case ProductType.bAmmoniaLoss:
                    case ProductType.yAmmoniaLoss:
                    case ProductType.zPlusOne:
                    case ProductType.zDot:
                        if (isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    // rna specific ions
                    case ProductType.aWaterLoss:
                    case ProductType.aBaseLoss:
                    case ProductType.bBaseLoss:
                    case ProductType.cWaterLoss:
                    case ProductType.cBaseLoss:
                    case ProductType.dWaterLoss:
                    case ProductType.dBaseLoss:
                    case ProductType.wWaterLoss:
                    case ProductType.wBaseLoss:
                    case ProductType.xWaterLoss:
                    case ProductType.xBaseLoss:
                    case ProductType.yBaseLoss:
                    case ProductType.zWaterLoss:
                    case ProductType.zBaseLoss:
                        if (!isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    // unsupported ions due to :
                    case ProductType.Y:
                    case ProductType.Ycore:
                    case ProductType.D:
                        break;

                    // default case
                    default:
                        yield return new FragmentViewModel(false, product);
                        break;
                }
            }
        }

        private void SetUseForFragmentsBasedUponDissociationType(DissociationType dissociationType, bool isProtein)
        {
            ProductType[] dissociationTypeProducts;
            try
            {
                dissociationTypeProducts = isProtein ?
                    Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[dissociationType].ToArray()
                    : Omics.Fragmentation.Oligo.DissociationTypeCollection.GetRnaProductTypesFromDissociationType(dissociationType).ToArray();
            }
            catch (Exception)
            {
                _selectedDissociationType = DissociationType.HCD;
                dissociationTypeProducts = isProtein ?
                    Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD].ToArray()
                    : Omics.Fragmentation.Oligo.DissociationTypeCollection.GetRnaProductTypesFromDissociationType(DissociationType.HCD).ToArray();
            }

            PossibleProducts.ForEach(product => product.Use = dissociationTypeProducts.Contains(product.ProductType));
        }

        public List<MatchedFragmentIon> MatchIonsWithNewTypes(MsDataScan ms2Scan, SpectrumMatchFromTsv psmToRematch)
        {
            if (psmToRematch.FullSequence.Contains('|'))
                return psmToRematch.MatchedIons;

            IBioPolymerWithSetMods bioPolymer = psmToRematch.ToBioPolymerWithSetMods();

            List<Product> terminalProducts = new List<Product>();
            Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = _productsToUse.ToList(); 
            bioPolymer.Fragment(DissociationType.Custom, FragmentationTerminus.Both, terminalProducts);

            List<Product> internalProducts = new List<Product>();
            if (UseInternalIons && bioPolymer is PeptideWithSetModifications) // internal ions are not currently implemented for RNA
            {
                Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = _productsToUse.ToList();
                bioPolymer.FragmentInternally(DissociationType.Custom, MinInternalIonLength, internalProducts);
            }
            var allProducts = terminalProducts.Concat(internalProducts).ToList();

            // TODO: Adjust decon params for when RNA gets incorporated
            var commonParams = new CommonParameters();
            if (Math.Abs(commonParams.ProductMassTolerance.Value - ProductIonMassTolerance) > 0.00001)
                commonParams.ProductMassTolerance = new PpmTolerance(ProductIonMassTolerance);

            var specificMass = new Ms2ScanWithSpecificMass(ms2Scan, psmToRematch.PrecursorMz,
                psmToRematch.PrecursorCharge, psmToRematch.FileNameWithoutExtension, commonParams);

            return MetaMorpheusEngine.MatchFragmentIons(specificMass, allProducts, commonParams, false)
                .Union(psmToRematch.MatchedIons.Where(p => _productsToUse.Contains(p.NeutralTheoreticalProduct.ProductType)))
                .ToList();
        }
    }
}
