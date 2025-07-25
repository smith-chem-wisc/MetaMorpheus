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
using Transcriptomics;
using Transcriptomics.Digestion;

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
            PossibleProducts = [.. GetPossibleProducts()];

            IEnumerable<DissociationType> values;
            if (isProtein)
            {
                values = Enum.GetValues<DissociationType>()
                    .Where(p => p != DissociationType.AnyActivationType 
                    && Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType.TryGetValue(p, out var prod) 
                    && prod.Count != 0);
                SelectedDissociationType = DissociationType.HCD;
            }
            else
            {
                values = Enum.GetValues<DissociationType>()
                    .Where(p => p != DissociationType.AnyActivationType
                                && Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType.TryGetValue(p, out var prod)
                                && prod.Count != 0);
                SelectedDissociationType = DissociationType.CID;
            }
            DissociationTypes = [.. values];

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
                SetUseForFragmentsBasedUponDissociationType(value);
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

        private IEnumerable<FragmentViewModel> GetPossibleProducts()
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
                        if (!_isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    case ProductType.x:
                    case ProductType.y:
                    case ProductType.yWaterLoss:
                    case ProductType.z:
                        yield return new FragmentViewModel(false, product);
                        break;
                    case ProductType.w:
                        if (!_isProtein)
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
                        if (_isProtein)
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
                        if (!_isProtein)
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

        private void SetUseForFragmentsBasedUponDissociationType(DissociationType dissociationType)
        {
            ProductType[] dissociationTypeProducts;
            try
            {
                dissociationTypeProducts = _isProtein ?
                    Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[dissociationType].ToArray()
                    : Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType[dissociationType].ToArray();
            }
            catch (Exception)
            {
                _selectedDissociationType = DissociationType.HCD;
                dissociationTypeProducts = _isProtein ?
                    Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD].ToArray()
                    : Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD].ToArray();
            }

            PossibleProducts.ForEach(product => product.Use = dissociationTypeProducts.Contains(product.ProductType));
        }

        public List<MatchedFragmentIon> MatchIonsWithNewTypes(MsDataScan ms2Scan, SpectrumMatchFromTsv smToRematch)
        {
            if (smToRematch.FullSequence.Contains('|'))
                return smToRematch.MatchedIons;

            IBioPolymerWithSetMods bioPolymer = smToRematch.ToBioPolymerWithSetMods();

            // temp patch until fragment is no longer dependent upon parent. 
            if (!smToRematch.IsPeptide())
            {
                var nucleolyticOligo = bioPolymer as NucleolyticOligo;
                if (nucleolyticOligo != null)
                {
                    // Convert GlobalVariables.AllModsKnown (IEnumerable<Modification>) to the expected type: IDictionary<int, List<Modification>>?
                    // If you have no modifications, pass null; otherwise, build the dictionary as needed.
                    // If you need to map modifications, do so here. For now, pass null as in most usages.
                    var rna = new RNA(bioPolymer.BaseSequence);

                    // Use reflection to set the protected setter of NucleicAcid
                    var prop = typeof(NucleolyticOligo).GetProperty("NucleicAcid");
                    if (prop != null)
                    {
                        prop.SetValue(nucleolyticOligo, rna);
                    }
                }
            }

            List<Product> terminalProducts = new List<Product>();
            smToRematch.ProductsFromDissociationType()[DissociationType.Custom] = _productsToUse.ToList(); 
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

            var specificMass = new Ms2ScanWithSpecificMass(ms2Scan, smToRematch.PrecursorMz,
                smToRematch.PrecursorCharge, smToRematch.FileNameWithoutExtension, commonParams);


            var newMatches = MetaMorpheusEngine.MatchFragmentIons(specificMass, allProducts, commonParams, false);
            var existingMatches = smToRematch.MatchedIons.Where(p => _productsToUse.Contains(p.NeutralTheoreticalProduct.ProductType));
            var uniqueMatches = newMatches.Concat(existingMatches)
                .Distinct(MatchedFragmentIonComparer)
                .ToList();
            return uniqueMatches;
        }

        public static readonly IEqualityComparer<MatchedFragmentIon> MatchedFragmentIonComparer = new MatchedFragmentIonEqualityComparer();

        private class MatchedFragmentIonEqualityComparer : IEqualityComparer<MatchedFragmentIon>
        {
            public bool Equals(MatchedFragmentIon x, MatchedFragmentIon y)
            {
                if (ReferenceEquals(x, y)) return true;
                if (x is null || y is null) return false;

                return Math.Round(x.Mz, 1).Equals(Math.Round(y.Mz, 1))
                       && x.Charge == y.Charge
                       && Math.Round(x.Intensity).Equals(Math.Round(y.Intensity))
                       && x.NeutralTheoreticalProduct.FragmentNumber == y.NeutralTheoreticalProduct.FragmentNumber
                       && x.NeutralTheoreticalProduct.ProductType == y.NeutralTheoreticalProduct.ProductType
                       && x.NeutralTheoreticalProduct.SecondaryProductType == y.NeutralTheoreticalProduct.SecondaryProductType
                       && x.NeutralTheoreticalProduct.SecondaryFragmentNumber == y.NeutralTheoreticalProduct.SecondaryFragmentNumber;
            }

            public int GetHashCode(MatchedFragmentIon obj)
            {
                if (obj is null) return 0;
                int hash = 17;
                hash = hash * 23 + Math.Round(obj.Mz, 1).GetHashCode();
                hash = hash * 23 + obj.Charge.GetHashCode();
                hash = hash * 23 + Math.Round(obj.Intensity).GetHashCode();
                hash = hash * 23 + obj.NeutralTheoreticalProduct.FragmentNumber.GetHashCode();
                hash = hash * 23 + obj.NeutralTheoreticalProduct.ProductType.GetHashCode();
                hash = hash * 23 + (obj.NeutralTheoreticalProduct.SecondaryProductType?.GetHashCode() ?? 0);
                if (obj.NeutralTheoreticalProduct.SecondaryFragmentNumber != null)
                    hash = hash * 23 + obj.NeutralTheoreticalProduct.SecondaryFragmentNumber.GetHashCode();

                return hash;
            }
        }
    }
}
