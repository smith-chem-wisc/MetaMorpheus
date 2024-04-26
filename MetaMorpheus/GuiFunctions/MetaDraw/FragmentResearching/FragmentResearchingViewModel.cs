using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using EngineLayer;
using MassSpectrometry;
using Omics;
using Omics.Fragmentation;
using Org.BouncyCastle.Asn1.Cms;
using Proteomics.ProteolyticDigestion;

namespace GuiFunctions
{

    // TODO: Option to use prexisting fragmentation method instead of custom

    /// <summary>
    /// Class representing all GUI interactions that occur within the Additional Fragment Types section of the MetaDraw GUI
    /// </summary>
    public class FragmentResearchingViewModel : BaseViewModel
    {

        private readonly bool _isProtein;

        public FragmentResearchingViewModel(DissociationType selectedDissociationType = DissociationType.HCD, bool isProtein = true)
        {
            _isProtein = isProtein;
            DissociationTypes = new ObservableCollection<DissociationType>(Enum.GetValues<DissociationType>());
            PossibleProducts = new ObservableCollection<FragmentViewModel>(GetPossibleProducts(_isProtein));
            SelectedDissociationType = selectedDissociationType;
        }

        private ObservableCollection<FragmentViewModel> _possibleProducts;
        public ObservableCollection<FragmentViewModel> PossibleProducts
        {
            get => _possibleProducts;
            set { _possibleProducts = value; OnPropertyChanged(nameof(PossibleProducts)); }
        }

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

        private IEnumerable<FragmentViewModel> GetPossibleProducts(bool isProtein)
        {
            foreach (var product in Enum.GetValues<ProductType>())
            {
                switch (product)
                {
                    // primary ions
                    case ProductType.a:
                    case ProductType.b:
                    case ProductType.c:
                        yield return new FragmentViewModel(false, product);
                        break;
                    case ProductType.d:
                        if (!isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    case ProductType.x:
                    case ProductType.y:
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
                    case ProductType.bWaterLoss:
                    case ProductType.bBaseLoss:
                    case ProductType.cWaterLoss:
                    case ProductType.cBaseLoss:
                    case ProductType.dWaterLoss:
                    case ProductType.dBaseLoss:
                    case ProductType.wWaterLoss:
                    case ProductType.wBaseLoss:
                    case ProductType.xWaterLoss:
                    case ProductType.xBaseLoss:
                    case ProductType.yWaterLoss:
                    case ProductType.yBaseLoss:
                    case ProductType.zWaterLoss:
                    case ProductType.zBaseLoss:
                        if (!isProtein)
                            yield return new FragmentViewModel(false, product);
                        break;

                    // unsupported ions due to :
                    case ProductType.Y:
                    case ProductType.Ycore:
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
            catch (Exception e)
            {
                MessageBox.Show($"Ruh Roh Raggy - Default to HCD\n{e.Message}");

                dissociationTypeProducts = isProtein ?
                    Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD].ToArray()
                    : Omics.Fragmentation.Oligo.DissociationTypeCollection.GetRnaProductTypesFromDissociationType(DissociationType.HCD).ToArray();

                Debugger.Break();
            }

            PossibleProducts.ForEach(product => product.Use = dissociationTypeProducts.Contains(product.ProductType));
        }

        public List<MatchedFragmentIon> MatchIonsWithNewTypes(MsDataScan ms2Scan, SpectralMatch psmToRematch)
        {
            IBioPolymerWithSetMods bioPolymer = /*_isProtein ? */
                new PeptideWithSetModifications(psmToRematch.FullSequence, GlobalVariables.AllModsKnownDictionary);
            /*: new OligoWithSetMods(psmToRematch.FullSequence, GlobalVariables.AllRNAModsKnownDictionary);*/

            List<Product> allProducts = new List<Product>();
            bioPolymer.Fragment(SelectedDissociationType, FragmentationTerminus.Both, allProducts);

            // TODO: Adjust decon params for when RNA gets incorporated
            var commonParams = new CommonParameters();
            var specificMass = new Ms2ScanWithSpecificMass(ms2Scan, psmToRematch.ScanPrecursorMonoisotopicPeakMz,
                psmToRematch.ScanPrecursorCharge, psmToRematch.FullFilePath, commonParams);

            // TOCHECK: if the matchAllCharges is the correct boolean here
            return MetaMorpheusEngine.MatchFragmentIons(specificMass, allProducts, commonParams, true);
        }
    }
}
