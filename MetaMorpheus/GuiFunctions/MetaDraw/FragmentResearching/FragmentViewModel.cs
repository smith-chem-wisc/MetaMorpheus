using Omics.Fragmentation;

namespace GuiFunctions
{
    public class FragmentViewModel : BaseViewModel
    {
        public ProductType ProductType { get; }
        public string TypeString => ProductType.ToString();
        private bool use;

        public bool Use
        {
            get => use;
            set { use = value; OnPropertyChanged(nameof(Use)); }
        }

        public FragmentViewModel(bool use, ProductType type)
        {
            Use = use;
            ProductType = type;
        }
    }

    public class FragmentModel : FragmentViewModel
    {
        public static FragmentModel Instance => new FragmentModel();

        public FragmentModel() : base(true, ProductType.a)
        {

        }
    }
}
