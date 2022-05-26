using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class IonTypeForTreeViewModel : BaseViewModel
    {
        #region Public Properties

        public string GroupName { get; set; }
        public ObservableCollection<IonForTreeViewModel> Ions { get; }
        public MyEnumerator GetEnumerator()
        {
            return new MyEnumerator(this);
        }

        #endregion

        #region Constructor

        public IonTypeForTreeViewModel(string groupName, IEnumerable<ProductType> ions, bool beta)
        {
            GroupName = groupName;
            Ions = new ObservableCollection<IonForTreeViewModel>();

            foreach (var ion in ions)
            {
                Ions.Add(new IonForTreeViewModel(ion, beta));
            }
        }

        /// <summary>
        /// static instance with example data for design time editing
        /// </summary>
        public static IonTypeForTreeViewModel Instance => new IonTypeForTreeViewModel();

        public IonTypeForTreeViewModel()
        {
            GroupName = "Design Time Ions";
            Ions = new ObservableCollection<IonForTreeViewModel>();
            var ionTypes = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToList();
            foreach (var ionType in ionTypes)
            {
                Ions.Add(new IonForTreeViewModel(ionType, false));
            }
        }

        #endregion

    }
 
    public class MyEnumerator
    {
        int nIndex;
        IonTypeForTreeViewModel collection;
        public MyEnumerator(IonTypeForTreeViewModel coll)
        {
            collection = coll;
            nIndex = -1;
        }

        public bool MoveNext()
        {
            nIndex++;
            return (nIndex < collection.Ions.Count);
        }

        public IonForTreeViewModel Current => collection.Ions[nIndex];
    }
}
