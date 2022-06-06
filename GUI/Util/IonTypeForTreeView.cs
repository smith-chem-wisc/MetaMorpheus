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
    public class IonTypeForTreeView : BaseView
    {
        #region Public Properties

        public string GroupName { get; set; }
        public ObservableCollection<IonForTreeView> Ions { get; }
        public MyEnumerator GetEnumerator()
        {
            return new MyEnumerator(this);
        }

        #endregion

        #region Constructor

        public IonTypeForTreeView(string groupName, IEnumerable<ProductType> ions, bool beta)
        {
            GroupName = groupName;
            Ions = new ObservableCollection<IonForTreeView>();

            foreach (var ion in ions)
            {
                Ions.Add(new IonForTreeView(ion, beta));
            }
        }

        #endregion

    }
 
    public class MyEnumerator
    {
        int nIndex;
        IonTypeForTreeView collection;
        public MyEnumerator(IonTypeForTreeView coll)
        {
            collection = coll;
            nIndex = -1;
        }

        public bool MoveNext()
        {
            nIndex++;
            return (nIndex < collection.Ions.Count);
        }

        public IonForTreeView Current => collection.Ions[nIndex];
    }
}
