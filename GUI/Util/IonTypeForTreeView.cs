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
    public class IonTypeForTreeView : INotifyPropertyChanged
    {
        #region Public Properties

        public string GroupName { get; set; }
        public ObservableCollection<IonForTreeView> Ions { get; }

        #endregion

        #region Constructor

        public IonTypeForTreeView(string groupName)
        {
            var products = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToList();
            var common = products.Where(p => p.ToString().Equals("a") || p.ToString().Equals("b") || p.ToString().Equals("c")
                                          || p.ToString().Equals("x") || p.ToString().Equals("y") || p.ToString().Equals("zDot"));
            var lessCommon = products.Where(p => !p.ToString().Equals("a") || !p.ToString().Equals("b") || !p.ToString().Equals("c")
                                          || !p.ToString().Equals("x") || !p.ToString().Equals("y") || !p.ToString().Equals("zDot"));

            GroupName = groupName;
            Ions = new ObservableCollection<IonForTreeView>();
            foreach (var name in common)
            {
                IonForTreeView ion = new IonForTreeView(name);
                Ions.Add(ion);
            }
        }

        #endregion

        #region Commands

        #endregion

        #region PropertyChanged

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion

    }
}
