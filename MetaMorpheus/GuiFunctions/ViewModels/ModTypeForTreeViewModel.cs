using GuiFunctions;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Windows.Media;

namespace GuiFunctions
{
    public class ModTypeForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        private bool? _isChecked;

        #endregion

        #region Public Properties
        
        public bool Expanded { get; set; }
        public string DisplayName { get; }
        public ObservableCollection<ModForTreeViewModel> Children { get; }
        public Brush Background { get; }
        public MyEnumerator GetEnumerator()
        {
            return new MyEnumerator(this);
        }
        public bool? Use
        {
            get
            {
                return _isChecked;
            }
            set
            {
                _isChecked = value;
                if (value.HasValue)
                    foreach (var child in Children)
                        child.Use = (value.Value);
                OnPropertyChanged(nameof(Use));
            }
        }

        #endregion

        #region Constructor

        public ModTypeForTreeViewModel(string displayName, bool bad)
        {
            Children = new ObservableCollection<ModForTreeViewModel>();
            Expanded = false;
            DisplayName = displayName;
            if (bad)
                Background = new SolidColorBrush(Colors.Red);
            else
                Background = new SolidColorBrush(Colors.Transparent);
        }

        #endregion

        public void VerifyCheckState()
        {
            bool? state = null;
            for (int i = 0; i < Children.Count; ++i)
            {
                bool current = Children[i].Use;
                if (i == 0)
                    state = current;
                else if (state != current)
                {
                    state = null;
                    break;
                }
            }
            Use = state;
        }

        public class MyEnumerator
        {
            int nIndex;
            ModTypeForTreeViewModel collection;
            public MyEnumerator(ModTypeForTreeViewModel coll)
            {
                collection = coll;
                nIndex = -1;
            }

            public bool MoveNext()
            {
                nIndex++;
                return (nIndex < collection.Children.Count);
            }

            public ModForTreeViewModel Current => collection.Children[nIndex];
        }

    }
}