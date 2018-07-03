using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Runtime.CompilerServices;

namespace MetaMorpheusGUI
{
    public class ForTreeView : INotifyPropertyChanged
    {
        private string _status;
        private int _progress;
        private bool _isIndeterminate;

        public ForTreeView(string displayName, string id)
        {
            DisplayName = displayName;
            Id = id;
            Children = new ObservableCollection<ForTreeView>();
        }

        public event PropertyChangedEventHandler PropertyChanged;

        public ObservableCollection<ForTreeView> Children { get; private set; }

        public string Status
        {
            get
            {
                return _status;
            }
            set
            {
                _status = value;
                OnPropertyChanged();
            }
        }

        public int Progress
        {
            get
            {
                return _progress;
            }
            set
            {
                _progress = value;
                OnPropertyChanged();
            }
        }

        public string DisplayName { get; }
        public string Id { get; }

        public bool IsIndeterminate
        {
            get
            {
                return _isIndeterminate;
            }
            set
            {
                _isIndeterminate = value;
                OnPropertyChanged();
            }
        }

        protected void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}