using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Runtime.CompilerServices;

namespace MetaMorpheusGUI
{
    public class ForTreeView : INotifyPropertyChanged
    {
        private string _Status;
        private int _Progress;
        private bool _IsIndeterminate;

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
            get { return _Status; }
            set
            {
                _Status = value;
                OnPropertyChanged();
            }
        }

        public int Progress
        {
            get { return _Progress; }
            set
            {
                _Progress = value;
                OnPropertyChanged();
            }
        }

        public string DisplayName { get; }
        public string Id { get; }

        public bool IsIndeterminate
        {
            get { return _IsIndeterminate; }
            set
            {
                _IsIndeterminate = value;
                OnPropertyChanged();
            }
        }

        protected void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}