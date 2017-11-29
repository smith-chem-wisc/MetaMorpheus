using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Runtime.CompilerServices;

namespace MetaMorpheusGUI
{
    public class ForTreeView : INotifyPropertyChanged
    {
        #region Private Fields

        private string status;
        private int progress;
        private bool isIndeterminate;

        #endregion Private Fields

        #region Public Constructors

        public ForTreeView(string displayName, string id)
        {
            DisplayName = displayName;
            Id = id;
            Children = new ObservableCollection<ForTreeView>();
        }

        #endregion Public Constructors

        #region Public Events

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion Public Events

        #region Public Properties

        public ObservableCollection<ForTreeView> Children { get; private set; }

        public string Status
        {
            get { return status; }
            set
            {
                status = value;
                OnPropertyChanged();
            }
        }

        public int Progress
        {
            get { return progress; }
            set
            {
                progress = value;
                OnPropertyChanged();
            }
        }

        public string DisplayName { get; }
        public string Id { get; }

        public bool IsIndeterminate
        {
            get { return isIndeterminate; }
            set
            {
                isIndeterminate = value;
                OnPropertyChanged();
            }
        }

        #endregion Public Properties

        #region Protected Methods

        protected void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

        #endregion Protected Methods
    }
}