using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.ComponentModel;

namespace MetaMorpheusGUI
{
    public class ExperimentalDesignForDataGrid : INotifyPropertyChanged
    {
        public ExperimentalDesignForDataGrid(string fullFilePathWithExtension)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            FileNameWithExtension = Path.GetFileName(fullFilePathWithExtension);
        }

        public string FullFilePathWithExtension { get; private set; }
        public string FileNameWithExtension { get; private set; }
        public string Condition { get { return _condition; } set { _condition = value; RaisePropertyChanged("Condition"); } }
        public string Biorep { get { return _biorep; } set { _biorep = value; RaisePropertyChanged("Biorep"); } }
        public string Fraction { get { return _fraction; } set { _fraction = value; RaisePropertyChanged("Fraction"); } }
        public string Techrep { get { return _techrep; } set { _techrep = value; RaisePropertyChanged("Techrep"); } }

        private string _condition;
        private string _biorep;
        private string _fraction;
        private string _techrep;

        public event PropertyChangedEventHandler PropertyChanged;
        private void RaisePropertyChanged(string propertyName)
        {
            if (this.PropertyChanged != null)
            {
                this.PropertyChanged(this, new PropertyChangedEventArgs(propertyName));
            }
        }
    }
}
