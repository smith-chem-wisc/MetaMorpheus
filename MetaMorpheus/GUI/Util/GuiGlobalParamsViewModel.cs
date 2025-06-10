using GuiFunctions;
using System.Windows;

namespace MetaMorpheusGUI.Util
{
    /// <summary>
    /// Class to represent Global parameters that can change at runtime. 
    /// Only one instance of this class should exist, and that is found in <see cref="UpdateGUISettings.Globals"/>
    /// </summary>
    public class GuiGlobalParamsViewModel(GuiGlobalParams parameters) : BaseViewModel
    {
        public bool IsRnaMode
        {
            get => parameters.IsRnaMode;
            set { parameters.IsRnaMode = value; OnPropertyChanged(nameof(IsRnaMode)); }
        }

        public bool ExposeRnaToUser
        {
            get => parameters.ExposeRnaToUser;
            set { parameters.ExposeRnaToUser = value; OnPropertyChanged(nameof(ExposeRnaToUser)); }
        }
    }
}
