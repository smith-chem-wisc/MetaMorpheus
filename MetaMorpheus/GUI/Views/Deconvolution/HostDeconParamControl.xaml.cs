using System.Windows;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for HostDeconParamControl.xaml
    /// </summary>
    public partial class HostDeconParamControl : UserControl
    {
        public HostDeconParamControl()
        {
            InitializeComponent();
        }

        public static readonly DependencyProperty GroupBoxHeaderProperty =
            DependencyProperty.Register(nameof(GroupBoxHeader), typeof(string), typeof(HostDeconParamControl), new PropertyMetadata("MS1 Deconvolution"));

        public string GroupBoxHeader
        {
            get => (string)GetValue(GroupBoxHeaderProperty);
            set => SetValue(GroupBoxHeaderProperty, value);
        }

        public static readonly DependencyProperty ShowMs2SectionProperty =
            DependencyProperty.Register(nameof(ShowMs2Section), typeof(bool), typeof(HostDeconParamControl), new PropertyMetadata(true));

        public bool ShowMs2Section
        {
            get => (bool)GetValue(ShowMs2SectionProperty);
            set => SetValue(ShowMs2SectionProperty, value);
        }
    }
}
