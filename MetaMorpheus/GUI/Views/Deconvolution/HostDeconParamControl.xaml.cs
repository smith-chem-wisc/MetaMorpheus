using System;
using System.ComponentModel;
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
            DependencyPropertyDescriptor
                .FromProperty(ContentControl.ContentProperty, typeof(ContentControl))
                .AddValueChanged(PrecursorSpecificParams, (_, _) => PropagateFilePickerOverride());
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

        public static readonly DependencyProperty ShowGlobalPrecursorControlsProperty =
            DependencyProperty.Register(nameof(ShowGlobalPrecursorControls), typeof(bool), typeof(HostDeconParamControl), new PropertyMetadata(true));

        public bool ShowGlobalPrecursorControls
        {
            get => (bool)GetValue(ShowGlobalPrecursorControlsProperty);
            set => SetValue(ShowGlobalPrecursorControlsProperty, value);
        }

        public static readonly DependencyProperty FilePickerOverrideProperty =
            DependencyProperty.Register(nameof(FilePickerOverride), typeof(Func<string>), typeof(HostDeconParamControl),
                new PropertyMetadata(null, OnFilePickerOverrideChanged));

        public Func<string?>? FilePickerOverride
        {
            get => (Func<string?>?)GetValue(FilePickerOverrideProperty);
            set => SetValue(FilePickerOverrideProperty, value);
        }

        public FromFileDeconParamsControl? FromFilePrecursorControl =>
            PrecursorSpecificParams.Content as FromFileDeconParamsControl;

        private static void OnFilePickerOverrideChanged(DependencyObject d, DependencyPropertyChangedEventArgs e) =>
            ((HostDeconParamControl)d).PropagateFilePickerOverride();

        private void PropagateFilePickerOverride()
        {
            if (PrecursorSpecificParams?.Content is FromFileDeconParamsControl child)
                child.FilePickerOverride = FilePickerOverride;
        }
    }
}
