using System.Runtime.CompilerServices;

// Type forwarding, so extracting these into GuiFunctions.Core is not a binary-breaking change.
//
// The types moved assembly but kept their namespaces, which is enough for C# callers - they recompile
// against GuiFunctions.Core through the project reference. It is not enough for XAML: every markup
// file says clr-namespace:GuiFunctions;assembly=GuiFunctions, and the WPF markup compiler resolves
// that against GuiFunctions.dll specifically. Without these forwards MetaDraw.xaml fails with
// MC3050 "Cannot find the type 'LoadingProgressViewModel'".
//
// Forwarding every public type rather than only the three XAML currently names, so a future markup
// or reflection reference to any of them keeps working.

[assembly: TypeForwardedTo(typeof(GuiFunctions.BaseViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.ClassicDeconParamsViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.CustomMdacMode))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.Databases.DownloadUniProtDatabaseFunctions))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.DeconHostModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.DeconHostViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.DeconParamsModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.DeconParamsViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.DeconvolutedSpeciesViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.FragmentViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.GptmdFilterViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.GuiGlobalParams))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.IsoDecDeconParamsViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.LoadingProgressViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.LoadingStepViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MassDifferenceAcceptorSelectionModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MassDifferenceAcceptorSelectionViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MassDifferenceAcceptorTypeViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient.ColorGradientType))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.BioPolymerCoverageResultModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.BioPolymerCoverageType))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.ChimericSpectralMatchModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.HungarianAlgorithm))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.MetaDrawTabViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MetaDraw.PlotModelStatParameters))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MultipleDeconParamsViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.MzLibExtensions))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.SelectableNotchViewModel))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.Util.CyclicalQueue<>))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.Util.ModeSwitchRequestEventArgs))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.Util.ModeSwitchResult))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.Util.SpectralMatchComparer))]
[assembly: TypeForwardedTo(typeof(GuiFunctions.XmlReaderWriter))]
