using GuiFunctions;
using System;
using System.IO;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media.Imaging;
using System.Windows.Media;
using System.Collections.Generic;
using System.Drawing;
using Brushes = System.Windows.Media.Brushes;
using ImageFormat = System.Drawing.Imaging.ImageFormat;
using Point = System.Windows.Point;
using Size = System.Windows.Size;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChimeraAnalysisTabView.xaml
    /// </summary>
    public partial class ChimeraAnalysisTabView : UserControl
    {
        public ChimeraAnalysisTabView()
        {
            InitializeComponent();
        }

        private void ChermicDataGrid_OnSelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            if (chermicDataGrid.SelectedItem == null || sender == null || dataContext == null)
            {
                return;
            }

            ChimeraGroupViewModel chimeraGroup = dataContext.SelectedChimeraGroup;
            if (chimeraGroup == null)
            {
                MessageBox.Show("No chimera group found for this PSM");
                return;
            }
            dataContext.Ms1ChimeraPlot = new Ms1ChimeraPlot(ms1ChimeraOverlaPlot, chimeraGroup);
            dataContext.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(ms2ChimeraPlot, chimeraGroup);
            dataContext.ChimeraDrawnSequence = new ChimeraDrawnSequence(chimeraSequenceCanvas, chimeraGroup, dataContext);
        }
    }
}
