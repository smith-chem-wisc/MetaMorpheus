using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using MassSpectrometry;

namespace TaskLayer.MbrAnalysis
{
    internal class RecoveredMs2ScanWithSpecificMass : Ms2ScanWithSpecificMass
    {
        public readonly bool DeconvolutedOrTheoreticalPrecursor;
        public RecoveredMs2ScanWithSpecificMass(MsDataScan dataScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath,
            CommonParameters commonParam, bool deconvolutedOrTheoreticalPrecursor, IsotopicEnvelope[] neutralExperimentalFragments = null) :
            base(dataScan, precursorMonoisotopicPeakMz, precursorCharge, fullFilePath, commonParam,
                neutralExperimentalFragments)
        {
            DeconvolutedOrTheoreticalPrecursor = deconvolutedOrTheoreticalPrecursor;
        }
    }
}
