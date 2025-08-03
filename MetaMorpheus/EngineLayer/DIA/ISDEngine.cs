using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISDEngine : DIAEngine
    {
        public ISDEngine(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(DIAparameters, dataFile, commonParameters, fileSpecificParameters, nestedIds)
        {
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            return new MetaMorpheusEngineResults(this);
        }

        public static void PreprocessIsdScans(MsDataScan[] isdScans, List<double> isdEnergies)
        {
            int numberOfScansPerCycle = isdEnergies.Count;
            foreach(var scan in isdScans)
            {
                scan.SetMsnOrder(2);
            }
        }
    }
}
