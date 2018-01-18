using Chemistry;
using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class NeoSearchTask : MetaMorpheusTask
    {

        public bool AggregateTargetDecoyFiles { get; set; }
        public bool GenerateSplicedPeptides { get; set; }
        public bool AggregateNormalSplicedFiles { get; set; }

        public NeoSearchTask() : base(MyTask.Neo)
        {
            NeoParameters = new NeoParameters();

            CommonParameters = new CommonParameters
            {
                DoPrecursorDeconvolution = false,
                PrecursorMassTolerance = null,
                ProductMassTolerance = null
            };
            CommonParameters.DigestionParams.MinPeptideLength = 8;
            CommonParameters.DigestionParams.MaxPeptideLength = 13;
            CommonParameters.DigestionParams.Protease = GlobalEngineLevelSettings.ProteaseDictionary["non-specific"];
            CommonParameters.DigestionParams.MaxMissedCleavages = 12;
        }

        #region Public Properties

        public NeoParameters NeoParameters { get; set; }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettingsList)
        {
            myTaskResults = new MyTaskResults(this);
            //DOSTUFF
            return myTaskResults;
        }

        public NeoSearchTask Clone()
        {
            return (NeoSearchTask)this.MemberwiseClone();
        }
        #endregion Public Properties
    }
}
