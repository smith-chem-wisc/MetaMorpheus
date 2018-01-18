using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class NeoParameters
    {
        #region Public Constructors

        public NeoParameters()
        {
            Calibrate = true;
            GPTMD = true;
            TargetSearch = true;
            DecoySearch = true;
            SearchNTerminus = true;
            SearchCTerminus = true;
            MaxMissedConsecutiveCleavages = 2;
            MaxMissedCleavages = 5;
            MaxCandidatesPerSpectrum = 2000;
            MinDistanceAllowed = 1;
            MaxDistanceAllowed = 25;
            NormalCis = true;
            ReverseCis = true;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Calibrate { get; set; }
        public double? PrecursorTolerancePPM { get; set; }
        public double? ProductTolerancePPM { get; set; }
        public bool GPTMD { get; set; }
        public bool TargetSearch { get; set; }
        public string TargetFilePath { get; set; }
        public bool DecoySearch { get; set; }
        public string DecoyFilePath { get; set; }

        public bool SearchNTerminus { get; set; }
        public string NFilePath { get; set; }
        public bool SearchCTerminus { get; set; }
        public string CFilePath { get; set; }

        public int MaxMissedConsecutiveCleavages { get; set; }
        public int MaxMissedCleavages { get; set; }
        public int MaxCandidatesPerSpectrum { get; set; }

        public int MinDistanceAllowed { get; set; }
        public int MaxDistanceAllowed { get; set; }
        public bool NormalCis { get; set; }
        public bool ReverseCis { get; set; }

        #endregion Public Properties
    }
}
