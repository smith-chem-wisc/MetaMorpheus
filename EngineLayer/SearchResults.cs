using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class SearchResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public SearchResults(PsmParent[][] psms, MetaMorpheusEngine searchParams) : base(searchParams)
        {
            Psms = new List<PsmParent>[psms.Length];
            for (int j = 0; j < psms.Length; j++)
                Psms[j] = psms[j].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();
        }

        #endregion Public Constructors

        #region Public Properties

        public List<PsmParent>[] Psms { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(base.ToString());
            return sb.ToString();
        }

        #endregion Public Methods

    }
}