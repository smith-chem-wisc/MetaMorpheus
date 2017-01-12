using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public abstract class SearchMode
    {
        public SearchMode(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
        }

        public string FileNameAddition { get; internal set; }

        public abstract bool Accepts(double scanPrecursorMass, double peptideMass);
        internal abstract IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass);

        public override string ToString()
        {
            return FileNameAddition + " " + SearchModeString();
        }

        internal abstract string SearchModeString();
    }
}