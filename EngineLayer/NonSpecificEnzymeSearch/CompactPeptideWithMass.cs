using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using Proteomics;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    class CompactPeptideWithMass
    {
        #region Private Fields

        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        #endregion Private Fields

        #region Public Constructors

        public CompactPeptideWithMass(CompactPeptide compactPeptide, double precursorMass)
        {
            this.compactPeptide = compactPeptide;
            this.precursorMass = precursorMass;
        }

        public CompactPeptideWithMass(double[] CTerminalMasses, double[] NTerminalMasses, double MonoisotopicMassIncludingFixedMods)
        {
            this.CTerminalMasses = CTerminalMasses;
            this.NTerminalMasses = NTerminalMasses;
            this.MonoisotopicMassIncludingFixedMods = MonoisotopicMassIncludingFixedMods;
        }

        #endregion Public Constructors

        #region Public Properties

        public double[] CTerminalMasses { get; }
        public double[] NTerminalMasses { get; }
        public double MonoisotopicMassIncludingFixedMods { get; set; }
        public readonly CompactPeptide compactPeptide;
        public readonly double precursorMass;

        #endregion Public Properties
    }
}
