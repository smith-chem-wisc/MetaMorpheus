﻿using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class SinglePpmAroundZeroSearchMode : SearchMode
    {

        #region Private Fields

        private readonly double ppmTolerance;

        #endregion Private Fields

        #region Public Constructors

        public SinglePpmAroundZeroSearchMode(double ppmTolerance) : base(ppmTolerance + "ppmAroundZero")
        {
            this.ppmTolerance = ppmTolerance;
        }

        #endregion Public Constructors

        #region Public Methods

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return Math.Abs((scanPrecursorMass - peptideMass) / (peptideMass) * 1e6) < ppmTolerance ? 0 : -1;
        }

        public override IEnumerable<Tuple<DoubleRange, int>> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            var diff = ppmTolerance / 1e6 * peptideMonoisotopicMass;
            yield return new Tuple<DoubleRange, int>(new DoubleRange(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff), 0);
        }

        public override string ToString()
        {
            return FileNameAddition + " ppmAroundZero " + ppmTolerance;
        }

        #endregion Public Methods

    }
}