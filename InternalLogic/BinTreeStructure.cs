﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class BinTreeStructure
    {

        #region Private Fields

		private const int minAdditionalInBin = 1;

        #endregion Private Fields

        #region Public Properties

        public List<Bin> FinalBins { get; private set; }

        #endregion Public Properties

        #region Internal Methods

        internal void GenerateBins(List<NewPsmWithFdr> targetAndDecoyMatches, double dc)
        {
            List<double> listOfMassShifts = targetAndDecoyMatches.Select(b => b.thisPSM.ScanPrecursorMass - b.thisPSM.PeptideMonoisotopicMass).OrderBy(b => b).ToList();
            double minMassShift = listOfMassShifts.Min();
            double maxMassShift = listOfMassShifts.Max();

			Console.WriteLine("Need those to be the same: "+string.Join(",", listOfMassShifts));

            int[] p = new int[listOfMassShifts.Count];

            int firstIndex = 0;
            int lastIndex = 0;
            for (int i = 0; i < listOfMassShifts.Count; i++)
            {
                var thisMassShift = listOfMassShifts[i];

                while (thisMassShift - listOfMassShifts[firstIndex] > dc)
                    firstIndex++;
                while (lastIndex + 1 < listOfMassShifts.Count && listOfMassShifts[lastIndex + 1] - thisMassShift <= dc)
                    lastIndex++;

                p[i] = lastIndex - firstIndex;
            }

            int maxP = p.Max();
            double[] sigma = new double[listOfMassShifts.Count];

            for (int i = 0; i < listOfMassShifts.Count; i++)
            {
                var thisMassShift = listOfMassShifts[i];
                var thisP = p[i];
                if (thisP == maxP)
                    sigma[i] = Math.Max(maxMassShift - thisMassShift, thisMassShift - minMassShift);
                else
                {
                    // SIGMA IS THE DISTANCE TO THE CLOSEST MASS SHIFT THAT HAS A HIGHER P VALUE THAN ITSELF

                    sigma[i] = getSigma(thisMassShift, thisP, i, listOfMassShifts, p);
                }
            }

            var listokbin = new List<OkBin>();
            for (int i = 0; i < sigma.Count(); i++)
                listokbin.Add(new OkBin(listOfMassShifts[i], sigma[i], p[i]));


			Console.WriteLine("listokbin: " + string.Join(",", listokbin.Select(b=>b.massShift)));
			Console.WriteLine("listokbin sigma: " + string.Join(",", listokbin.Select(b => b.sigma)));
			Console.WriteLine("listokbin p: " + string.Join(",", listokbin.Select(b => b.p)));

            var prelimBins = new HashSet<double>();
            foreach (OkBin okbin in listokbin.OrderByDescending(b => b.p))
            {
                if (okbin.sigma < dc || okbin.p < minAdditionalInBin)
                    continue;
                bool add = true;
                foreach (double a in prelimBins)
                {
                    if (Math.Abs(okbin.massShift - a) <= dc)
                    {
                        add = false;
                        break;
                    }
                }
                if (add)
                    prelimBins.Add(okbin.massShift);
            }


			Console.WriteLine("prelimBins: " + string.Join(",", prelimBins));


            var forFinalBins = new Dictionary<double, List<double>>();
            foreach (double ok in prelimBins)
                forFinalBins.Add(ok, new List<double>());
            foreach (double a in listOfMassShifts)
                foreach (double b in prelimBins)
                    if (Math.Abs(a - b) <= dc)
                        forFinalBins[b].Add(a);

            FinalBins = forFinalBins.Select(b => new Bin(b.Value.Average())).ToList();

            for (int i = 0; i < targetAndDecoyMatches.Count; i++)
            {
                foreach (Bin bin in FinalBins)
                    if (Math.Abs(targetAndDecoyMatches[i].thisPSM.ScanPrecursorMass - targetAndDecoyMatches[i].thisPSM.PeptideMonoisotopicMass - bin.MassShift) <= dc)
                        bin.Add(targetAndDecoyMatches[i]);
            }

            FinalBins = FinalBins.Where(b => b.Count > 1).ToList();


			Console.WriteLine("Final bins: " + string.Join(",", FinalBins.Select(b=>b.MassShift)));


        }

        #endregion Internal Methods

        #region Private Methods

        private static double getSigma(double thisMassShift, int thisP, int i, List<double> listOfMassShifts, int[] p)
        {
            int currentDown = i - 1;
            int currentUp = i + 1;
            int lookingAtP;
            double distDown, distUp;
            while (true)
            {
                distDown = currentDown == -1 ? double.MaxValue : thisMassShift - listOfMassShifts[currentDown];
                distUp = currentUp == listOfMassShifts.Count ? double.MaxValue : listOfMassShifts[currentUp] - thisMassShift;
                if (distDown < distUp)
                {
                    lookingAtP = p[currentDown];
                    if (lookingAtP > thisP)
                    {
                        return distDown;
                    }
                    currentDown--;
                }
                else
                {
                    lookingAtP = p[currentUp];
                    if (lookingAtP > thisP)
                    {
                        return distUp;
                    }
                    currentUp++;
                }
            }
        }

        #endregion Private Methods

    }
}