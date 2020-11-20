using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public static class SpectralLibrarySearchFunction
    {
      
        //1024TestDoneCompareFunction
        public static double matchedSpectraCompare(List<MatchedFragmentIon> standardSpectra, List<MatchedFragmentIon> spectraToCompare)
        {

            double[] mz1 = standardSpectra.Select(b => b.Mz).ToArray();
            double intensitySum1 = standardSpectra.Select(b => b.Intensity).Sum();
            double[] intensity1 = standardSpectra.Select(b => b.Intensity / intensitySum1).ToArray();
            //Console.WriteLine(mz1.Length + "  " + intensity1.Length);
            Array.Sort(mz1, intensity1);

            double[] mz2 = spectraToCompare.Select(b => b.Mz).ToArray();
            double intensitySum2 = spectraToCompare.Select(b => b.Intensity).Sum();
            double[] intensity2 = spectraToCompare.Select(b => b.Intensity / intensitySum2).ToArray();
            Array.Sort(mz2, intensity2);
            //Console.WriteLine(mz2.Length + "  " + intensity2.Length);

            var commonNumbers = mz1.Union(mz2).ToArray();
            double min = commonNumbers.Min();
            double max = commonNumbers.Max();
            int roundMin = (int)min;
            int roundMax = (int)max + 1;
            //Console.WriteLine(roundMin + "  " + roundMax);

            //convert spectra to vectors
            List<double> vector1 = new List<double>();
            List<double> vector2 = new List<double>();

            int i = 0; //iterate through mz1
            int k = 0; //iterate through bin
            double oneMz = mz1[0];
            double oneIntensity = intensity1[0];
            //find where peaks match
            while (roundMin + k * 0.5 < roundMax)
            {
                List<double> x1 = new List<double>();
                while (i < mz1.Length && roundMin + k * 0.5 <= oneMz && oneMz < roundMin + k * 0.5 + 0.5)
                {
                    x1.Add(oneIntensity);
                    i++;
                    if (i != mz1.Length)
                    {
                        oneMz = mz1[i];
                        oneIntensity = intensity1[i];
                    }
                }
                vector1.Add(x1.Sum());
                k++;
            }

            int j = 0; //iterate through mz2
            int n = 0; //iterate through bin
            double twoMz = mz2[0];
            double twoIntensity = intensity2[0];
            while (roundMin + n * 0.5 < roundMax)
            {
                List<double> x2 = new List<double>();
                while (j < mz2.Length && roundMin + n * 0.5 <= twoMz && twoMz < roundMin + n * 0.5 + 0.5)
                {
                    x2.Add(twoIntensity);
                    j++;
                    if (j != mz2.Length)
                    {
                        twoMz = mz2[j];
                        twoIntensity = intensity2[j];
                    }
                }
                vector2.Add(x2.Sum());
                n++;
            }

            //numerator of dot product
            double numerator = 0;
            for (i = 0; i < vector1.Count; i++)
            {
                numerator += vector1[i] * vector2[i];
            }

            //denominator of dot product
            double denominator = Math.Sqrt(vector1.Sum(x => x * x)) * Math.Sqrt(vector2.Sum(x => x * x));

            var score = numerator / denominator;
            return score;
        }

    }
}
