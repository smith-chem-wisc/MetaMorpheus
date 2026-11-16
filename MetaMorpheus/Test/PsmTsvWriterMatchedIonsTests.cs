using Chemistry;
using EngineLayer;
using NUnit.Framework;
using Omics.Fragmentation;
using Readers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Test
{
    /// <summary>
    /// AddMatchedIonsData formats every matched ion into five columns of every PSM row. It now formats each value
    /// straight into presized builders instead of through intermediate strings, and must write exactly the same text,
    /// in whatever culture the process runs.
    /// </summary>
    [TestFixture]
    internal static class PsmTsvWriterMatchedIonsTests
    {
        [Test]
        [TestCase("")]
        [TestCase("en-US")]
        [TestCase("de-DE")]
        [TestCase("fr-FR")]
        public static void MatchedIonColumnsAreUnchanged(string cultureName)
        {
            CultureInfo original = CultureInfo.CurrentCulture;
            try
            {
                CultureInfo.CurrentCulture = new CultureInfo(cultureName);
                var random = new Random(31);
                ProductType[] types = { ProductType.b, ProductType.y, ProductType.c, ProductType.zDot, ProductType.M, ProductType.D };

                for (int trial = 0; trial < 3000; trial++)
                {
                    List<MatchedFragmentIon> ions = null;
                    if (trial % 50 != 0) // every 50th PSM has no ions at all
                    {
                        ions = new List<MatchedFragmentIon>();
                        int count = trial % 49 == 0 ? 0 : random.Next(1, 400);
                        for (int i = 0; i < count; i++)
                        {
                            ProductType type = types[random.Next(types.Length)];
                            ProductType? secondary = random.Next(8) == 0 ? types[random.Next(types.Length)] : null;
                            double neutralMass = random.Next(200) == 0 ? 0 : random.NextDouble() * 20000;
                            var product = new Product(type, FragmentationTerminus.Both, neutralMass, random.Next(1, 60), random.Next(1, 60),
                                random.Next(10) == 0 ? 97.976896 : 0, secondary, secondary == null ? 0 : random.Next(1, 60));
                            int charge = random.Next(1, 30);
                            double mz = neutralMass.ToMz(charge) + (random.NextDouble() - 0.5) * 0.05;
                            double intensity = random.Next(100) == 0 ? double.NaN : Math.Pow(10, random.NextDouble() * 10);
                            ions.Add(new MatchedFragmentIon(product, mz, intensity, charge));
                        }
                    }

                    var expected = new Dictionary<string, string>();
                    PreviousAddMatchedIonsData(expected, ions);
                    var actual = new Dictionary<string, string>();
                    PsmTsvWriter.AddMatchedIonsData(actual, ions);

                    Assert.That(actual.Keys, Is.EqualTo(expected.Keys), "same columns in the same order");
                    foreach (var column in expected)
                    {
                        Assert.That(actual[column.Key], Is.EqualTo(column.Value), $"trial {trial}, column {column.Key}");
                    }
                }
            }
            finally
            {
                CultureInfo.CurrentCulture = original;
            }
        }

        /// <summary>
        /// Verbatim copy of PsmTsvWriter.AddMatchedIonsData before it formatted straight into presized builders
        /// </summary>
        private static void PreviousAddMatchedIonsData(Dictionary<string, string> s, List<MatchedFragmentIon> matchedIons)
        {
            bool nullPsm = matchedIons == null;

            StringBuilder seriesStringBuilder = new StringBuilder();
            StringBuilder mzStringBuilder = new StringBuilder();
            StringBuilder fragmentDaErrorStringBuilder = new StringBuilder();
            StringBuilder fragmentPpmErrorStringBuilder = new StringBuilder();
            StringBuilder fragmentIntensityStringBuilder = new StringBuilder();
            List<StringBuilder> stringBuilders = new List<StringBuilder> { seriesStringBuilder, mzStringBuilder, fragmentDaErrorStringBuilder, fragmentPpmErrorStringBuilder, fragmentIntensityStringBuilder };

            if (!nullPsm)
            {
                // using ", " instead of "," improves human readability
                const string delimiter = ", ";

                var matchedIonsGroupedByProductType = matchedIons.GroupBy(x => new { x.NeutralTheoreticalProduct.ProductType, x.NeutralTheoreticalProduct.SecondaryProductType }).ToList();

                foreach (var productType in matchedIonsGroupedByProductType)
                {
                    var products = productType.OrderBy(p => p.NeutralTheoreticalProduct.FragmentNumber)
                        .ToList();

                    stringBuilders.ForEach(p => p.Append("["));

                    for (int i = 0; i < products.Count; i++)
                    {
                        MatchedFragmentIon ion = products[i];
                        string ionLabel;

                        double massError = ion.Mz.ToMass(ion.Charge) - ion.NeutralTheoreticalProduct.NeutralMass;
                        double ppmMassError = massError / ion.NeutralTheoreticalProduct.NeutralMass * 1e6;

                        ionLabel = ion.Annotation;

                        // append ion label
                        seriesStringBuilder.Append(ionLabel);

                        // append experimental m/z
                        mzStringBuilder.Append(ionLabel + ":" + ion.Mz.ToString("F5"));

                        // append absolute mass error
                        fragmentDaErrorStringBuilder.Append(ionLabel + ":" + massError.ToString("F5"));

                        // append ppm mass error
                        fragmentPpmErrorStringBuilder.Append(ionLabel + ":" + ppmMassError.ToString("F2"));

                        // append fragment ion intensity
                        fragmentIntensityStringBuilder.Append(ionLabel + ":" + ion.Intensity.ToString("F0"));

                        // append delimiter ", "
                        if (i < products.Count - 1)
                        {
                            stringBuilders.ForEach(p => p.Append(delimiter));
                        }
                    }

                    // append product type delimiter
                    stringBuilders.ForEach(p => p.Append("];"));
                }
            }

            // save ion series strings to output dictionary
            s[SpectrumMatchFromTsvHeader.MatchedIonSeries] = nullPsm ? " " : seriesStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonMzRatios] = nullPsm ? " " : mzStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa] = nullPsm ? " " : fragmentDaErrorStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonMassDiffPpm] = nullPsm ? " " : fragmentPpmErrorStringBuilder.ToString().TrimEnd(';');
            s[SpectrumMatchFromTsvHeader.MatchedIonIntensities] = nullPsm ? " " : fragmentIntensityStringBuilder.ToString().TrimEnd(';');

            // number of matched ions
            s[SpectrumMatchFromTsvHeader.MatchedIonCounts] = nullPsm ? " " : matchedIons.Count.ToString();
        }
    }
}
