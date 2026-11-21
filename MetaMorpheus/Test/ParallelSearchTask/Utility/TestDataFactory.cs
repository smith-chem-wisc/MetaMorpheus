using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer.ParallelSearch.Analysis;

namespace Test.ParallelSearchTask.Utility;

/// <summary>
/// Factory for creating test data for statistical analysis testing
/// Provides methods to create mock AggregatedAnalysisResults with controlled distributions
/// </summary>
public static class TestDataFactory
{
    /// <summary>
    /// Random number generator with fixed seed for reproducibility
    /// </summary>
    private static readonly Random Random = new Random(42);

    #region Basic Result Creation

    /// <summary>
    /// Create a basic result with specified counts
    /// </summary>
    public static TransientDatabaseMetrics CreateBasicResult(
        string databaseName,
        int psmCount = 0,
        int peptideCount = 0,
        int proteinGroupCount = 0)
    {
        return new TransientDatabaseMetrics
        {
            DatabaseName = databaseName,
            TargetPsmsFromTransientDbAtQValueThreshold = psmCount,
            TargetPeptidesFromTransientDbAtQValueThreshold = peptideCount,
            TargetProteinGroupsFromTransientDbAtQValueThreshold = proteinGroupCount,
            TransientProteinCount = 100,
            TransientPeptideCount = 250,
            PsmBacterialUnambiguousTargets = psmCount / 2,
            PsmBacterialUnambiguousDecoys = 5,
            PeptideBacterialUnambiguousTargets = peptideCount / 2,
            PeptideBacterialUnambiguousDecoys = 3,
            ProteinGroupBacterialUnambiguousTargets = proteinGroupCount / 2,
            ProteinGroupBacterialUnambiguousDecoys = 2
        };
    }

    /// <summary>
    /// Create a result with realistic organism-specific counts and score distributions
    /// </summary>
    public static TransientDatabaseMetrics CreateRealisticResult(
        string databaseName,
        int targetPsms,
        int decoyPsms,
        int targetPeptides,
        int decoyPeptides,
        int targetProteinGroups = 0,
        int decoyProteinGroups = 0,
        bool includeScores = true,
        double meanTargetScore = 100.0,
        double meanDecoyScore = 50.0, 
        int transientProteinCount = 100,
        int transientPeptideCount = 250)
    {
        var result = new TransientDatabaseMetrics
        {
            DatabaseName = databaseName,
            TargetPsmsFromTransientDbAtQValueThreshold = targetPsms,
            TargetPeptidesFromTransientDbAtQValueThreshold = targetPeptides,
            TargetProteinGroupsFromTransientDbAtQValueThreshold = targetProteinGroups,
            TransientProteinCount = transientProteinCount,
            TransientPeptideCount = transientPeptideCount,

            // Organism-specific counts
            PsmBacterialUnambiguousTargets = targetPsms / 2,
            PsmBacterialUnambiguousDecoys = decoyPsms / 2,
            PeptideBacterialUnambiguousTargets = targetPeptides / 2,
            PeptideBacterialUnambiguousDecoys = decoyPeptides / 2,
            ProteinGroupBacterialUnambiguousTargets = targetProteinGroups / 2,
            ProteinGroupBacterialUnambiguousDecoys = decoyProteinGroups / 2
        };

        if (includeScores)
        {
            result.PsmBacterialUnambiguousTargetScores = GenerateScores(targetPsms, meanTargetScore, 20.0);
            result.PsmBacterialUnambiguousDecoyScores = GenerateScores(decoyPsms, meanDecoyScore, 15.0);
            result.PeptideBacterialUnambiguousTargetScores = GenerateScores(targetPeptides, meanTargetScore, 20.0);
            result.PeptideBacterialUnambiguousDecoyScores = GenerateScores(decoyPeptides, meanDecoyScore, 15.0);
        }

        return result;
    }

    #endregion

    #region Scenario-Based Result Creation

    /// <summary>
    /// Create a collection of results with known statistical properties
    /// Suitable for testing that high-signal organisms are detected
    /// </summary>
    public static List<TransientDatabaseMetrics> CreateHighSignalScenario(int count = 50)
    {
        var results = new List<TransientDatabaseMetrics>();

        // Create mostly low-signal databases
        for (int i = 0; i < count - 5; i++)
        {
            results.Add(CreateRealisticResult(
                databaseName: $"LowSignal_{i}",
                targetPsms: Random.Next(10, 30),
                decoyPsms: Random.Next(8, 25),
                targetPeptides: Random.Next(5, 15),
                decoyPeptides: Random.Next(4, 12),
                targetProteinGroups: Random.Next(2, 8),
                decoyProteinGroups: Random.Next(1, 6)
            ));
        }

        // Add a few high-signal databases that should be detected
        for (int i = 0; i < 5; i++)
        {
            results.Add(CreateRealisticResult(
                databaseName: $"HighSignal_{i}",
                targetPsms: Random.Next(200, 500),
                decoyPsms: Random.Next(5, 15),
                targetPeptides: Random.Next(100, 250),
                decoyPeptides: Random.Next(2, 8),
                targetProteinGroups: Random.Next(50, 100),
                decoyProteinGroups: Random.Next(1, 3),
                meanTargetScore: 150.0,
                meanDecoyScore: 50.0
            ));
        }

        return results;
    }

    /// <summary>
    /// Create results following a Gaussian distribution
    /// </summary>
    public static List<TransientDatabaseMetrics> CreateGaussianDistribution(
        int count = 100,
        double mean = 50,
        double stdDev = 10)
    {
        var results = new List<TransientDatabaseMetrics>();

        for (int i = 0; i < count; i++)
        {
            int psmCount = Math.Max(0, (int)(mean + (GaussianRandom() * stdDev)));
            results.Add(CreateBasicResult(
                databaseName: $"Organism_{i:D3}",
                psmCount: psmCount,
                peptideCount: psmCount / 2,
                proteinGroupCount: psmCount / 5
            ));
        }

        return results;
    }

    /// <summary>
    /// Create results with overdispersion (variance > mean)
    /// Suitable for testing Negative Binomial distribution
    /// </summary>
    public static List<TransientDatabaseMetrics> CreateOverdispersedDistribution(
        int count = 100,
        double mean = 20,
        double overdispersionFactor = 2.0)
    {
        var results = new List<TransientDatabaseMetrics>();
        double variance = mean * overdispersionFactor;

        for (int i = 0; i < count; i++)
        {
            // Use negative binomial sampling
            int psmCount = SampleNegativeBinomial(mean, variance);
            results.Add(CreateBasicResult(
                databaseName: $"Organism_{i:D3}",
                psmCount: psmCount,
                peptideCount: psmCount / 2,
                proteinGroupCount: psmCount / 5
            ));
        }

        return results;
    }

    /// <summary>
    /// Create null scenario where all databases have similar low signal
    /// None should be statistically significant
    /// </summary>
    public static List<TransientDatabaseMetrics> CreateNullScenario(int count = 50)
    {
        var results = new List<TransientDatabaseMetrics>();

        for (int i = 0; i < count; i++)
        {
            results.Add(CreateRealisticResult(
                databaseName: $"NullOrganism_{i}",
                targetPsms: Random.Next(8, 15),
                decoyPsms: Random.Next(7, 14),
                targetPeptides: Random.Next(4, 8),
                decoyPeptides: Random.Next(3, 7),
                targetProteinGroups: Random.Next(1, 4),
                decoyProteinGroups: Random.Next(1, 3),
                meanTargetScore: 60.0,
                meanDecoyScore: 55.0
            ));
        }

        return results;
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Create edge case scenarios for testing robustness
    /// </summary>
    public static class EdgeCases
    {
        /// <summary>
        /// Single observation
        /// </summary>
        public static List<TransientDatabaseMetrics> SingleObservation()
        {
            return new List<TransientDatabaseMetrics>
            {
                CreateBasicResult("OnlyOne", 10, 5, 2)
            };
        }

        /// <summary>
        /// All zeros
        /// </summary>
        public static List<TransientDatabaseMetrics> AllZeros(int count = 10)
        {
            var results = new List<TransientDatabaseMetrics>();
            for (int i = 0; i < count; i++)
            {
                results.Add(CreateBasicResult($"Zero_{i}", 0, 0, 0));
            }
            return results;
        }

        /// <summary>
        /// All same value
        /// </summary>
        public static List<TransientDatabaseMetrics> AllSameValue(int count = 10, int value = 20)
        {
            var results = new List<TransientDatabaseMetrics>();
            for (int i = 0; i < count; i++)
            {
                results.Add(CreateBasicResult($"Same_{i}", value, value / 2, value / 5));
            }
            return results;
        }

        /// <summary>
        /// Single outlier
        /// </summary>
        public static List<TransientDatabaseMetrics> SingleOutlier(int count = 50)
        {
            var results = CreateGaussianDistribution(count - 1, 20, 5);
            results.Add(CreateBasicResult("Outlier", 500, 250, 100));
            return results;
        }

        /// <summary>
        /// Empty list
        /// </summary>
        public static List<TransientDatabaseMetrics> Empty()
        {
            return new List<TransientDatabaseMetrics>();
        }

        /// <summary>
        /// Very large numbers (test numerical stability)
        /// </summary>
        public static List<TransientDatabaseMetrics> LargeNumbers(int count = 10)
        {
            var results = new List<TransientDatabaseMetrics>();
            for (int i = 0; i < count; i++)
            {
                results.Add(CreateBasicResult($"Large_{i}", 1000000, 500000, 100000));
            }
            return results;
        }

        /// <summary>
        /// Mix of very small and very large values
        /// </summary>
        public static List<TransientDatabaseMetrics> MixedScale(int count = 20)
        {
            var results = new List<TransientDatabaseMetrics>();
            for (int i = 0; i < count / 2; i++)
            {
                results.Add(CreateBasicResult($"Small_{i}", 1, 0, 0));
            }
            for (int i = 0; i < count / 2; i++)
            {
                results.Add(CreateBasicResult($"Large_{i}", 10000, 5000, 1000));
            }
            return results;
        }
    }

    #endregion

    #region Helper Methods

    /// <summary>
    /// Generate normally distributed scores
    /// </summary>
    private static double[] GenerateScores(int count, double mean, double stdDev)
    {
        if (count == 0) return Array.Empty<double>();

        var scores = new double[count];
        for (int i = 0; i < count; i++)
        {
            scores[i] = Math.Max(0, mean + (GaussianRandom() * stdDev));
        }
        return scores;
    }

    /// <summary>
    /// Box-Muller transform for Gaussian random numbers
    /// </summary>
    private static double GaussianRandom()
    {
        double u1 = 1.0 - Random.NextDouble();
        double u2 = 1.0 - Random.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
    }

    /// <summary>
    /// Sample from negative binomial distribution using gamma-Poisson mixture
    /// </summary>
    private static int SampleNegativeBinomial(double mean, double variance)
    {
        if (variance <= mean)
        {
            // Fallback to Poisson if not overdispersed
            return SamplePoisson(mean);
        }

        // Parameterize using mean and variance
        double r = (mean * mean) / (variance - mean);
        double p = mean / variance;

        // Use gamma-Poisson mixture
        double lambda = SampleGamma(r, (1 - p) / p);
        return SamplePoisson(lambda);
    }

    /// <summary>
    /// Sample from Poisson distribution
    /// </summary>
    private static int SamplePoisson(double lambda)
    {
        if (lambda < 30)
        {
            // Use Knuth's algorithm for small lambda
            double L = Math.Exp(-lambda);
            int k = 0;
            double p = 1.0;

            do
            {
                k++;
                p *= Random.NextDouble();
            } while (p > L);

            return k - 1;
        }
        else
        {
            // Use normal approximation for large lambda
            return Math.Max(0, (int)(lambda + Math.Sqrt(lambda) * GaussianRandom()));
        }
    }

    /// <summary>
    /// Sample from Gamma distribution using Marsaglia and Tsang's method
    /// </summary>
    private static double SampleGamma(double shape, double scale)
    {
        if (shape < 1)
        {
            // Use Johnk's generator for shape < 1
            double u, v, x, y;
            do
            {
                u = Random.NextDouble();
                v = Random.NextDouble();
                x = Math.Pow(u, 1.0 / shape);
                y = Math.Pow(v, 1.0 / (1.0 - shape));
            } while (x + y > 1);

            double e = -Math.Log(Random.NextDouble());
            return scale * e * x / (x + y);
        }
        else
        {
            // Marsaglia and Tsang's method
            double d = shape - 1.0 / 3.0;
            double c = 1.0 / Math.Sqrt(9.0 * d);

            while (true)
            {
                double x = GaussianRandom();
                double v = 1.0 + c * x;

                if (v <= 0) continue;

                v = v * v * v;
                double u = Random.NextDouble();
                double x2 = x * x;

                if (u < 1.0 - 0.0331 * x2 * x2)
                    return scale * d * v;

                if (Math.Log(u) < 0.5 * x2 + d * (1.0 - v + Math.Log(v)))
                    return scale * d * v;
            }
        }
    }

    #endregion

    #region Fisher Exact Test Data

    /// <summary>
    /// Create data suitable for Fisher's Exact Test
    /// Returns results with controlled ambiguous vs unambiguous ratios
    /// </summary>
    public static List<TransientDatabaseMetrics> CreateFisherTestData(
        int nullCount = 40,
        int enrichedCount = 10)
    {
        var results = new List<TransientDatabaseMetrics>();

        // Null organisms: similar ratios of ambiguous/unambiguous
        for (int i = 0; i < nullCount; i++)
        {
            results.Add(new TransientDatabaseMetrics
            {
                DatabaseName = $"Null_{i}",
                PsmBacterialUnambiguousTargets = 10,
                PsmBacterialUnambiguousDecoys = 2,
                PsmBacterialAmbiguous = 15,
                TransientProteinCount = 100
            });
        }

        // Enriched organisms: high ratio of unambiguous
        for (int i = 0; i < enrichedCount; i++)
        {
            results.Add(new TransientDatabaseMetrics
            {
                DatabaseName = $"Enriched_{i}",
                PsmBacterialUnambiguousTargets = 50,
                PsmBacterialUnambiguousDecoys = 2,
                PsmBacterialAmbiguous = 5,
                TransientProteinCount = 100
            });
        }

        return results;
    }

    #endregion

    #region Kolmogorov-Smirnov Test Data

    /// <summary>
    /// Create data with shifted score distributions for K-S test
    /// </summary>
    public static List<TransientDatabaseMetrics> CreateKSTestData(
        int nullCount = 40,
        int shiftedCount = 10,
        double nullMean = 50.0,
        double shiftedMean = 100.0)
    {
        var results = new List<TransientDatabaseMetrics>();

        // Null: target and decoy scores similar
        for (int i = 0; i < nullCount; i++)
        {
            results.Add(new TransientDatabaseMetrics
            {
                DatabaseName = $"Null_{i}",
                PsmBacterialUnambiguousTargetScores = GenerateScores(20, nullMean, 15.0),
                PsmBacterialUnambiguousDecoyScores = GenerateScores(20, nullMean - 5.0, 15.0),
                TransientProteinCount = 100
            });
        }

        // Shifted: target scores significantly higher
        for (int i = 0; i < shiftedCount; i++)
        {
            results.Add(new TransientDatabaseMetrics
            {
                DatabaseName = $"Shifted_{i}",
                PsmBacterialUnambiguousTargetScores = GenerateScores(30, shiftedMean, 20.0),
                PsmBacterialUnambiguousDecoyScores = GenerateScores(15, nullMean, 15.0),
                TransientProteinCount = 100
            });
        }

        return results;
    }

    #endregion
}
