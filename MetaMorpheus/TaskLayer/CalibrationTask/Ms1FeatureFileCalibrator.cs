using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    /// <summary>
    /// Gives the masses in an external MS1 feature file the same correction that calibration
    /// applied to the spectra.
    /// </summary>
    /// <remarks>
    /// A feature file is produced from the uncalibrated run, and mzLib's FromFile source takes its
    /// masses as given. Calibration multiplies every MS1 m/z by <c>1 - error</c>, with the error
    /// smoothed per scan. Without the same correction, the FromFile precursors in a calibrated file
    /// keep their pre-calibration offset, are searched at the tightened post-calibration tolerance,
    /// and can survive deduplication as a second copy of a classic precursor.
    ///
    /// Each feature is corrected by the factor of the MS1 scan nearest its apex. The correction is
    /// relative, and neutral mass is linear in m/z, so one factor serves every charge state of the
    /// feature (to within the proton term, about z x 1.007 x error, far below any tolerance).
    /// </remarks>
    public static class Ms1FeatureFileCalibrator
    {
        /// <summary>True when <paramref name="path"/> is an <c>_ms1.feature</c> file, the only feature format this can rewrite.</summary>
        public static bool CanCalibrate(string path) =>
            !string.IsNullOrWhiteSpace(path)
            && path.EndsWith(SupportedFileType.Ms1Feature.GetFileExtension(), StringComparison.OrdinalIgnoreCase);

        /// <summary>
        /// Combines successive calibration rounds into one multiplicative factor per MS1 scan.
        /// Every round must cover the same MS1 scans in the same order, which holds because each
        /// round recalibrates the previous round's output scan for scan.
        /// </summary>
        public static (double[] RetentionTimes, double[] Factors) Combine(
            IReadOnlyList<IReadOnlyList<(double RetentionTime, double RelativeError)>> rounds)
        {
            if (rounds == null || rounds.Count == 0)
                return (Array.Empty<double>(), Array.Empty<double>());

            double[] rts = rounds[0].Select(c => c.RetentionTime).ToArray();
            double[] factors = Enumerable.Repeat(1.0, rts.Length).ToArray();
            foreach (var round in rounds)
            {
                if (round.Count != rts.Length)
                    throw new ArgumentException("calibration rounds cover different numbers of MS1 scans", nameof(rounds));
                for (int i = 0; i < rts.Length; i++)
                    factors[i] *= 1 - round[i].RelativeError;
            }
            return (rts, factors);
        }

        /// <summary>The factor of the MS1 scan nearest <paramref name="retentionTime"/> (minutes); 1 when there are no scans.</summary>
        public static double FactorAt(double[] retentionTimes, double[] factors, double retentionTime)
        {
            if (retentionTimes.Length == 0)
                return 1.0;
            int idx = Array.BinarySearch(retentionTimes, retentionTime);
            if (idx >= 0)
                return factors[idx];
            idx = ~idx;
            if (idx == 0) return factors[0];
            if (idx == retentionTimes.Length) return factors[^1];
            return retentionTime - retentionTimes[idx - 1] <= retentionTimes[idx] - retentionTime
                ? factors[idx - 1]
                : factors[idx];
        }

        /// <summary>
        /// Writes a copy of <paramref name="sourcePath"/> to <paramref name="destinationPath"/> with
        /// every feature's mass multiplied by the combined calibration factor at its apex.
        /// </summary>
        /// <param name="rtInSeconds">
        /// The source's retention times are in seconds (mzLib detects and records this on load, see
        /// <c>FromFileDeconvolutionParameters.RetentionTimeNormalizedFromSeconds</c>). The written
        /// file keeps the source's units; only the apex lookup converts them.
        /// </param>
        /// <returns>The number of features written.</returns>
        public static int WriteCalibratedCopy(string sourcePath, string destinationPath,
            IReadOnlyList<IReadOnlyList<(double RetentionTime, double RelativeError)>> rounds, bool rtInSeconds)
        {
            var (rts, factors) = Combine(rounds);
            var source = new Ms1FeatureFile(sourcePath);
            List<Ms1Feature> rows = source.Results;
            double toMinutes = rtInSeconds ? 1.0 / 60.0 : 1.0;
            foreach (var row in rows)
                row.Mass *= FactorAt(rts, factors, row.RetentionTimeApex * toMinutes);

            Directory.CreateDirectory(Path.GetDirectoryName(Path.GetFullPath(destinationPath)));
            var calibrated = new Ms1FeatureFile { Results = rows, Software = source.Software };
            calibrated.WriteResults(destinationPath);
            return rows.Count;
        }
    }
}
