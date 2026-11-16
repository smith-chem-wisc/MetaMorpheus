#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Test.DiaLibrarySearch;

/// <summary>
/// A synthetic DIA run and the spectral library it was drawn from, for testing the DIA library search end to end
/// with a known answer.
/// <para>
/// Every target has a reversed-sequence decoy with the same precursor m/z and iRT but different fragment m/z, as a
/// real reversed decoy would have. The run holds only MS2 scans: one per isolation window per cycle. The chosen
/// ("planted") precursors elute as Gaussians centred on <see cref="TrueRtMinutes"/> of their library iRT. That map is
/// deliberately nonlinear, so a search that treats iRT as minutes looks in the wrong place. Every scan also carries
/// seeded random noise peaks.
/// </para>
/// </summary>
[ExcludeFromCodeCoverage]
internal sealed class SyntheticDiaRun
{
    public const double FirstWindowLowMz = 400;
    public const double WindowWidth = 40;
    public const int WindowCount = 10;
    public const double CycleMinutes = 0.05;
    public const double RunMinutes = 10;
    public const double ElutionSigmaMinutes = 0.06;
    public const int FragmentsPerPrecursor = 6;

    public List<MslLibraryEntry> Library { get; }
    public MsDataScan[] Scans { get; }

    /// <summary>Full sequences of the library entries (targets or decoys) whose fragments were planted in the run.</summary>
    public IReadOnlySet<string> PlantedSequences { get; }

    /// <summary>
    /// The true run retention time, in minutes, of a library iRT. It is nonlinear and increasing over the library's
    /// iRT range [-20, 120], and maps that range onto roughly [0.5, 9.9] minutes.
    /// </summary>
    public static double TrueRtMinutes(double irt) => 0.5 + 0.045 * (irt + 20) + 0.00016 * (irt + 20) * (irt + 20);

    /// <summary>A stable bucket for choosing a fraction of entries. string.GetHashCode is randomized per process.</summary>
    public static int Bucket(MslLibraryEntry entry, int buckets) => Bucket(entry.FullSequence, buckets);

    /// <inheritdoc cref="Bucket(MslLibraryEntry, int)"/>
    public static int Bucket(string fullSequence, int buckets) => fullSequence.Sum(c => c) % buckets;

    private SyntheticDiaRun(List<MslLibraryEntry> library, MsDataScan[] scans, IReadOnlySet<string> planted)
    {
        Library = library;
        Scans = scans;
        PlantedSequences = planted;
    }

    /// <param name="targetCount">Number of target precursors in the library. Each also gets one decoy.</param>
    /// <param name="plant">Chooses which library entries elute in the run, by full sequence.</param>
    /// <param name="withDecoys">False builds a library with no decoys at all.</param>
    /// <param name="abundance">Multiplies a planted entry's elution peak; null leaves every peak at the same height.</param>
    /// <param name="withMs1">Adds an MS1 scan at the start of each cycle: noise, and each planted precursor's M0-M3 isotopes
    /// (Poisson, lambda = neutral mass / 1800, about averagine) on the same elution profile.</param>
    /// <param name="withCharge3Siblings">Adds each peptide's 3+ precursor (same fragments, same iRT) where its m/z is inside
    /// the run's windows, for targets and decoys alike. It elutes with its 2+ sibling when planted.</param>
    /// <param name="ppmOffset">Shifts every recorded m/z by this many ppm, as a miscalibrated instrument would.</param>
    public static SyntheticDiaRun Build(int targetCount, Func<MslLibraryEntry, bool> plant, bool withDecoys = true,
        int seed = 42, double noisePeaksPerScan = 60, Func<MslLibraryEntry, double>? abundance = null, bool withMs1 = false, double ppmOffset = 0, bool withCharge3Siblings = false)
    {
        var random = new Random(seed);
        var library = new List<MslLibraryEntry>();
        const string residues = "ACDEFGHILMNPQSTVWY";

        for (int i = 0; i < targetCount; i++)
        {
            // A distinct, non-palindromic tryptic-looking sequence, so its reverse is a different sequence
            string core = new string(Enumerable.Range(0, 9).Select(_ => residues[random.Next(residues.Length)]).ToArray());
            string sequence = $"{core}{i:D3}K".Replace("0", "G").Replace("1", "A").Replace("2", "S").Replace("3", "T")
                .Replace("4", "V").Replace("5", "L").Replace("6", "N").Replace("7", "D").Replace("8", "E").Replace("9", "Q");
            double precursorMz = FirstWindowLowMz + 1 + random.NextDouble() * (WindowCount * WindowWidth - 2);
            double irt = -15 + random.NextDouble() * 130;

            library.Add(Entry(sequence, precursorMz, irt, isDecoy: false, random));
            if (withDecoys)
                library.Add(Entry(new string(sequence.Reverse().ToArray()), precursorMz, irt, isDecoy: true, random));
            if (withCharge3Siblings)
                foreach (var twoPlus in library.Skip(library.Count - (withDecoys ? 2 : 1)).ToList())
                    if (Charge3Sibling(twoPlus) is { } threePlus)
                        library.Add(threePlus);
        }

        var planted = library.Where(plant).ToList();
        var scans = new List<MsDataScan>();
        int scanNumber = 1;
        for (double rt = 0; rt < RunMinutes; rt += CycleMinutes)
        {
            if (withMs1)
            {
                var ms1Peaks = new List<(double Mz, double Intensity)>();
                for (int n = 0; n < 200; n++)
                    ms1Peaks.Add((FirstWindowLowMz + random.NextDouble() * WindowCount * WindowWidth, 50 + random.NextDouble() * 950));
                foreach (var entry in planted)
                {
                    double elution = 1e6 * (abundance?.Invoke(entry) ?? 1) * Math.Exp(-0.5 * Math.Pow((rt - TrueRtMinutes(entry.RetentionTime)) / ElutionSigmaMinutes, 2));
                    if (elution < 1)
                        continue;
                    double lambda = (entry.PrecursorMz - 1.007276) * entry.ChargeState / 1800;
                    double p = Math.Exp(-lambda);
                    for (int k = 0; k <= 3; k++, p *= lambda / k)
                        ms1Peaks.Add((entry.PrecursorMz + k * 1.0033548 / entry.ChargeState, elution * p));
                }
                var ms1Ordered = ms1Peaks.OrderBy(p => p.Mz).ToArray();
                var ms1Spectrum = new MzSpectrum(ms1Ordered.Select(p => p.Mz * (1 + ppmOffset * 1e-6)).ToArray(), ms1Ordered.Select(p => p.Intensity).ToArray(), false);
                scans.Add(new MsDataScan(ms1Spectrum, scanNumber, 1, true, Polarity.Positive, rt, new MzRange(FirstWindowLowMz, FirstWindowLowMz + WindowCount * WindowWidth),
                    "synthetic", MZAnalyzerType.Orbitrap, ms1Spectrum.SumOfAllY, 10, null, $"scan={scanNumber}"));
                scanNumber++;
            }
            for (int w = 0; w < WindowCount; w++)
            {
                double low = FirstWindowLowMz + w * WindowWidth;
                var peaks = new List<(double Mz, double Intensity)>();
                for (int n = 0; n < noisePeaksPerScan; n++)
                    peaks.Add((150 + random.NextDouble() * 1350, 50 + random.NextDouble() * 950));

                foreach (var entry in planted.Where(e => e.PrecursorMz >= low && e.PrecursorMz < low + WindowWidth))
                {
                    double elution = 1e5 * (abundance?.Invoke(entry) ?? 1) * Math.Exp(-0.5 * Math.Pow((rt - TrueRtMinutes(entry.RetentionTime)) / ElutionSigmaMinutes, 2));
                    if (elution < 1)
                        continue;
                    foreach (var fragment in entry.MatchedFragmentIons)
                        peaks.Add((fragment.Mz, elution * fragment.Intensity));
                }

                var ordered = peaks.OrderBy(p => p.Mz).ToArray();
                var spectrum = new MzSpectrum(ordered.Select(p => p.Mz * (1 + ppmOffset * 1e-6)).ToArray(), ordered.Select(p => p.Intensity).ToArray(), false);
                scans.Add(Scan(spectrum, scanNumber, rt, isolationMz: low + WindowWidth / 2));
                scanNumber++;
            }
        }

        return new SyntheticDiaRun(library, scans.ToArray(), planted.Select(e => e.FullSequence).ToHashSet());
    }

    /// <summary>
    /// A scan as the synthetic run builds them. An MS2 scan isolates <see cref="WindowWidth"/> around
    /// <paramref name="isolationMz"/>; with no isolation m/z it has no isolation range.
    /// </summary>
    public static MsDataScan Scan(MzSpectrum spectrum, int scanNumber, double rtMinutes, double? isolationMz, int msnOrder = 2) =>
        new(spectrum, scanNumber, msnOrder, true, Polarity.Positive, rtMinutes, new MzRange(150, 1500), "synthetic",
            MZAnalyzerType.Orbitrap, spectrum.SumOfAllY, 10, null, $"scan={scanNumber}",
            isolationMZ: isolationMz, isolationWidth: isolationMz is null ? null : WindowWidth,
            dissociationType: msnOrder == 2 ? DissociationType.HCD : null);

    /// <summary>Writes the library as a .msl, validated first, and returns its path.</summary>
    public string WriteLibrary(string directory)
    {
        var problems = MslWriter.ValidateEntries(Library);
        if (problems.Count > 0)
            throw new InvalidOperationException("Synthetic library is invalid: " + string.Join("; ", problems));
        string path = Path.Combine(directory, $"synthetic_{Guid.NewGuid():N}.msl");
        MslLibrary.Save(path, Library);
        return path;
    }

    /// <summary>The same peptide at 3+: its m/z from the 2+ m/z, the same fragments, or null when outside the run's windows.</summary>
    private static MslLibraryEntry? Charge3Sibling(MslLibraryEntry twoPlus)
    {
        const double proton = 1.007276;
        double mz = ((twoPlus.PrecursorMz - proton) * 2 + 3 * proton) / 3;
        if (mz < FirstWindowLowMz || mz >= FirstWindowLowMz + WindowCount * WindowWidth)
            return null;
        var sibling = new MslLibraryEntry
        {
            FullSequence = twoPlus.FullSequence,
            BaseSequence = twoPlus.BaseSequence,
            PrecursorMz = mz,
            ChargeState = 3,
            RetentionTime = twoPlus.RetentionTime,
            IsDecoy = twoPlus.IsDecoy,
        };
        foreach (var f in twoPlus.MatchedFragmentIons)
            sibling.MatchedFragmentIons.Add(new MslFragmentIon
            {
                Mz = f.Mz, Intensity = f.Intensity, ProductType = f.ProductType, FragmentNumber = f.FragmentNumber,
                ResiduePosition = f.ResiduePosition, Charge = f.Charge,
            });
        return sibling;
    }

    private static MslLibraryEntry Entry(string sequence, double precursorMz, double irt, bool isDecoy, Random random)
    {
        var entry = new MslLibraryEntry
        {
            FullSequence = sequence,
            BaseSequence = sequence,
            PrecursorMz = precursorMz,
            ChargeState = 2,
            RetentionTime = irt,
            IsDecoy = isDecoy,
        };
        for (int f = 0; f < FragmentsPerPrecursor; f++)
        {
            entry.MatchedFragmentIons.Add(new MslFragmentIon
            {
                Mz = (float)(200 + random.NextDouble() * 1000),
                Intensity = (float)(0.1 + random.NextDouble() * 0.9),
                ProductType = ProductType.y,
                FragmentNumber = f + 2,
                ResiduePosition = sequence.Length - f - 2,
                Charge = 1,
            });
        }
        return entry;
    }
}
