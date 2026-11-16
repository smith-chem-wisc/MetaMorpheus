#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// What a run's confident targets say of the library's iRT, by sequence composition. A predicted iRT is fixed for some other column, and on both benchmark runs what is left after the run's RT->iRT
/// curve has a shape: K, R and H elute 1.5 to 3 iRT earlier than predicted. Calibration fits this on its anchors, and the search expects each precursor where the fit says, in the window it reads as
/// well as in the error it scores. Fitted on 1,000 anchors, the median error of the targets a search scores too low goes from 2.47 to 1.59 iRT on PXD022589 HF-X, and from 2.22 to 1.25 on PXD005573 1 h.
/// </summary>
/// <remarks>
/// Only targets fit it. A permuted decoy keeps its target's composition and length, so those terms move it as its target; the terms on the first two residues read the decoy's own, as they would
/// any peptide's. Charge, oxidation and further length terms were tried and read the same out of sample.
/// </remarks>
public sealed class DiaIrtCompositionModel
{
    /// <summary>With fewer anchors the fit moves nothing: it has 25 terms, and calibration usually has 2,000 anchors.</summary>
    public const int MinimumAnchors = 200;

    /// <summary>How far a search looks past the run's iRT range, since a precursor's expected iRT can lie this far from its library value.</summary>
    public const double QueryMargin = 20;

    private const string Acids = "ACDEFGHIKLMNPQRSTVWY";
    private const string Basic = "KRH";
    private const string Acidic = "DE";
    private const double Penalty = 1;

    /// <summary>Anchors further than this many standard deviations from the first fit do not take part in the second: a calibration's anchors are confident, but some sit on the wrong peak.</summary>
    private const double OutlierSds = 3;

    private static readonly int TermCount = Acids.Length + 5;

    private readonly double[]? _weights;

    private DiaIrtCompositionModel(double[]? weights) => _weights = weights;

    /// <summary>A confident target.</summary>
    /// <param name="BaseSequence">Its plain sequence.</param>
    /// <param name="Residual">Its apex iRT on the run's curve, less its library iRT.</param>
    public sealed record Anchor(string BaseSequence, double Residual);

    /// <summary>How far from its library iRT this run's curve should find a precursor.</summary>
    public double Offset(string baseSequence)
    {
        if (_weights is null)
            return 0;
        var terms = new double[TermCount];
        FillTerms(baseSequence, terms);
        return Dot(_weights, terms);
    }

    /// <summary>A ridge fit of the residuals on composition, length and the first two residues, then again without the anchors the first fit could not explain.</summary>
    public static DiaIrtCompositionModel Fit(IReadOnlyCollection<Anchor> anchors)
    {
        ArgumentNullException.ThrowIfNull(anchors);
        if (anchors.Count < MinimumAnchors)
            return new DiaIrtCompositionModel(null);
        var rows = anchors.Select(a =>
        {
            var terms = new double[TermCount];
            FillTerms(a.BaseSequence, terms);
            return (Terms: terms, a.Residual);
        }).ToList();
        double[] weights = Solve(rows);
        double[] left = rows.Select(r => r.Residual - Dot(weights, r.Terms)).ToArray();
        double mean = left.Average();
        double sd = Math.Sqrt(left.Sum(v => (v - mean) * (v - mean)) / left.Length);
        var kept = rows.Where((_, i) => Math.Abs(left[i] - mean) <= OutlierSds * sd).ToList();
        return new DiaIrtCompositionModel(Solve(kept));
    }

    private static double Dot(double[] a, double[] b)
    {
        double sum = 0;
        for (int i = 0; i < a.Length; i++)
            sum += a[i] * b[i];
        return sum;
    }

    private static double[] Solve(List<(double[] Terms, double Residual)> rows)
    {
        var xtx = new double[TermCount, TermCount];
        var xty = new double[TermCount];
        foreach (var (terms, residual) in rows)
            for (int i = 0; i < TermCount; i++)
            {
                xty[i] += terms[i] * residual;
                for (int j = 0; j < TermCount; j++)
                    xtx[i, j] += terms[i] * terms[j];
            }
        for (int i = 0; i < TermCount; i++)
            xtx[i, i] += Penalty;
        return Matrix<double>.Build.DenseOfArray(xtx).Solve(Vector<double>.Build.DenseOfArray(xty)).ToArray();
    }

    private static void FillTerms(string baseSequence, double[] terms)
    {
        foreach (char residue in baseSequence)
        {
            int index = Acids.IndexOf(residue);
            if (index >= 0)
                terms[index]++;
        }
        int next = Acids.Length;
        terms[next++] = 10.0 / Math.Max(1, baseSequence.Length);
        terms[next++] = baseSequence.Length > 0 && Acidic.Contains(baseSequence[0]) ? 1 : 0;
        terms[next++] = baseSequence.Length > 0 && Basic.Contains(baseSequence[0]) ? 1 : 0;
        terms[next++] = baseSequence.Length > 1 && Basic.Contains(baseSequence[1]) ? 1 : 0;
        // The intercept
        terms[next] = 1;
    }
}
