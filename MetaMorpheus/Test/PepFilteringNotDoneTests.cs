using EngineLayer;
using EngineLayer.SpectrumMatch;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// PEP filtering requested (a PEP q-value threshold below the q-value threshold) but not possible: the filter says
    /// what it did instead, and never drops every match because PEP was not computed. Also the q-value notch filter a
    /// caller can ask for by name.
    /// </summary>
    [TestFixture]
    public static class PepFilteringNotDoneTests
    {
        // What the GUI sends when its PEP box is ticked.
        private static CommonParameters GuiPepBox => new CommonParameters(qValueThreshold: 1.0, pepQValueThreshold: 0.01);

        /// <summary>
        /// <paramref name="count"/> target matches; the first <paramref name="confident"/> have q, notch and PEP
        /// q-value 0.001, the rest 0.5. A PEP q-value of 2 is FdrInfo's "never computed".
        /// </summary>
        private static List<SpectralMatch> Matches(int count, int confident, bool pepComputed)
        {
            var peptide = new Protein("PEPTIDEK", "target1")
                .Digest(new DigestionParams(protease: "trypsin", minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First();
            var dataScan = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false),
                0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null,
                null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.HCD, 0, null);
            var scan = new Ms2ScanWithSpecificMass(dataScan, 2, 0, "File", new CommonParameters());

            var psms = new List<SpectralMatch>();
            for (int i = 0; i < count; i++)
            {
                var psm = new PeptideSpectralMatch(peptide, 0, 10, i, scan, new CommonParameters(), new List<MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                double q = i < confident ? 0.001 : 0.5;
                psm.SetFdrValues(i + 1, 0, q, i + 1, 0, q, pep: pepComputed ? q : 0, pepQValue: pepComputed ? q : 2);
                psms.Add(psm);
            }
            return psms;
        }

        [Test]
        public static void PepNotComputed_FiltersByQValueAndSaysSo_InsteadOfDroppingEveryMatch()
        {
            // 100 to 1000 matches: FdrAnalysisEngine trains PEP only above 1000 per digestion agent, so every PEP
            // q-value is still 2. Before, the PEP q-value filter was chosen anyway and kept none of them.
            var filtered = FilteredPsms.Filter(Matches(300, 200, pepComputed: false), GuiPepBox);

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.QValue));
            Assert.That(filtered.FilterThreshold, Is.EqualTo(0.01));
            Assert.That(filtered.FilteringNotPerformed, Is.False);
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(filtered.PepFilteringNotDone,
                Is.EqualTo("PEP filtering was requested, but PEP was not computed for these 300 PSMs, so results were filtered by q-value."));
        }

        [Test]
        public static void FewerThan100_IsNotFilteredAndSaysSo()
        {
            var filtered = FilteredPsms.Filter(Matches(50, 20, pepComputed: false), GuiPepBox);

            Assert.That(filtered.FilteringNotPerformed, Is.True);
            Assert.That(filtered.FilterThreshold, Is.EqualTo(1));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(50), "behaviour unchanged: nothing is filtered");
            Assert.That(filtered.PepFilteringNotDone,
                Is.EqualTo("PEP filtering was requested, but there were fewer than 100 PSMs (50), so results were not filtered."));
        }

        [Test]
        public static void PepComputed_FiltersOnPepAndSaysNothing()
        {
            var filtered = FilteredPsms.Filter(Matches(300, 150, pepComputed: true), GuiPepBox);

            Assert.That(filtered.FilterType, Is.EqualTo(FilterType.PepQValue));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(150));
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void DefaultSettings_KeepTheQValueFilterAndSayNothing()
        {
            var filtered = FilteredPsms.Filter(Matches(300, 200, pepComputed: false), new CommonParameters());

            Assert.That((filtered.FilterType, filtered.FilterThreshold, filtered.FilteringNotPerformed), Is.EqualTo((FilterType.QValue, 0.01, false)));
            Assert.That(filtered.FilteredPsmsList, Has.Count.EqualTo(200));
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void QValueNotch_ByName_FiltersOnTheNotchAlone()
        {
            // A notch q-value below the threshold and a whole-search q-value above it: the q-value filter (q AND notch)
            // drops the match, the notch filter keeps it.
            var psms = Matches(300, 0, pepComputed: false);
            foreach (var psm in psms.Take(120)) psm.SetFdrValues(1, 0, 0.5, 1, 0, 0.001, 0, 2);

            var byQValue = FilteredPsms.Filter(psms, new CommonParameters());
            var byNotch = FilteredPsms.Filter(psms, new CommonParameters(), filterType: FilterType.QValueNotch);

            Assert.That(byQValue.FilteredPsmsList, Is.Empty);
            Assert.That(byNotch.FilteredPsmsList, Has.Count.EqualTo(120));
            Assert.That((byNotch.FilterType, byNotch.FilterThreshold, byNotch.GetFilterTypeString()), Is.EqualTo((FilterType.QValueNotch, 0.01, "q-value notch")));
            Assert.That(byNotch.GetFilterValue(psms[0]), Is.EqualTo(0.001), "the value handed to quantification is the notch q-value");
            Assert.That(byNotch.TargetPsmsAboveThreshold, Is.EqualTo(120));
        }

        [Test]
        public static void QValueNotch_ByName_HasNoCountRule()
        {
            // The fewer-than-100 rule belongs to the inferred PEP filter; a filter asked for by name always filters.
            var filtered = FilteredPsms.Filter(Matches(50, 20, pepComputed: false), GuiPepBox, filterType: FilterType.QValueNotch);

            Assert.That(filtered.FilteringNotPerformed, Is.False);
            Assert.That(filtered.FilterThreshold, Is.EqualTo(1.0), "a filter asked for by name is held to its own threshold");
            Assert.That(filtered.PepFilteringNotDone, Is.Null);
        }

        [Test]
        public static void PepQValue_ByName_WithoutPepQValues_Throws()
        {
            var ex = Assert.Throws<MetaMorpheusException>(() =>
                FilteredPsms.Filter(Matches(300, 200, pepComputed: false), GuiPepBox, filterType: FilterType.PepQValue));
            Assert.That(ex!.Message, Does.Contain("none of these 300 PSMs has a PEP q-value"));

            // With PEP q-values it filters on them, at the PEP q-value threshold.
            var filtered = FilteredPsms.Filter(Matches(300, 150, pepComputed: true), GuiPepBox, filterType: FilterType.PepQValue);
            Assert.That((filtered.FilterType, filtered.FilterThreshold, filtered.FilteredPsmsList.Count), Is.EqualTo((FilterType.PepQValue, 0.01, 150)));
        }

        [Test]
        public static void GetFilterValue_FollowsTheFilterType()
        {
            var psm = Matches(1, 0, pepComputed: false)[0];
            psm.SetFdrValues(1, 0, 0.3, 1, 0, 0.2, 0.05, 0.1);

            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.QValue).GetFilterValue(psm), Is.EqualTo(0.3));
            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.QValueNotch).GetFilterValue(psm), Is.EqualTo(0.2));
            Assert.That(FilteredPsms.Filter(new[] { psm }, new CommonParameters(), filterType: FilterType.PepQValue).GetFilterValue(psm), Is.EqualTo(0.1));
        }
    }
}
