// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// REPO:     mzLib
// LOCATION: mzLib/Test/Dia/DiaSequenceParserTests.cs

using MassSpectrometry.Dia;
using NUnit.Framework;
using Org.BouncyCastle.Asn1;

namespace Test.Dia
{
    /// <summary>
    /// Unit tests for <see cref="DiaSequenceParser"/>.
    /// All tests are self-contained with no external data files.
    /// </summary>
    [TestFixture]
    public class DiaSequenceParserTests
    {
        [Test]
        public void GetBaseSequence_NoMods_ReturnsIdentical()
        {
            const string seq = "PEPTIDEK";
            Assert.That(DiaSequenceParser.GetBaseSequence(seq), Is.EqualTo(seq));
        }

        [Test]
        public void GetBaseSequence_SingleMod_StripsCorrectly()
        {
            Assert.That(DiaSequenceParser.GetBaseSequence("PEPTM[Oxidation]IDEK"),
                Is.EqualTo("PEPTMIDEK"));
        }

        [Test]
        public void GetBaseSequence_TwoMods_BothStripped()
        {
            Assert.That(DiaSequenceParser.GetBaseSequence("M[Oxidation]PEPTM[Oxidation]"),
                Is.EqualTo("MPEPTM"));
        }

        [Test]
        public void Parse_SingleMod_CorrectResidueIndex()
        {
            DiaSequenceParser.Parse("PEPTM[Oxidation]IDEK",
                out string baseSeq,
                out var mods);

            Assert.That(baseSeq, Is.EqualTo("PEPTMIDEK"));
            Assert.That(mods, Has.Count.EqualTo(1));
            Assert.That(mods[0].Position, Is.EqualTo(5));
            Assert.That(mods[0].Name, Is.EqualTo("Oxidation"));
        }

        [Test]
        public void Parse_NTerminalMod_IndexIsZero()
        {
            DiaSequenceParser.Parse("[Acetyl]PEPTIDE",
                out string baseSeq,
                out var mods);

            Assert.That(baseSeq, Is.EqualTo("PEPTIDE"));
            Assert.That(mods, Has.Count.EqualTo(1));
            Assert.That(mods[0].Position, Is.EqualTo(0));
            Assert.That(mods[0].Name, Is.EqualTo("Acetyl"));
        }

        [Test]
        public void FormatModificationSummary_SingleMod_CorrectString()
        {
            string summary = DiaSequenceParser.FormatModificationSummary("PEPTM[Oxidation]IDEK");
            Assert.That(summary, Is.EqualTo("Oxidation@5"));
        }

        [Test]
        public void FormatModificationSummary_NoMods_ReturnsEmpty()
        {
            string summary = DiaSequenceParser.FormatModificationSummary("PEPTIDEK");
            Assert.That(summary, Is.EqualTo(string.Empty));
        }
    }
}