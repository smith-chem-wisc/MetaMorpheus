using InternalLogicEngineLayer;
using NUnit.Framework;
using System;
using System.Text;
using System.Collections.Generic;
using OldInternalLogic;
using MassSpectrometry;
using Spectra;

namespace Test
{
    [TestFixture]
    public class GPTMDengineTest
    {

		[Test]
		public void TestGPTMDengine()
		{

			List<NewPsmWithFDR> allResultingIdentifications = null;
			List<MorpheusModification> gptmdModifications = new List<MorpheusModification>() { new MorpheusModification("name", ModificationType.AminoAcidResidue, 'N', 42, null, null, '\0', double.NaN, false, null) };
			IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
			double tol = 0.1;
			bool isotopeErrors = false;
			var engine = new GPTMDEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
			Assert.That(() => engine.Run(),
			Throws.TypeOf<EngineValidationException>()
				.With.Property("Message").EqualTo("allResultingIdentifications cannot be null"));


			allResultingIdentifications = new List<NewPsmWithFDR>();
			engine = new GPTMDEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
			var res = (GPTMDResults)engine.Run();
			Assert.AreEqual(0, res.mods.Count);


			ParentSpectrumMatch newPsm = new TestParentSpectrumMatch(588.22520189093 + 42);
			Protein parentProtein = new Protein("NNNNN", "accession", null, new Dictionary<int, List<MorpheusModification>>(), null, null, null, null, null, 0, false);
			PeptideWithPossibleModifications modPep = new PeptideWithPossibleModifications(1,5, parentProtein,0,"ugh");
			var twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
			HashSet<PeptideWithSetModifications> peptidesWithSetModifications = new HashSet<PeptideWithSetModifications>(){new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss)};
			Tolerance fragmentTolerance = null;
			IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = null;
			var thisPSM = new PSMwithTargetDecoyKnown(newPsm, peptidesWithSetModifications, fragmentTolerance, myMsDataFile);
			var newPsmWithFDR = new NewPsmWithFDR(thisPSM, 1, 0, 0);
			allResultingIdentifications.Add(newPsmWithFDR);

			engine = new GPTMDEngine(allResultingIdentifications, isotopeErrors, gptmdModifications, combos, tol);
			res = (GPTMDResults)engine.Run();
			Assert.AreEqual(1, res.mods.Count);
			Assert.AreEqual(5, res.mods["accession"].Count);

        }

		class TestParentSpectrumMatch : ParentSpectrumMatch
		{
			public override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
			{
				throw new NotImplementedException();
			}
			public TestParentSpectrumMatch(double scanPrecursorMass)
			{
				this.scanPrecursorMass= scanPrecursorMass;
			}
		}
	}
}