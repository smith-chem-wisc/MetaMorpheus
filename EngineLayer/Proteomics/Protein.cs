using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class Protein
    {

        #region Private Fields

        private string fullDescription;

        #endregion Private Fields

        #region Public Constructors

        public Protein(string baseSequence, string accession, Dictionary<int, List<MetaMorpheusModification>> oneBasedPossibleLocalizedModifications, int[] beginPositions, int[] endPositions, string[] bigPeptideTypes, string name, string fullName, int offset, bool isDecoy, bool isContaminant)
        {
            BaseSequence = baseSequence;
            Accession = accession;
            OneBasedPossibleLocalizedModifications = oneBasedPossibleLocalizedModifications;
            OneBasedBeginPositions = beginPositions;
            OneBasedEndPositions = endPositions;
            BigPeptideTypes = bigPeptideTypes;
            Name = name;
            FullName = fullName;
            Offset = offset;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
        }

        #endregion Public Constructors

        #region Public Properties

        public int[] OneBasedBeginPositions { get; private set; }
        public int[] OneBasedEndPositions { get; private set; }
        public string[] BigPeptideTypes { get; private set; }
        public Dictionary<int, List<MetaMorpheusModification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public string Accession { get; private set; }
        public string BaseSequence { get; private set; }
        public bool IsDecoy { get; private set; }

        public int Length
        {
            get
            {
                return BaseSequence.Length;
            }
        }

        public string FullDescription
        {
            get
            {
                if (fullDescription == null)
                {
                    fullDescription = Accession + "|" + Name + "|" + FullName;
                }
                return fullDescription;
            }
        }

        public string Name { get; private set; }

        public string FullName { get; private set; }

        public int Offset { get; private set; }
        public bool IsContaminant { get; set; }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                return BaseSequence[zeroBasedIndex];
            }
        }

        #endregion Public Indexers

        #region Public Methods

        public IEnumerable<PeptideWithPossibleModifications> Digest(Protease protease, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior)
        {
            if (protease.CleavageSpecificity != CleavageSpecificity.None)
            {
                // these are the 1-based residue indices the protease cleaves AFTER
                List<int> oneBasedIndicesToCleaveAfter = protease.GetDigestionSiteIndices(BaseSequence);
                // Cleave after 0, or before the first one
                oneBasedIndicesToCleaveAfter.Insert(0, 0);
                // Cleave after Length index
                oneBasedIndicesToCleaveAfter.Add(Length);

                if (protease.CleavageSpecificity == CleavageSpecificity.Full)
                {
                    for (int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
                    {
                        for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missed_cleavages - 1; i++)
                        {
                            // Retain!
                            if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || i != 0 || this[0] != 'M')
                            {
                                yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], this, missed_cleavages, "full");
                            }
                            // Cleave!
                            if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && i == 0 && this[0] == 'M')
                            {
                                yield return new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], this, missed_cleavages, "full:M cleaved");
                            }
                        }

                        // Also digest using the chain peptide start/end indices
                        for (int chainPeptideIndex = 0; chainPeptideIndex < OneBasedBeginPositions.Length; chainPeptideIndex++)
                        {
                            if (OneBasedBeginPositions[chainPeptideIndex] != 1 || OneBasedEndPositions[chainPeptideIndex] != Length)
                            {
                                int i = 0;
                                while (oneBasedIndicesToCleaveAfter[i] < OneBasedBeginPositions[chainPeptideIndex])
                                    i++;
                                // Start peptide
                                if (i + missed_cleavages < oneBasedIndicesToCleaveAfter.Count && oneBasedIndicesToCleaveAfter[i + missed_cleavages] <= OneBasedEndPositions[chainPeptideIndex])
                                    yield return new PeptideWithPossibleModifications(OneBasedBeginPositions[chainPeptideIndex], oneBasedIndicesToCleaveAfter[i + missed_cleavages], this, missed_cleavages, BigPeptideTypes[chainPeptideIndex] + " start");

                                while (oneBasedIndicesToCleaveAfter[i] < OneBasedEndPositions[chainPeptideIndex])
                                    i++;
                                // End
                                if (i - missed_cleavages - 1 >= 0 && oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 >= OneBasedBeginPositions[chainPeptideIndex])
                                    yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, OneBasedEndPositions[chainPeptideIndex], this, missed_cleavages, BigPeptideTypes[chainPeptideIndex] + " end");
                            }
                        }
                    }
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
            else  // protease.CleavageSpecificity == CleavageSpecificity.None
            {
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || this[0] != 'M')
                {
                    yield return new PeptideWithPossibleModifications(1, Length, this, 0, "full");
                }
                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && this[0] == 'M')
                {
                    yield return new PeptideWithPossibleModifications(2, Length, this, 0, "full:M cleaved");
                }

                // Also digest using the chain peptide start/end indices
                for (int chainPeptideIndex = 0; chainPeptideIndex < OneBasedBeginPositions.Length; chainPeptideIndex++)
                {
                    yield return new PeptideWithPossibleModifications(OneBasedBeginPositions[chainPeptideIndex], OneBasedEndPositions[chainPeptideIndex], this, 0, BigPeptideTypes[chainPeptideIndex]);
                }
            }
        }

        public override bool Equals(object obj)
        {
            var q = obj as Protein;
            return q != null && q.Accession.Equals(Accession);
        }

        public override int GetHashCode()
        {
            return Accession.GetHashCode();
        }

        #endregion Public Methods

    }
}