using System;
using System.Collections.Generic;

namespace OldInternalLogic
{
    public class Protein
    {
        #region Private Fields

        private string fullDescription;

        #endregion Private Fields

        #region Public Constructors

        public Protein(string baseSequence, string accession, string dataset_abbreviation, Dictionary<int, List<MorpheusModification>> oneBasedPossibleLocalizedModifications, int[] beginPositions, int[] endPositions, string[] bigPeptideTypes, string name, string fullName, int offset, bool isDecoy)
        {
            this.BaseSequence = baseSequence;
            this.Accession = accession;
            this.dataset_abbreviation = dataset_abbreviation;
            this.OneBasedPossibleLocalizedModifications = oneBasedPossibleLocalizedModifications;
            this.oneBasedBeginPositions = beginPositions;
            this.oneBasedEndPositions = endPositions;
            this.bigPeptideTypes = bigPeptideTypes;
            this.name = name;
            this.fullName = fullName;
            this.offset = offset;
            this.isDecoy = isDecoy;
        }

        #endregion Public Constructors

        #region Public Properties

        public int[] oneBasedBeginPositions { get; private set; }
        public int[] oneBasedEndPositions { get; private set; }
        public string[] bigPeptideTypes { get; private set; }
        public Dictionary<int, List<MorpheusModification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public string Accession { get; private set; }
        public string BaseSequence { get; private set; }
        public string dataset_abbreviation { get; private set; }
        public bool isDecoy { get; private set; }

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
                    fullDescription = Accession + "|" + name + "|" + fullName;
                }
                return fullDescription;
            }
        }

        public string name { get; private set; }

        public string fullName { get; private set; }

        public int offset { get; private set; }

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
            //p.po.RTBoutput("Digesting " + this.BaseSequence);
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
                        //p.po.RTBoutput("missed_cleavages: " + missed_cleavages);
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
                                //p.po.RTBoutput("Start of REG cleave!:");
                                //p.po.RTBoutput(" Starting index: " + (2));
                                //p.po.RTBoutput("Ending index: " + (oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1]));
                                yield return new PeptideWithPossibleModifications(2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1], this, missed_cleavages, "full:M cleaved");
                            }
                        }

                        // Also digest using the chain peptide start/end indices

                        for (int chainPeptideIndex = 0; chainPeptideIndex < oneBasedBeginPositions.Length; chainPeptideIndex++)
                        {
                            if (oneBasedBeginPositions[chainPeptideIndex] != 1 || oneBasedEndPositions[chainPeptideIndex] != Length)
                            {
                                int i = 0;
                                while (oneBasedIndicesToCleaveAfter[i] < oneBasedBeginPositions[chainPeptideIndex])
                                    i++;
                                // Start peptide
                                //p.po.RTBoutput("Start of chain:");
                                //p.po.RTBoutput(" Starting index: " + protein.oneBasedBeginPositions[chainPeptideIndex]);
                                //p.po.RTBoutput(" Ending index: " + oneBasedIndicesToCleaveAfter[i + missed_cleavages ]);
                                if (i + missed_cleavages < oneBasedIndicesToCleaveAfter.Count && oneBasedIndicesToCleaveAfter[i + missed_cleavages] <= oneBasedEndPositions[chainPeptideIndex])
                                    yield return new PeptideWithPossibleModifications(oneBasedBeginPositions[chainPeptideIndex], oneBasedIndicesToCleaveAfter[i + missed_cleavages], this, missed_cleavages, this.bigPeptideTypes[chainPeptideIndex] + " start");

                                while (oneBasedIndicesToCleaveAfter[i] < oneBasedEndPositions[chainPeptideIndex])
                                    i++;
                                // End
                                //p.po.RTBoutput("End of chain:");
                                //p.po.RTBoutput(" Starting index: " + (oneBasedIndicesToCleaveAfter[i - missed_cleavages-1]+1));
                                //p.po.RTBoutput(" Ending index: " + protein.oneBasedEndPositions[chainPeptideIndex]);
                                if (i - missed_cleavages - 1 >= 0 && oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 >= oneBasedBeginPositions[chainPeptideIndex])
                                    yield return new PeptideWithPossibleModifications(oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, oneBasedEndPositions[chainPeptideIndex], this, missed_cleavages, this.bigPeptideTypes[chainPeptideIndex] + " end");
                            }
                        }
                    }
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
            else
            {
                throw new NotImplementedException();
            }
            //else  // protease.CleavageSpecificity == CleavageSpecificity.Semi || protease.CleavageSpecificity == CleavageSpecificity.SemiN || protease.CleavageSpecificity == CleavageSpecificity.SemiC
            //{
            //    if (protease.CleavageSpecificity == CleavageSpecificity.Semi || protease.CleavageSpecificity == CleavageSpecificity.SemiN)
            //    {
            //        for (int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
            //        {
            //            for (int i = 0; i < indices.Count - missed_cleavages - 1; i++)
            //            {
            //                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || indices[i] + 1 != 0 || this[0] != 'M')
            //                {
            //                    // conditional ensures that we are generating peptides at their lowest missed cleavage state
            //                    for (int length = indices[i + missed_cleavages + 1] - indices[i]; length > (indices[i + missed_cleavages + 1] - indices[i]) - (indices[i + missed_cleavages + 1] - indices[(i + missed_cleavages + 1) - 1]); length--)
            //                    {
            //                        if ((indices[i] + 1 + 1) + length - 1 <= Length)
            //                        {
            //                            yield return new PeptideWithPossibleModifications(indices[i] + this.OneBasedStartResidueInProtein, indices[i + missed_cleavages + 1] + this.OneBasedStartResidueInProtein, protein, PeptideDescription);
            //                        }
            //                    }
            //                }

            //                if (initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && indices[i] + 1 == 0 && this[0] == 'M')
            //                {
            //                    // conditional ensures that we are generating peptides at their lowest missed cleavage state
            //                    for (int length = indices[i + missed_cleavages + 1] - indices[i]; length > (indices[i + missed_cleavages + 1] - indices[i]) - (indices[i + missed_cleavages + 1] - indices[(i + missed_cleavages + 1) - 1]); length--)
            //                    {
            //                        if ((indices[i] + 1 + 1 + 1) + length - 1 <= Length)
            //                        {
            //                            yield return new PeptideWithPossibleModifications(indices[i] + this.OneBasedStartResidueInProtein, indices[i + missed_cleavages + 1] + this.OneBasedStartResidueInProtein, protein, PeptideDescription);
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }
            //    if (protease.CleavageSpecificity == CleavageSpecificity.Semi || protease.CleavageSpecificity == CleavageSpecificity.SemiC)
            //    {
            //        for (int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
            //        {
            //            for (int i = 0; i < indices.Count - missed_cleavages - 1; i++)
            //            {
            //                // handling for initiator methionine not required

            //                // - (protease.CleavageSpecificity == CleavageSpecificity.Semi ? 1 : 0) ensures that we don't repeat the same peptides we generated above in the SemiN digestion
            //                // conditional ensures that we are generating peptides at their lowest missed cleavage state
            //                for (int length = indices[i + missed_cleavages + 1] - indices[i] - (protease.CleavageSpecificity == CleavageSpecificity.Semi ? 1 : 0); length > (indices[i + missed_cleavages + 1] - indices[i]) - (indices[i + 1] - indices[i]); length--)
            //                {
            //                    if ((indices[i + missed_cleavages + 1] + 1) - length + 1 >= 1)
            //                    {
            //                        yield return new PeptideWithPossibleModifications(indices[i] + this.OneBasedStartResidueInProtein, indices[i + missed_cleavages + 1] + this.OneBasedStartResidueInProtein, protein, PeptideDescription);
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}
            //}
        }

        public override bool Equals(object obj)
        {
            Protein q = obj as Protein;
            return q != null && q.Accession.Equals(Accession);
        }

        public override int GetHashCode()
        {
            return Accession.GetHashCode();
        }

        #endregion Public Methods
    }
}