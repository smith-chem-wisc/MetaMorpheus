using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class Protease
    {
        #region Public Constructors

        public Protease(string name, IEnumerable<string> sequencesInducingCleavage, IEnumerable<string> sequencesPreventingCleavage, TerminusType cleavageTerminus, CleavageSpecificity cleavageSpecificity, string psiMSAccessionNumber, string psiMSName, string siteRegexp)
        {
            Name = name;
            SequencesInducingCleavage = sequencesInducingCleavage;
            SequencesPreventingCleavage = sequencesPreventingCleavage;
            CleavageTerminus = cleavageTerminus;
            CleavageSpecificity = cleavageSpecificity;
            PsiMsAccessionNumber = psiMSAccessionNumber;
            PsiMsName = psiMSName;
            SiteRegexp = siteRegexp;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Name { get; }
        public TerminusType CleavageTerminus { get; }
        public IEnumerable<string> SequencesInducingCleavage { get; }
        public IEnumerable<string> SequencesPreventingCleavage { get; }
        public CleavageSpecificity CleavageSpecificity { get; }
        public string PsiMsAccessionNumber { get; }
        public string PsiMsName { get; }
        public string SiteRegexp { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object obj)
        {
            var a = obj as Protease;
            return a != null
                && a.Name.Equals(Name);
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        #endregion Public Methods

        #region Digestion Methods

        /// <summary>
        /// Gets the indices after which this protease will cleave a given protein sequence
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <returns></returns>
        internal List<int> GetDigestionSiteIndices(string proteinSequence)
        {
            var indices = new List<int>();

            for (int i = 0; i < proteinSequence.Length - 1; i++)
            {
                foreach (string c in SequencesInducingCleavage)
                {
                    if (SequenceInducesCleavage(proteinSequence, i, c)
                        && !SequencesPreventingCleavage.Any(nc => SequencePreventsCleavage(proteinSequence, i, nc)))
                    {
                        indices.Add(i + 1);
                    }
                }
            }

            indices.Insert(0, 0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide
            indices.Add(proteinSequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide

            return indices;
        }

        /// <summary>
        /// Checks if select subsequence of protein matches sequence that induces cleavage
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <param name="proteinSequenceIndex"></param>
        /// <param name="sequenceInducingCleavage"></param>
        /// <returns></returns>
        private bool SequenceInducesCleavage(string proteinSequence, int proteinSequenceIndex, string sequenceInducingCleavage)
        {
            return (CleavageTerminus != TerminusType.N 
                    && proteinSequenceIndex - sequenceInducingCleavage.Length + 1 >= 0 
                    && proteinSequence.Substring(proteinSequenceIndex - sequenceInducingCleavage.Length + 1, sequenceInducingCleavage.Length).Equals(sequenceInducingCleavage, StringComparison.OrdinalIgnoreCase))
                || (CleavageTerminus == TerminusType.N && proteinSequenceIndex + 1 + sequenceInducingCleavage.Length <= proteinSequence.Length 
                    && proteinSequence.Substring(proteinSequenceIndex + 1, sequenceInducingCleavage.Length).Equals(sequenceInducingCleavage, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Checks if select subsequence of protein matches sequence preventing cleavage
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <param name="proteinSequenceIndex"></param>
        /// <param name="sequencePreventingCleavage"></param>
        /// <returns></returns>
        private bool SequencePreventsCleavage(string proteinSequence, int proteinSequenceIndex, string sequencePreventingCleavage)
        {
            return (CleavageTerminus != TerminusType.N 
                    && proteinSequenceIndex + 1 + sequencePreventingCleavage.Length <= proteinSequence.Length 
                    && proteinSequence.Substring(proteinSequenceIndex + 1, sequencePreventingCleavage.Length).Equals(sequencePreventingCleavage, StringComparison.OrdinalIgnoreCase))
                || (CleavageTerminus == TerminusType.N 
                    && proteinSequenceIndex - sequencePreventingCleavage.Length + 1 >= 0 
                    && proteinSequence.Substring(proteinSequenceIndex - sequencePreventingCleavage.Length + 1, sequencePreventingCleavage.Length).Equals(sequencePreventingCleavage, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Gets intervals of a protein sequence that will result from digestion by this protease.
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        internal List<Peptide> GetDigestionIntervals(Protein protein, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int minPeptidesLength, int maxPeptidesLength)
        {
            List<Peptide> intervals = new List<Peptide>();

            // proteolytic cleavage in one spot
            if (CleavageSpecificity == CleavageSpecificity.SingleN || CleavageSpecificity == CleavageSpecificity.SingleC)
            {
                bool maxTooBig = protein.Length + maxPeptidesLength < 0; //when maxPeptidesLength is too large, it becomes negative and causes issues
                for (int proteinStart = 1; proteinStart <= protein.Length; proteinStart++)
                {
                    if (CleavageSpecificity == CleavageSpecificity.SingleN && OkayMinLength(protein.Length - proteinStart + 1, minPeptidesLength))
                    {
                        //need Math.Max if max length is int.MaxLength, since +proteinStart will make it negative
                        intervals.Add(new Peptide(protein, proteinStart, maxTooBig ? protein.Length : Math.Min(protein.Length, proteinStart + maxPeptidesLength), 0, "SingleN"));
                    }

                    if (CleavageSpecificity == CleavageSpecificity.SingleC && OkayMinLength(proteinStart, minPeptidesLength))
                    {
                        intervals.Add(new Peptide(protein, Math.Max(1, proteinStart - maxPeptidesLength), proteinStart, 0, "SingleC"));
                    }
                }
            }
            else if (CleavageSpecificity == CleavageSpecificity.None)
            {
                // retain methionine
                if ((initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M') && OkayLength(protein.Length, minPeptidesLength, maxPeptidesLength))
                {
                    intervals.Add(new Peptide(protein, 1, protein.Length, 0, "full"));
                }

                // cleave methionine
                if ((initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M') && OkayLength(protein.Length - 1, minPeptidesLength, maxPeptidesLength))
                {
                    intervals.Add(new Peptide(protein, 2, protein.Length, 0, "full:M cleaved"));
                }

                // Also digest using the proteolysis product start/end indices
                intervals.AddRange(
                    protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue
                        && OkayLength(proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1, minPeptidesLength, maxPeptidesLength))
                    .Select(proteolysisProduct =>
                        new Peptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, 0, proteolysisProduct.Type)));
            }

            // Full proteolytic cleavage
            else if (CleavageSpecificity == CleavageSpecificity.Full)
            {
                intervals.AddRange(FullDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptidesLength, maxPeptidesLength));
            }

            // Cleavage rules for nonspecific search
            else if (CleavageSpecificity == CleavageSpecificity.Semi)
            {
                intervals.AddRange(SemiProteolyticDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptidesLength, maxPeptidesLength));
            }
            else
            {
                throw new NotImplementedException();
            }

            return intervals;
        }

        /// <summary>
        /// Gets protein intervals for digestion by this specific protease.
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private IEnumerable<Peptide> FullDigestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptidesLength, int maxPeptidesLength)
        {
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(protein.BaseSequence);
            for (int missed_cleavages = 0; missed_cleavages <= maximumMissedCleavages; missed_cleavages++)
            {
                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missed_cleavages - 1; i++)
                {
                    if (Retain(i, initiatorMethionineBehavior, protein[0])
                        && OkayLength(oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - oneBasedIndicesToCleaveAfter[i], minPeptidesLength, maxPeptidesLength))
                    {
                        yield return new Peptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1],
                            missed_cleavages, "full");
                    }
                    if (Cleave(i, initiatorMethionineBehavior, protein[0])
                        && OkayLength(oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1] - 1, minPeptidesLength, maxPeptidesLength))
                    {
                        yield return new Peptide(protein, 2, oneBasedIndicesToCleaveAfter[i + missed_cleavages + 1],
                            missed_cleavages, "full:M cleaved");
                    }
                }

                // Also digest using the proteolysis product start/end indices
                foreach (var proteolysisProduct in protein.ProteolysisProducts)
                {
                    if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                    {
                        int i = 0;
                        while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedBeginPosition)
                        {
                            i++;
                        }

                        bool startPeptide = i + missed_cleavages < oneBasedIndicesToCleaveAfter.Count
                            && oneBasedIndicesToCleaveAfter[i + missed_cleavages] <= proteolysisProduct.OneBasedEndPosition
                            && proteolysisProduct.OneBasedBeginPosition.HasValue
                            && OkayLength(oneBasedIndicesToCleaveAfter[i + missed_cleavages] - proteolysisProduct.OneBasedBeginPosition.Value + 1, minPeptidesLength, maxPeptidesLength);
                        if (startPeptide)
                        {
                            yield return new Peptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[i + missed_cleavages],
                                missed_cleavages, proteolysisProduct.Type + " start");
                        }

                        while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedEndPosition)
                        {
                            i++;
                        }

                        bool end = i - missed_cleavages - 1 >= 0
                            && oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 >= proteolysisProduct.OneBasedBeginPosition
                            && proteolysisProduct.OneBasedEndPosition.HasValue
                            && OkayLength(proteolysisProduct.OneBasedEndPosition.Value - oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1 - 1, minPeptidesLength, maxPeptidesLength);
                        if (end)
                        {
                            yield return new Peptide(protein, oneBasedIndicesToCleaveAfter[i - missed_cleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value,
                                missed_cleavages, proteolysisProduct.Type + " end");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Gets the protein intervals based on nonspecific digestion rules
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private IEnumerable<Peptide> SemiProteolyticDigestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptidesLength, int maxPeptidesLength)
        {
            List<Peptide> intervals = new List<Peptide>();
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(protein.BaseSequence);
            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1; i++)
            {
                bool retain = Retain(i, initiatorMethionineBehavior, protein[0]);
                bool cleave = Cleave(i, initiatorMethionineBehavior, protein[0]);
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1];

                if (retain)
                {
                    intervals.AddRange(FixedTermini(oneBasedIndicesToCleaveAfter[i], cTerminusProtein, protein, cleave, minPeptidesLength, maxPeptidesLength));
                }

                if (cleave)
                {
                    intervals.AddRange(FixedTermini(1, cTerminusProtein, protein, cleave, minPeptidesLength, maxPeptidesLength));
                }
            }

            //finish C-term of protein caused by loop being "i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1"
            int last = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndexSemi = maximumMissedCleavages < last ? maximumMissedCleavages : last;
            //Fringe C-term peptides
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                //fixedN
                int nTerminusProtein = oneBasedIndicesToCleaveAfter[last - i];
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[last] - 1;
                for (int j = cTerminusProtein; j > nTerminusProtein; j--)
                {
                    if (OkayLength(j - nTerminusProtein, minPeptidesLength, maxPeptidesLength))
                    {
                        intervals.Add(new Peptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, "semiN"));
                    }
                }
            }

            //Fringe N-term peptides
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                //fixedC
                int nTerminusProtein = oneBasedIndicesToCleaveAfter[0];
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[i];
                for (int j = nTerminusProtein + 1; j < cTerminusProtein; j++)
                {
                    if (OkayLength(cTerminusProtein - j, minPeptidesLength, maxPeptidesLength))
                    {
                        intervals.Add(new Peptide(protein, j + 1, cTerminusProtein, cTerminusProtein - j, "semiC"));
                    }
                }
            }

            // Also digest using the proteolysis product start/end indices
            // This should only be things where the proteolysis is not K/R and the
            foreach (var proteolysisProduct in protein.ProteolysisProducts)
            {
                if (proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue
                    && (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)) //if at least one side is not a terminus
                {
                    int i = 0;
                    while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedBeginPosition)//"<" to prevent additions if same index as residues
                    {
                        i++; //last position in protein is an index to cleave after
                    }

                    // Start peptide
                    for (int j = proteolysisProduct.OneBasedBeginPosition.Value; j < oneBasedIndicesToCleaveAfter[i]; j++)
                    {
                        if (OkayLength(j - proteolysisProduct.OneBasedBeginPosition + 1, minPeptidesLength, maxPeptidesLength))
                        {
                            intervals.Add(new Peptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, j,
                                j - proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.Type + " start"));
                        }
                    }
                    while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedEndPosition) //"<" to prevent additions if same index as residues, since i-- is below
                    {
                        i++;
                    }

                    //Now that we've obtained an index to cleave after that is past the proteolysis product
                    //we need to backtrack to get the index to cleave that is immediately before the the proteolysis product
                    //to do this, we will do i--
                    //In the nitch case that the proteolysis product is already an index to cleave
                    //no new peptides will be generated using this, so we will forgo i--
                    //this makes peptides of length 0, which are not generated due to the for loop
                    //removing this if statement will result in crashes from c-terminal proteolysis product end positions
                    if (oneBasedIndicesToCleaveAfter[i] != proteolysisProduct.OneBasedEndPosition)
                    {
                        i--;
                    }

                    // End
                    for (int j = oneBasedIndicesToCleaveAfter[i] + 1; j < proteolysisProduct.OneBasedEndPosition.Value; j++)
                    {
                        if (OkayLength(proteolysisProduct.OneBasedEndPosition - j + 1, minPeptidesLength, maxPeptidesLength))
                        {
                            intervals.Add(new Peptide(protein, j, proteolysisProduct.OneBasedEndPosition.Value,
                                proteolysisProduct.OneBasedEndPosition.Value - j, proteolysisProduct.Type + " end"));
                        }
                    }
                }
            }
            return intervals;
        }

        /// <summary>
        /// Get protein intervals for fixed termini. For semi-proteolytic cleavage.
        /// </summary>
        /// <param name="nTerminusProtein"></param>
        /// <param name="cTerminusProtein"></param>
        /// <param name="protein"></param>
        /// <param name="cleave"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private static IEnumerable<Peptide> FixedTermini(int nTerminusProtein, int cTerminusProtein, Protein protein, bool cleave, int minPeptidesLength, int maxPeptidesLength)
        {
            List<Peptide> intervals = new List<Peptide>();
            if (OkayLength(cTerminusProtein - nTerminusProtein, minPeptidesLength, maxPeptidesLength))
            {
                intervals.Add(new Peptide(protein, nTerminusProtein + 1, cTerminusProtein,
                    cTerminusProtein - nTerminusProtein, "semi" + (cleave ? ":M cleaved" : ""))); // Maximum sequence length
            }

            // Fixed termini at each internal index
            IEnumerable<int> internalIndices = Enumerable.Range(nTerminusProtein + 1, cTerminusProtein - nTerminusProtein - 1);
            IEnumerable<Peptide> fixedCTermIntervals =
                internalIndices
                .Where(j => OkayLength(cTerminusProtein - j, minPeptidesLength, maxPeptidesLength))
                .Select(j => new Peptide(protein, j + 1, cTerminusProtein, cTerminusProtein - j, "semiC" + (cleave ? ":M cleaved" : "")));
            IEnumerable<Peptide> fixedNTermIntervals =
                internalIndices
                .Where(j => OkayLength(j - nTerminusProtein, minPeptidesLength, maxPeptidesLength))
                .Select(j => new Peptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, "semiN" + (cleave ? ":M cleaved" : "")));

            return intervals.Concat(fixedCTermIntervals).Concat(fixedNTermIntervals);
        }

        /// <summary>
        /// Retain N-terminal residue?
        /// </summary>
        /// <param name="oneBasedCleaveAfter"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="nTerminus"></param>
        /// <returns></returns>
        internal static bool Retain(int oneBasedCleaveAfter, InitiatorMethionineBehavior initiatorMethionineBehavior, char nTerminus)
        {
            return oneBasedCleaveAfter != 0 // this only pertains to the n-terminus
                || initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave
                || nTerminus != 'M';
        }

        /// <summary>
        /// Cleave N-terminal residue?
        /// </summary>
        /// <param name="oneBasedCleaveAfter"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="nTerminus"></param>
        /// <returns></returns>
        internal static bool Cleave(int oneBasedCleaveAfter, InitiatorMethionineBehavior initiatorMethionineBehavior, char nTerminus)
        {
            return oneBasedCleaveAfter == 0 // this only pertains to the n-terminus
                && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain
                && nTerminus == 'M';
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum and maximum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        internal static bool OkayLength(int? peptideLength, int? minPeptidesLength, int? maxPeptidesLength)
        {
            return OkayMinLength(peptideLength, minPeptidesLength) && OkayMaxLength(peptideLength, maxPeptidesLength);
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="minPeptidesLength"></param>
        /// <returns></returns>
        internal static bool OkayMinLength(int? peptideLength, int? minPeptidesLength)
        {
            return !minPeptidesLength.HasValue || peptideLength >= minPeptidesLength;
        }

        /// <summary>
        /// Is length of given peptide okay, given maximum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        internal static bool OkayMaxLength(int? peptideLength, int? maxPeptidesLength)
        {
            return !maxPeptidesLength.HasValue || peptideLength <= maxPeptidesLength;
        }

        #endregion Digestion Methods
    }
}