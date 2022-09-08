using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;

namespace EngineLayer.MbrAnalysis
{
    public class SpectralComparison
    {
        public string SpectraSource;
        public string LibrarySpectraSource;
        public string AcceptorSequence;
        public string DonorSequence;
        public double? Homology;
        public Dictionary<string, double?> ScoreDictionary;
        public SpectralSimilarity Similarity;

        public SpectralComparison(string spectraSource, string librarySpectraSource, string fullSequence, 
            SpectralSimilarity similarity, string donorSequence = null, double? homology = null)
        {
            SpectraSource = spectraSource;
            LibrarySpectraSource = librarySpectraSource;
            AcceptorSequence = fullSequence;
            DonorSequence = donorSequence;
            ScoreDictionary = MakeScoreDictionary(similarity);
            Similarity = similarity;
            Homology = homology;
        }

        private Dictionary<string, double?> MakeScoreDictionary(SpectralSimilarity similarity)
        {
            Dictionary<string, double?> scoreDictionary = new Dictionary<string, double?>();
            scoreDictionary.Add("Cosine", similarity.CosineSimilarity());
            scoreDictionary.Add("Spectral Contrast", similarity.SpectralContrastAngle());
            return scoreDictionary;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(SpectraSource);
            sb.Append('\t');
            sb.Append(LibrarySpectraSource);
            sb.Append('\t');
            sb.Append(AcceptorSequence);
            sb.Append('\t');
            if(DonorSequence != null)
            {
                sb.Append(DonorSequence);
                sb.Append('\t');
            }
            if (Homology != null)
            {
                sb.Append(((double)Homology).ToString("F4"));
                sb.Append('\t');
            }
            foreach (KeyValuePair<string, double?> pair in ScoreDictionary)
            {
                // Currently, column names are hard coded in SpectralScoreDistribution
                string score = (pair.Value == null) ? "-1" : ((double)pair.Value).ToString("F4");
                sb.Append(score);
                sb.Append('\t');
            }


            foreach(var x in Similarity.theoreticalXArray)
            {
                sb.Append(x.ToString() + ", ");
            }
            sb.Append('\t');
            foreach (var y in Similarity.experimentalXArray)
            {
                sb.Append(y.ToString() + ", ");
            }
            sb.Append('\t');
            foreach (var x in Similarity.theoreticalXArray)
            {
                foreach (var y in Similarity.experimentalXArray)
                {
                    if (Math.Abs(x - y) < 0.1)
                    {
                        sb.Append(y.ToString() + ", ");
                    }
                        
                }
            }

            return sb.ToString().Trim();
        }
    }
}
