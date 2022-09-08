using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using System.IO;

namespace EngineLayer.MbrAnalysis
{
    public abstract class SpectralScoreDistribution
    {
        public List<SpectralComparison> SpectralComparisonList { get; private set; }
        public string OutputFolder { get; set; }
        public string DistributionType { get; private set; }

        protected SpectralScoreDistribution(string outputFolder)
        {
            OutputFolder = outputFolder; 
            System.IO.Directory.CreateDirectory(outputFolder);
            SpectralComparisonList = new();
        }

        protected void SetDistributionType(string type)
        {
            DistributionType = type;
        }

        public void AddComparison(string spectraSource, string librarySpectraSource, string fullSequence,
            SpectralSimilarity similarity, string donorSequence = null, double? homology = null)
        {
            SpectralComparisonList.Add(new SpectralComparison(spectraSource, librarySpectraSource,
                fullSequence, similarity, donorSequence, homology));
        }

        public void WriteResults()
        {
            string fileName = DistributionType + "Distribution.tsv";
            string outputFile = Path.Combine(OutputFolder, fileName);
            using (var output = new StreamWriter(outputFile))
            {
                output.WriteLine(ScoreDistributionHeader);
                foreach (SpectralComparison sc in SpectralComparisonList)
                {
                    output.WriteLine(sc.ToString());
                }
            }
        }

        public string ScoreDistributionHeader 
        {
            get
            {
                StringBuilder sb = new();
                sb.Append("Spectra File");
                sb.Append('\t');
                sb.Append("Library Spectra File");
                sb.Append('\t');
                sb.Append("Acceptor Sequence");
                sb.Append('\t');
                if(DistributionType.Equals("TrueNegative")) 
                {
                    sb.Append("Donor Sequence");
                    sb.Append('\t');
                    sb.Append("Sequence Homology");
                    sb.Append('\t');
                }
                // Cosine Angle and Spectral Contrast are hard coded and they shouldn't be.
                // This should be extensible any similarity measurement
                sb.Append("Cosine Angle");
                sb.Append('\t');
                sb.Append("Spectral Contrast");
                return sb.ToString();
            }

        }

    }
}
