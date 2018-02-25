using System.Text;

namespace EngineLayer.Neo
{
    public class PsmTsvLine
    {
        #region Public Constructors

        public PsmTsvLine(string[] line, int scanNumber, double score, string baseSequence, string fullSequence, string accession, string proteinName, string geneName, string DCT, string target, string decoy, string q)
        {
            this.line = line;
            this.scanNumber = scanNumber;
            this.score = score;
            this.baseSequence = baseSequence;
            this.fullSequence = fullSequence;
            this.accession = accession;
            this.proteinName = proteinName;
            this.geneName = geneName;
            this.DCT = DCT;
            this.target = target;
            this.decoy = decoy;
            this.q = q;
            neoType = null;
        }

        #endregion Public Constructors

        #region Public Enums

        public enum NeoType { Normal, Spliced, DecoySpliced };

        #endregion Public Enums

        #region Public Properties

        public string[] line { get; set; }
        public int scanNumber { get; set; }
        public double score { get; set; }
        public string baseSequence { get; set; }
        public string fullSequence { get; set; }
        public string accession { get; set; }
        public string proteinName { get; set; }
        public string geneName { get; set; }
        public string DCT { get; set; }
        public string target { get; set; }
        public string decoy { get; set; }
        public string q { get; set; }
        public NeoType? neoType { get; set; }

        #endregion Public Properties

        #region Public Methods

        public PsmTsvLine AggregateLine(PsmTsvLine secondary)
        {
            baseSequence += "|" + secondary.baseSequence;
            fullSequence += "|" + secondary.fullSequence;
            accession += "|" + secondary.accession;
            proteinName += "|" + secondary.proteinName;
            geneName += "|" + secondary.geneName;
            return this;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < line.Length; i++)
            {
                if (i == ImportPsmtsv.baseIndex)
                    sb.Append(baseSequence + '\t');
                else if (i == ImportPsmtsv.fullIndex)
                    sb.Append(fullSequence + '\t');
                else if (i == ImportPsmtsv.accessionIndex)
                    sb.Append(accession + '\t');
                else if (i == ImportPsmtsv.proteinIndex)
                    sb.Append(proteinName + '\t');
                else if (i == ImportPsmtsv.geneIndex)
                    sb.Append(geneName + '\t');
                else if (i == ImportPsmtsv.targetIndex)
                    sb.Append(target + '\t');
                else if (i == ImportPsmtsv.decoyIndex)
                    sb.Append(decoy + '\t');
                else if (i == ImportPsmtsv.qIndex)
                    sb.Append(q + '\t');
                else
                    sb.Append(line[i] + '\t');
            }
            if (neoType != null)
                sb.Append(neoType.ToString());
            return sb.ToString();
        }

        #endregion Public Methods
    }
}