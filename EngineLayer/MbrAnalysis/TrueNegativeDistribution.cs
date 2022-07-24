using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.MbrAnalysis
{
    public class TrueNegativeDistribution : SpectralScoreDistribution
    {
        private int ModDictChar = 0; //Making this a byte and later converting to char could enable 100+ ptms in the dict

        Dictionary<string, char> ModDict = new();
        Dictionary<char, int> ResidueCount = new();

        private static char[] AminoAcids = new char[] 
        { 
            'A', 'C', 'D', 'E', 'F', 
            'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R',
            'S', 'T', 'V', 'W', 'Y'                                             
        };

        public TrueNegativeDistribution(string outputFolder) : base(outputFolder)
        {
            SetDistributionType("TrueNegative");
        }

        public double GetPercentHomology(string acceptorSequence, string donorSequence)
        {
            if (acceptorSequence == null | donorSequence == null) return -1;
            char[] acceptorArray = ConvertSequence(acceptorSequence);
            char[] donorArray = ConvertSequence(donorSequence);
            if (acceptorArray.Length != donorArray.Length || acceptorArray.Length < 1) return -1;

            int rawScore = 0;
            int maxScore = 0;
            int peptideLength = acceptorArray.Length;
            for (int i = 0; i < peptideLength; i++)
            {
                Dictionary<char, int> acceptorResidueCount = new();
                Dictionary<char, int> donorResidueCount = new();
                // b ion check
                for (int j = 0; j <= i; j++)
                {
                    if (!acceptorResidueCount.ContainsKey(acceptorArray[j])) 
                    {
                        acceptorResidueCount.Add(acceptorArray[j], 1);
                    } else
                    {
                        acceptorResidueCount[acceptorArray[j]]++;
                    }
                    if (!donorResidueCount.ContainsKey(donorArray[j]))
                    {
                        donorResidueCount.Add(donorArray[j], 1);
                    } else
                    {
                        donorResidueCount[donorArray[j]]++;
                    }
                }
                maxScore++;
                if (acceptorResidueCount.Count == donorResidueCount.Count)
                {
                    bool match = true;
                    foreach (char key in acceptorResidueCount.Keys)
                    {
                        if (!donorResidueCount.ContainsKey(key))
                        {
                            match = false;
                            break;
                        }
                        else if (acceptorResidueCount[key] != donorResidueCount[key])
                        {
                            match = false;
                            break;
                        }
                    }
                    if (match) rawScore++;
                }

                //y ion check
                acceptorResidueCount = new();
                donorResidueCount = new();
                int yPosition = peptideLength - (i + 1);
                for (int j = peptideLength-1; j >= yPosition; j--)
                {
                    if (!acceptorResidueCount.ContainsKey(acceptorArray[j]))
                    {
                        acceptorResidueCount.Add(acceptorArray[j], 1);
                    }
                    else
                    {
                        acceptorResidueCount[acceptorArray[j]]++;
                    }
                    if (!donorResidueCount.ContainsKey(donorArray[j]))
                    {
                        donorResidueCount.Add(donorArray[j], 1);
                    }
                    else
                    {
                        donorResidueCount[donorArray[j]]++;
                    }
                }
                maxScore++;
                if (acceptorResidueCount.Count == donorResidueCount.Count)
                {
                    bool match = true;
                    foreach (char key in acceptorResidueCount.Keys)
                    {
                        if (!donorResidueCount.ContainsKey(key))
                        {
                            match = false;
                            break;
                        }
                        else if (acceptorResidueCount[key] != donorResidueCount[key])
                        {
                            match = false;
                            break;
                        }
                    }
                    if (match) rawScore++;
                }
            }
            return (double)rawScore / (double)maxScore;
        }

        public char[] ConvertSequence(string peptideSequence)
        {
            char[] peptideArray = peptideSequence.ToCharArray();
            List<char> peptideSansMods = new();
            bool copyOn = true;
            List<int> bracketPositions = new();
            for (int i = 0; i < peptideArray.Length; i++)
            {
                if (peptideArray[i] == '[' | peptideArray[i] == ']')
                {
                    bracketPositions.Add(i);
                    if (copyOn) peptideSansMods.Add('X');
                    copyOn = !copyOn;
                }
                if (copyOn)
                {
                    peptideSansMods.Add(peptideArray[i]);
                }
            }
            if (bracketPositions.Count % 2 != 0)
            {
                return new char[] { }; // uneven number of brackets indicates unknown formatting was encountered
            }

            for (int i = 0; i < bracketPositions.Count ; i = i + 2)
            {
                string mod = peptideSequence.Substring(bracketPositions[i], bracketPositions[i+1] - bracketPositions[i] + 1);
                if (!ModDict.ContainsKey(mod))
                {
                    if (ModDictChar > 9) return new char[] { };
                    ModDict.Add(mod, ModDictChar.ToString().ToCharArray()[0]);
                    ModDictChar++;
                }
                int y = peptideSansMods.IndexOf('X');
                peptideSansMods[y - 1] = ModDict[mod];
                //Remove opening and closing brackets (Where opening bracket was replaced with 'X')
                peptideSansMods.RemoveAt(y);
                peptideSansMods.RemoveAt(peptideSansMods.IndexOf(']')); 
            }
            return peptideSansMods.ToArray();
        }

        public static string ScoreDistributionHeader
        {
            get
            {
                StringBuilder sb = new();
                sb.Append("Spectra File");
                sb.Append('\t');
                sb.Append("Library Spectra File");
                sb.Append('\t');
                sb.Append("Full Sequence");
                sb.Append('\t');
                sb.Append("Donor Sequence");
                sb.Append('\t');
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
