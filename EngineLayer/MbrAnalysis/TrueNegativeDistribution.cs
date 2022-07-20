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

        public TrueNegativeDistribution(string outputFolder) : base(outputFolder)
        {
            SetDistributionType("TrueNegative");
        }

        public double GetPercentHomology(string acceptorSequence, string donorSequence)
        {
            char[] acceptorArray = ConvertSequence(acceptorSequence);
            char[] donorArray = ConvertSequence(donorSequence);
            if (acceptorArray.Length != donorArray.Length || acceptorArray.Length < 1) return -1;

            int rawScore = 0;
            for (int i = 0; i < acceptorArray.Length; i++)
            {
                if (acceptorArray[i] == donorArray[i]) rawScore++;
            }
            return (double)rawScore / (double)acceptorArray.Length;
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
    }
}
