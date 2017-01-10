using System;
using System.Globalization;
using System.IO;

namespace OldInternalLogic
{
    public static class AminoAcidMasses
    {
        private static float[] MONOISOTOPIC_AMINO_ACID_MASSES;

        public static void LoadAminoAcidMasses()
        {
            MONOISOTOPIC_AMINO_ACID_MASSES = new float['Z' - 'A' + 1];
            for (int i = 0; i < MONOISOTOPIC_AMINO_ACID_MASSES.Length; i++)
            {
                MONOISOTOPIC_AMINO_ACID_MASSES[i] = float.NaN;
            }

            using (StreamReader amino_acids = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "amino_acids.tsv")))
            {
                amino_acids.ReadLine();

                while (amino_acids.Peek() != -1)
                {
                    string line = amino_acids.ReadLine();
                    string[] fields = line.Split('\t');
                    char one_letter_code = char.Parse(fields[0]);
                    float monoisotopic_mass = float.Parse(fields[1], CultureInfo.InvariantCulture);
                    MONOISOTOPIC_AMINO_ACID_MASSES[one_letter_code - 'A'] = monoisotopic_mass;
                }
            }
        }

        public static float GetMonoisotopicMass(char aminoAcid)
        {
            return MONOISOTOPIC_AMINO_ACID_MASSES[aminoAcid - 'A'];
        }
    }
}