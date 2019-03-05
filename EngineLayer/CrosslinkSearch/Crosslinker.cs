using System.IO;
using System.Collections.Generic;

namespace EngineLayer
{
    public class Crosslinker
    {
        public Crosslinker(string crosslinkerModSites = "K", string crosslinkerModSites2 ="K", string crosslinkerName = "DSSO", bool cleavable = true, double totalMass = 158.0038,
            double cleaveMassShort = 54.01056, double cleaveMassLong = 85.982635, double loopMass = 158.0038, double deadendMassH2O = 176.0143, double deadendMassNH2 = 175.0303, double deadendMassTris = 279.0777)
        {
            CrosslinkerModSites = crosslinkerModSites;
            CrosslinkerModSites2 = crosslinkerModSites2;
            CrosslinkerName = crosslinkerName;
            Cleavable = cleavable;
            TotalMass = totalMass;
            CleaveMassShort = cleaveMassShort;
            CleaveMassLong = cleaveMassLong;
            LoopMass = loopMass;
            DeadendMassH2O = deadendMassH2O;
            DeadendMassNH2 = deadendMassNH2;
            DeadendMassTris = deadendMassTris;
        }

        public Crosslinker()
        {
        }

        public string CrosslinkerModSites { get; set; }
        public string CrosslinkerModSites2 { get; set; }
        public string CrosslinkerName { get; set; }
        public bool Cleavable { get; set; }
        public double TotalMass { get; set; }
        public double CleaveMassShort { get; set; }
        public double CleaveMassLong { get; set; }
        public double LoopMass { get; set; }
        public double DeadendMassH2O { get; set; }
        public double DeadendMassNH2 { get; set; }
        public double DeadendMassTris { get; set; }

        public static IEnumerable<Crosslinker> LoadCrosslinkers(string CrosslinkerLocation)
        {
            using(StreamReader crosslinkers = new StreamReader(CrosslinkerLocation))
            {
                int lineCount = 0;
                while (crosslinkers.Peek() != -1)
                {
                    lineCount++;
                    string line = crosslinkers.ReadLine();
                    if (lineCount == 1)
                    {
                        continue;
                    }

                    yield return GetCrosslinker(line);
                }
            }
        }

        public static Crosslinker GetCrosslinker(string line)
        {
            var split = line.Split('\t');
            bool cleable = true;
            if (split[3] == "F")
            {
                cleable = false;
            }

            Crosslinker crosslinker = new Crosslinker(crosslinkerName: split[0], crosslinkerModSites : split[1], crosslinkerModSites2: split[2], 
                cleavable: cleable, totalMass: double.Parse(split[4]), cleaveMassShort: double.Parse(split[5]), cleaveMassLong: double.Parse(split[6]), 
                loopMass: double.Parse(split[4]),deadendMassH2O: double.Parse(split[7]), deadendMassNH2: double.Parse(split[8]), 
                deadendMassTris: double.Parse(split[9]));

            return crosslinker;
        }

        public override string ToString()
        {
            return CrosslinkerName;
        }

        public override bool Equals(object obj)
        {
            var a = obj as Crosslinker;
            return a != null
                && (a.CrosslinkerName == null && CrosslinkerName == null || a.CrosslinkerName.Equals(CrosslinkerName));
        }

        public override int GetHashCode()
        {
            return (CrosslinkerName ?? "").GetHashCode();
        }
    }
}