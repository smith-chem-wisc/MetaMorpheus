using System.Collections.Generic;
using System.IO;
using MassSpectrometry;

namespace EngineLayer
{
    public class Crosslinker
    {
        /// <summary>
        /// Do not use this constructor for anything. It exists so that the .toml can be read.
        /// </summary>
        public Crosslinker()
        {

        }

        public Crosslinker(string crosslinkerModSites, string crosslinkerModSites2, string crosslinkerName, bool cleavable, string dissociationTypes, double totalMass,
            double cleaveMassShort, double cleaveMassLong, double loopMass, double deadendMassH2O, double deadendMassNH2, double deadendMassTris)
        {
            CrosslinkerModSites = crosslinkerModSites;
            CrosslinkerModSites2 = crosslinkerModSites2;
            CrosslinkerName = crosslinkerName;
            Cleavable = cleavable;
            CleaveDissociationTypes = GetCleaveDissociationTypes(dissociationTypes);
            TotalMass = totalMass;
            CleaveMassShort = cleaveMassShort;
            CleaveMassLong = cleaveMassLong;
            LoopMass = loopMass;
            DeadendMassH2O = deadendMassH2O;
            DeadendMassNH2 = deadendMassNH2;
            DeadendMassTris = deadendMassTris;
        }

        public string CrosslinkerModSites { get; set; }
        public string CrosslinkerModSites2 { get; set; }
        public string CrosslinkerName { get; set; }
        public bool Cleavable { get; set; }
        public List<DissociationType> CleaveDissociationTypes { get; set; }
        public double TotalMass { get; set; }
        public double CleaveMassShort { get; set; }
        public double CleaveMassLong { get; set; }
        public double LoopMass { get; set; }
        public double DeadendMassH2O { get; set; }
        public double DeadendMassNH2 { get; set; }
        public double DeadendMassTris { get; set; }

        private List<DissociationType> GetCleaveDissociationTypes(string cleaveDissociationTypesInString)
        {
            List<DissociationType> cleaveDissociationTypes = new List<DissociationType>();
            foreach (var x in cleaveDissociationTypesInString.Split('|'))
            {
                switch (x)
                {
                    case "CID":
                        cleaveDissociationTypes.Add(DissociationType.CID);
                        break;
                    case "HCD":
                        cleaveDissociationTypes.Add(DissociationType.HCD);
                        break;
                    case "ETD":
                        cleaveDissociationTypes.Add(DissociationType.ETD);
                        break;
                    default:
                        break;
                }
            }
            return cleaveDissociationTypes;
        }

        public static IEnumerable<Crosslinker> LoadCrosslinkers(string CrosslinkerLocation)
        {
            using (StreamReader crosslinkers = new StreamReader(CrosslinkerLocation))
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

                    yield return ParseCrosslinkerFromString(line);
                }
            }
        }

        public static Crosslinker ParseCrosslinkerFromString(string line)
        {
            var split = line.Split('\t');
            bool cleavable = true;
            if (split[3] == "F")
            {
                cleavable = false;
            }
  
            Crosslinker crosslinker = new Crosslinker(
                crosslinkerName: split[0],
                crosslinkerModSites: split[1],
                crosslinkerModSites2: split[2],
                cleavable: cleavable,
                dissociationTypes: split[4],
                totalMass: double.Parse(split[5]),
                cleaveMassShort: double.Parse(split[6]),
                cleaveMassLong: double.Parse(split[7]),
                loopMass: double.Parse(split[5]),
                deadendMassH2O: double.Parse(split[8]),
                deadendMassNH2: double.Parse(split[9]),
                deadendMassTris: double.Parse(split[10]));

            return crosslinker;
        }

        public override string ToString()
        {
            return CrosslinkerName;
        }

        public static string DissociationTypes2String(List<DissociationType> dissociationTypes)
        {
            string x = "";
            foreach (var d in dissociationTypes)
            {
                x += d.ToString() + "|";
            }          
            return x.Remove(x.Length - 1);
        }

        public string ToString(bool writeCrosslinker)
        {
            if (writeCrosslinker)
            {               
                return (CrosslinkerName + "\t" + CrosslinkerModSites + "\t" + CrosslinkerModSites2 + "\t" + Cleavable + "\t" + DissociationTypes2String(CleaveDissociationTypes) + "\t" + TotalMass + "\t" + CleaveMassShort + "\t"
                    + CleaveMassLong + "\t" + DeadendMassH2O + "\t" + DeadendMassNH2 + "\t" + DeadendMassTris);
            }
            else
            {
                return CrosslinkerName;
            }
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