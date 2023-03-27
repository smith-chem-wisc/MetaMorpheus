using System.Collections.Generic;
using System.IO;
using MassSpectrometry;
using System.Globalization;
using Proteomics;

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
                    case "":
                        break;
                    case "CID":
                        cleaveDissociationTypes.Add(DissociationType.CID);
                        break;
                    case "HCD":
                        cleaveDissociationTypes.Add(DissociationType.HCD);
                        break;
                    case "ETD":
                        cleaveDissociationTypes.Add(DissociationType.ETD);
                        break;
                    case "ETHCD":
                        cleaveDissociationTypes.Add(DissociationType.EThcD);
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
                totalMass: double.Parse(split[5], CultureInfo.InvariantCulture),
                cleaveMassShort: double.Parse(split[6], CultureInfo.InvariantCulture),
                cleaveMassLong: double.Parse(split[7], CultureInfo.InvariantCulture),
                loopMass: double.Parse(split[5], CultureInfo.InvariantCulture),
                deadendMassH2O: double.Parse(split[8], CultureInfo.InvariantCulture),
                deadendMassNH2: double.Parse(split[9], CultureInfo.InvariantCulture),
                deadendMassTris: double.Parse(split[10], CultureInfo.InvariantCulture));

            return crosslinker;
        }

        /// <summary>
        /// Deadend and loop crosslink are changed to a type of modification and searched as a modified peptide.
        /// </summary>
        /// <param name="crosslinker"></param>
        public static void GenerateCrosslinkModifications(Crosslinker crosslinker, out Modification TrisDeadEnd, out Modification H2ODeadEnd, out Modification NH2DeadEnd, out Modification Loop)
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            TrisDeadEnd = new Modification(_originalId: "Tris Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: crosslinker.DeadendMassTris);
            H2ODeadEnd = new Modification(_originalId: "H2O Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: crosslinker.DeadendMassH2O);
            NH2DeadEnd = new Modification(_originalId: "NH2 Dead End", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: crosslinker.DeadendMassNH2);
            Loop = new Modification(_originalId: "Loop", _modificationType: "Crosslink", _locationRestriction: "Anywhere.", _target: motif, _monoisotopicMass: crosslinker.LoopMass);
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
            return x;
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