using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;

namespace EngineLayer
{

    public class GlycanDatabase
    {
        //Load Glycan. Generally, glycan-ions should be generated for N-Glycopepitdes which produce Y-ions; MS method couldn't produce o-glycan-ions.
        public static IEnumerable<Glycan> LoadGlycan(string filePath, bool toGenerateIons = true)
        {
            bool isKind = true;
            using (StreamReader lines = new StreamReader(filePath))
            {
                while(lines.Peek() != -1)
                {
                    string line = lines.ReadLine();
                    if (!line.Contains("HexNAc"))
                    {
                        isKind = false;
                    }
                    break;
                }
            }

            if (isKind)
            {
                return LoadKindGlycan(filePath, toGenerateIons);
            }
            else
            {
                return LoadStructureGlycan(filePath);
            }
        }

        //Load KindGlycan. Compatible with Byonic.
        public static IEnumerable<Glycan> LoadKindGlycan(string filePath, bool toGenerateIons)
        {
            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine().Split('\t').First();

                    if (!line.Contains("HexNAc"))
                    {
                        continue;
                    }

                    byte[] kind = new byte[
                        ] { 0, 0, 0, 0, 0, 0, 0, 0 };
                    var x = line.Split('(', ')');
                    int i = 0;
                    while (i < x.Length - 1)
                    {
                        kind[Glycan.NameCharDic[x[i]].Item2] = byte.Parse(x[i + 1]);
                        i = i + 2;
                    }
                    var mass = Glycan.GetMass(kind);

                    if (toGenerateIons)
                    {

                    }
                    else
                    {
                        var glycan = new Glycan(kind);
                        glycan.GlyId = id++;
                        yield return glycan;
                    }

                }
            }
        }

        //Load structured Glycan database.
        public static IEnumerable<Glycan> LoadStructureGlycan(string filePath)
        {
            using (StreamReader glycans = new StreamReader(filePath))
            {
                int id = 1;
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();
                    yield return Glycan.Struct2Glycan(line, id++);
                }
            }
        }

        #region LoadKindGlycan based on Structured Glycan

        public static IEnumerable<Glycan> LoadKindGlycan(string filePath, IEnumerable<Glycan> NGlycans)
        {
            var groupedGlycans = NGlycans.GroupBy(p => Glycan.GetKindString(p.Kind)).ToDictionary(p => p.Key, p => p.ToList());

            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine().Split('\t').First();

                    byte[] kind = new byte[8] { 0, 0, 0, 0, 0, 0, 0, 0 };
                    var x = line.Split('(', ')');
                    int i = 0;
                    while (i < x.Length - 1)
                    {
                        kind[Glycan.NameCharDic[x[i]].Item2] = byte.Parse(x[i + 1]);
                        i = i + 2;
                    }

                    var mass = Glycan.GetMass(kind);

                    var glycans = GetAllIonMassFromKind(kind, groupedGlycans);

                    var glycan = new Glycan(glycans.First().Struc, mass, kind, glycans.First().Ions, true);
                    glycan.GlyId = id++;
                    yield return glycan;
                }
            }
        }

        //Find glycans from structured glycan database
        public static List<Glycan> GetAllIonMassFromKind(byte[] kind, Dictionary<string, List<Glycan>> groupedGlycans)
        {
            var kindKey = Glycan.GetKindString(kind);
            List<Glycan> glycans = new List<Glycan>();

            groupedGlycans.TryGetValue(kindKey, out glycans);

            if (glycans == null)
            {
                //if not in the structured glycan database, find a smaller one.
                bool notFound = true;
                while (notFound)
                {
                    var childKinds = BuildChildKindKey(kind);
                    foreach (var child in childKinds)
                    {
                        var key = Glycan.GetKindString(child);
                        if (groupedGlycans.TryGetValue(key, out glycans))
                        {
                            notFound = false;
                            break;
                        }
                    }

                    if (notFound == true)
                    {
                        glycans = GetAllIonMassFromKind(childKinds[0], groupedGlycans);
                        notFound = false;
                    }
                }
            }

            return glycans;
        }

        private static List<byte[]> BuildChildKindKey(byte[] kind)
        {
            List<byte[]> childKinds = new List<byte[]>();
            for (int i = kind.Length - 1; i >= 0; i--)
            {
                if (kind[i] >= 1)
                {
                    var childKind = new byte[kind.Length];
                    Array.Copy(kind, childKind, kind.Length);
                    childKind[i]--;
                    childKinds.Add(childKind);
                }
            }
            return childKinds;
        }

        #endregion
    }
}
