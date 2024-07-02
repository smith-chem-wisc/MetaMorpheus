using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    // in our database, the N-glycan.gdb should be correct to the new format
    // the class for loading glycan database then creeat the glycan object.
    public static class GlycanDatabase 
    {
        
        /// <summary>
        /// Load Glycan from the database file. Generally, glycan-ions should be generated for N-Glycopepitdes which produce Y-ions; MS method couldn't produce o-glycan-ions
        /// </summary>
        /// <param name="filePath"> Database file path</param>
        /// <param name="ToGenerateIons"> Do we need to generate the glycanIon? </param>
        /// <param name="IsOGlycanSearch"></param>
        /// <returns> A glycan object collection </returns>
        public static IEnumerable<Glycan> LoadGlycan(string filePath, bool ToGenerateIons, bool IsOGlycanSearch)
        {
            bool isKind = true;
            using (StreamReader lines = new StreamReader(filePath))
            {
                while(lines.Peek() != -1)
                {
                    string line = lines.ReadLine();
                    if (!line.Contains("HexNAc"))  // use the first line to determine the format (kind / structure) of glycan database.
                    {
                        isKind = false;
                    }
                    break;
                }
            }

            if (isKind)
            {
                return LoadKindGlycan(filePath, ToGenerateIons, IsOGlycanSearch); // open the file of the kind format, example: HexNAc(2)Hex(5)NeuAc(1)Fuc(1)
            }
            else
            {
                return LoadStructureGlycan(filePath, IsOGlycanSearch);            // open the file of the structure format, example: (N(H(A))(A))
            }
        }


        /// <summary>
        /// Load composition format Glycan database, then convert to kind format followed by generating the glycan object.
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="ToGenerateIons"></param>
        /// <param name="IsOGlycanSearch"></param>
        /// <returns>The glycan collection </returns>
        public static IEnumerable<Glycan> LoadKindGlycan(string filePath, bool ToGenerateIons, bool IsOGlycanSearch)
        {
            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                while (lines.Peek() != -1)
                {
                    string line = lines.ReadLine().Split('\t').First();

                    if (!(line.Contains("HexNAc") || line.Contains("Hex"))) // Make sure the line is a glycan line. The line should contain HexNAc or Hex.
                    {
                        continue;
                    }

                    var kind = String2Kind(line);  // Convert the database string to kind[] format (byte array).

                    var glycan = new Glycan(kind); // Use the kind[] to create a glycan object.
                    glycan.GlyId = id++;
                    if (ToGenerateIons)
                    {
                        if (IsOGlycanSearch)
                        {
                            glycan.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        else
                        {
                            glycan.Ions = NGlycanCompositionFragments(kind);
                        }
                    }
                    yield return glycan;
                }
            }
        }

        /// <summary>
        /// Convert the glycan string to Kind array
        /// </summary>
        /// <param name="line"> ex. HexNAc(2)Hex(5)NeuAc(1)Fuc(1) </param>
        /// <returns> The glycan Kind List ex. [2, 5, 0, 0, 1, 0, 0, 0, 0, 1] </returns>
        public static byte[] String2Kind(string line) 
        {
            byte[] kind = new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            var x = line.Split(new char[] { '(', ')' });
            int i = 0;
            while (i < x.Length - 1)
            {
                kind[Glycan.NameCharDic[x[i]].Item2] = byte.Parse(x[i + 1]);
                i = i + 2;
            }

            return kind;
        }

        /// <summary>
        /// Load structured format Glycan database and generate the glycan object.
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="IsOGlycan"></param>
        /// <returns> The Glycan object collection </returns>
        public static IEnumerable<Glycan> LoadStructureGlycan(string filePath, bool IsOGlycan)
        {
            using (StreamReader glycans = new StreamReader(filePath))
            {
                int id = 1;
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();                         // Read the line from the database file. Ex. (N(H(A))(A))
                    yield return Glycan.Struct2Glycan(line, id++, IsOGlycan); // Directly convert the string to Glycan object.
                }
            }
        }

        //This function build fragments based on the general core of NGlyco fragments. 
        //From https://github.com/mobiusklein/glycopeptidepy/structure/fragmentation_strategy/glycan.py#L408
        //The fragment generation is not as good as structure based method. So it is better to use a structure based N-Glycan database.
        // The function is used to load the database from the different formats, but we don't use it now.
        public static List<GlycanIon> NGlycanCompositionFragments(byte[] kind)
        {
            int glycan_mass = Glycan.GetMass(kind);

            int core_count = 1;
            int iteration_count = 0;
            bool extended = true;
            bool extended_fucosylation = false;

            int fuc_count = kind[4];
            int xyl_count = kind[9];
            int hexnac_inaggregate = kind[0];
            int hexose_inaggregate = kind[1];

            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int base_hexnac = Math.Min(hexnac_inaggregate + 1, 3);
            for (int hexnac_count = 0; hexnac_count < base_hexnac; hexnac_count++)
            {
                if (hexnac_count == 0)
                {
                    GlycanIon glycanIon = new GlycanIon(null, 8303819, new byte[] { 0, (byte)hexnac_count, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, glycan_mass - 8303819);
                    glycanIons.Add(glycanIon);
                }
                else if (hexnac_count == 1)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);

                    glycanIons.Add(glycanIon);

                    if (iteration_count < fuc_count)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);

                        glycanIons.Add(fuc_glycanIon);
                    }
                }
                else if (hexnac_count == 2)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);
                    glycanIons.Add(glycanIon);

                    if (!extended_fucosylation)
                    {
                        if (iteration_count < fuc_count)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }
                    }
                    else if (fuc_count > 0)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                        glycanIons.Add(fuc_glycanIon);

                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                        {
                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                            glycanIons.Add(add_fuc_glycanIon);
                        }

                        if (iteration_count < xyl_count)
                        {
                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                            glycanIons.Add(xyl_fuc_glycanIon);
                        }
                    }

                    if (iteration_count < xyl_count)
                    {
                        GlycanIon xyl_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 0, 1, glycan_mass);
                        glycanIons.Add(xyl_glycanIon);
                    }


                    int min_hexose_inaggregate = Math.Min(hexose_inaggregate + 1, 4);
                    for (int hexose_count = 1; hexose_count <= min_hexose_inaggregate; hexose_count++)
                    {
                        GlycanIon hexose_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)hexnac_count, 0, 0, glycan_mass);
                        glycanIons.Add(hexose_glycanIon);

                        if (!extended_fucosylation)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }
                        else if (fuc_count > 0)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                            {
                                GlycanIon add_fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                                glycanIons.Add(add_fuc_glycanIon);
                            }

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }

                        if (iteration_count < xyl_count)
                        {
                            GlycanIon xyl_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 0, 1, glycan_mass);
                            glycanIons.Add(xyl_glycanIon);
                        }

                        if (hexose_count == 3 && hexnac_count >= 2 * core_count && extended)
                        {
                            for (int extra_hexnac_count = 0; extra_hexnac_count < hexnac_inaggregate - hexnac_count + 1; extra_hexnac_count++)
                            {
                                if (extra_hexnac_count + hexnac_count > hexnac_inaggregate)
                                {
                                    continue;
                                }

                                if (extra_hexnac_count > 0)
                                {
                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);

                                    if (!extended_fucosylation)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }
                                    else if (fuc_count > 0)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                                        {
                                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                                            glycanIons.Add(add_fuc_glycanIon);
                                        }

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }

                                    if (iteration_count < xyl_count)
                                    {
                                        GlycanIon xyl_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 0, 1, glycan_mass);
                                        glycanIons.Add(xyl_glycanIon);
                                    }

                                }

                                for (int extra_hexose_count = 1; extra_hexose_count < hexose_inaggregate - hexose_count + 1; extra_hexose_count++)
                                {
                                    if (extra_hexose_count + hexose_count > hexose_inaggregate)
                                    {
                                        continue;
                                    }

                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)(hexose_count + extra_hexose_count), (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);

                                    if (!extended_fucosylation)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }
                                    else if (fuc_count > 0)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                                        {
                                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, (byte)add_fuc_count, 0, glycan_mass);
                                            glycanIons.Add(add_fuc_glycanIon);
                                        }

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }

                                    if (iteration_count < xyl_count)
                                    {
                                        GlycanIon xyl_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 0, 1, glycan_mass);
                                        glycanIons.Add(xyl_glycanIon);
                                    }

                                }
                            }
                        }
                    }

                }


            }

            return glycanIons;
        }

        private static GlycanIon GenerateGlycanIon(byte hexose_count, byte hexnac_count, byte fuc_count, byte xyl_count, int glycan_mass)
        {
            byte[] ionKind = new byte[] { hexose_count, hexnac_count, 0, 0, fuc_count, 0, 0, 0, 0, xyl_count };

            int ionMass = Glycan.GetMass(ionKind);

            GlycanIon glycanIon = new GlycanIon(null, ionMass, ionKind, glycan_mass - ionMass);

            return glycanIon;
        }

        private static GlycanIon ExtendGlycanIon(GlycanIon glycanIon, byte hexose_count, byte hexnac_count, byte fuc_count, byte xyl_count, int glycan_mass)
        {
            byte[] ionKind = glycanIon.IonKind;
            ionKind[0] += hexose_count;
            ionKind[1] += hexnac_count;
            ionKind[4] += fuc_count;
            ionKind[9] += xyl_count;

            int ionMass = Glycan.GetMass(ionKind);

            GlycanIon extend_glycanIon = new GlycanIon(null, ionMass, ionKind, glycan_mass - ionMass);

            return extend_glycanIon;
        }

        //This function build fragments based on the general core of OGlyco fragments. 
        //From https://github.com/mobiusklein/glycopeptidepy/structure/fragmentation_strategy/glycan.py
        //The fragment generation is not as good as structure based method. So it is better to use a structure based O-Glycan database.
        // We don't use this function now, alternatively, we use the 'OGlycanCompositionCombinationChildIons'.
        public static List<GlycanIon> OGlycanCompositionFragments(byte[] kind)
        {
            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int glycan_mass = Glycan.GetMass(kind);

            int iteration_count = 0;
            bool extended = true;

            int fuc_count = kind[4];
            int hexnac_inaggregate = kind[0];
            int hexose_inaggregate = kind[1];

            for (int hexnac_count = 0; hexnac_count < 3; hexnac_count++)
            {
                if (hexnac_inaggregate < hexnac_count)
                {
                    continue;
                }


                if (hexnac_count >= 1)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);

                    glycanIons.Add(glycanIon);

                    if (iteration_count < fuc_count)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);

                        glycanIons.Add(fuc_glycanIon);
                    }

                    for (int hexose_count = 0; hexose_count < 2; hexose_count++)
                    {
                        if (hexose_inaggregate < hexose_count)
                        {
                            continue;
                        }

                        if (hexose_count > 0)
                        {
                            GlycanIon hexose_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)hexnac_count, 0, 0, glycan_mass);
                            glycanIons.Add(hexose_glycanIon);

                            if (iteration_count < fuc_count)
                            {
                                GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);

                                glycanIons.Add(fuc_glycanIon);
                            }
                        }

                        // After the core motif has been exhausted, speculatively add on the remaining core monosaccharides sequentially until exhausted.

                        if (extended && hexnac_inaggregate - hexnac_count >= 0)
                        {
                            for (int extra_hexnac_count = 0; extra_hexnac_count  < hexnac_inaggregate - hexnac_count + 1; extra_hexnac_count ++)
                            {
                                if (extra_hexnac_count > 0)
                                {
                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);


                                    if (iteration_count < fuc_count)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);

                                        glycanIons.Add(fuc_glycanIon);
                                    }

                                }

                                if (hexose_inaggregate > hexose_count && hexose_count > 0)
                                {
                                    for (int extra_hexose_count = 0; extra_hexose_count < hexose_inaggregate - hexose_count; extra_hexose_count++)
                                    {
                                        if (extra_hexose_count > 0 && extra_hexose_count + hexose_count >0)
                                        {

                                            GlycanIon new_glycanIon = GenerateGlycanIon((byte)(hexose_count + extra_hexose_count), (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                            glycanIons.Add(new_glycanIon);


                                            if (iteration_count < fuc_count)
                                            {
                                                GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);

                                                glycanIons.Add(fuc_glycanIon);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }

            }


            return glycanIons;
        }

        /// <summary>
        /// Generate some child ions based on the kind array. The kind array is the combination of the monosaccharides then filter by the rules.
        /// </summary>
        /// <param name="kind"> glycan Kind[]</param>
        /// <returns> The glycanIon collection </returns>
        public static List<GlycanIon> OGlycanCompositionCombinationChildIons(byte[] kind)
        {
            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int glycan_mass = Glycan.GetMass(kind);

            List<byte[]> _kinds = new List<byte[]>();
            HashSet<string> _keys = new HashSet<string>();

            _kinds.Add((byte[])kind.Clone());
            _GetCombinations(kind, _kinds, _keys);

            foreach (var k in _kinds)
            {
                //Rules to build OGlycan child ions. Filter the kind array which doesn't meet the rules.
                //At least one HexNAc
                if (k[1] == 0)
                {
                    continue;
                }

                //#Fucose <= #HexNAc. One Fucose modify one 
                if (k[4]!= 0 && k[4] > k[1] )
                {
                    continue;
                }

                //#NeuAc * 2 >= #Acetylation. One NeuAc can be modified with two Acetylation
                if (k[9]!= 0 && k[2]*2 < k[9])
                {
                    continue;
                }

                var ionMass = Glycan.GetMass(k);
                GlycanIon glycanIon = new GlycanIon(null, ionMass, k, glycan_mass - ionMass);
                glycanIons.Add(glycanIon);
            }

            return glycanIons.OrderBy(p=>p.IonMass).ToList();
        }

        /// <summary>
        /// Try to create all possible combinations from the glycan kind[]. And store the combination array in the _kinds list.
        /// </summary>
        /// <param name="kind"> ex. [2,2,0]</param>
        /// <param name="_kinds"></param>
        /// <param name="_keys"></param>
        private static void _GetCombinations(byte[] kind, List<byte[]> _kinds, HashSet<string> _keys) 
        {                                                                                            
            if (kind.Sum(p=>p) == 0)                                                                 
            {
                return; // if we don't have any monosaccharide, no need to generate the child ions.
            }
            else
            {
                for (int i = 0; i < kind.Length; i++) //traverse the kind array
                {
                    if (kind[i] >= 1)
                    {
                        byte[] akind = (byte[])kind.Clone();
                        akind[i]--;
                        if (akind.Sum(p => p) != 0)
                        {
                            if (!_keys.Contains(Glycan.GetKindString(akind)))
                            {
                                _keys.Add(Glycan.GetKindString(akind));
                                _kinds.Add((byte[])akind.Clone());
                                _GetCombinations(akind, _kinds, _keys);
                            }
                           
                        }

                    }
                }
            }
        }
    }
}
