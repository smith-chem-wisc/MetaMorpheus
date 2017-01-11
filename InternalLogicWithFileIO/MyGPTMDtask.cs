using InternalLogic;
using OldInternalLogic;
using System;
using System.Collections.Generic;

namespace InternalLogicWithFileIO
{
    public class MyGPTMDtask : MyTask
    {
        public MyGPTMDtask() : base(2)
        {
        }

        public override MyTaskResults DoTask(AllTasksParams po)
        {
            throw new NotImplementedException();// For each PSM, look at modifications
            //string line;
            //Dictionary<string, HashSet<Tuple<int, string>>> Mods = new Dictionary<string, HashSet<Tuple<int, string>>>();
            //OnStatusUpdate("Reading PSMs...");
            //OnOutput("Reading PSMs...");

            //int modsAdded = 0;
            ////using (StreamReader file = new StreamReader(ok.psmsTSVName))
            ////{
            ////    file.ReadLine();
            ////    while ((line = file.ReadLine()) != null)
            ////    {
            ////        string[] parts = line.Split('\t');
            ////        double massDiff = Convert.ToDouble(parts[19]);
            ////        foreach (MorpheusModification mod in GetMod(massDiff, isotopeErrors, allmods, combos))
            ////        {
            ////            string proteinRange = parts[14].Split('|')[1];
            ////            int proteinLength = Convert.ToInt32(new string(proteinRange.TakeWhile(char.IsDigit).ToArray()));
            ////            var proteinAcession = parts[13].Split('|')[0];
            ////            for (int i = 0; i < parts[12].Length; i++)
            ////            {
            ////                int indexInProtein = Convert.ToInt32(parts[15]) + i;

            ////                if (ModFits(mod, parts[12][i], i > 0 ? parts[12][i] : parts[11][0], i + 1, parts[12].Length, indexInProtein, proteinLength))
            ////                {
            ////                    if (!Mods.ContainsKey(proteinAcession))
            ////                        Mods[proteinAcession] = new HashSet<Tuple<int, string>>();
            ////                    var theTuple = new Tuple<int, string>(indexInProtein, mod.NameInXML);
            ////                    if (!Mods[proteinAcession].Contains(theTuple))
            ////                    {
            ////                        Mods[proteinAcession].Add(theTuple);
            ////                        modsAdded++;
            ////                    }
            ////                }
            ////            }
            ////        }
            ////    }
            ////}

            //OnOutput("Modifications added = " + modsAdded);

            //XmlWriterSettings xmlWriterSettings = new XmlWriterSettings()
            //{
            //    Indent = true,
            //    IndentChars = "  "
            //};

            //OnStatusUpdate("Writing XML");
            //OnOutput("Writing XML...");
            //using (XmlWriter writer = XmlWriter.Create(outputFileName, xmlWriterSettings))
            //{
            //    writer.WriteStartDocument();
            //    writer.WriteStartElement("uniprot");

            //    foreach (Protein protein in allProteins)
            //    {
            //        writer.WriteStartElement("entry");
            //        writer.WriteAttributeString("dataset", protein.dataset_abbreviation);
            //        writer.WriteStartElement("accession");
            //        writer.WriteString(protein.Accession);
            //        writer.WriteEndElement();
            //        writer.WriteStartElement("name");
            //        writer.WriteString(protein.name);
            //        writer.WriteEndElement();

            //        writer.WriteStartElement("protein");
            //        writer.WriteStartElement("recommendedName");
            //        writer.WriteStartElement("fullName");
            //        writer.WriteString(protein.fullName);
            //        writer.WriteEndElement();
            //        writer.WriteEndElement();
            //        writer.WriteEndElement();

            //        for (int i = 0; i < protein.bigPeptideTypes.Count(); i++)
            //        {
            //            writer.WriteStartElement("feature");
            //            writer.WriteAttributeString("type", protein.bigPeptideTypes[i]);
            //            writer.WriteStartElement("location");
            //            writer.WriteStartElement("begin");
            //            writer.WriteAttributeString("position", protein.oneBasedBeginPositions[i].ToString());
            //            writer.WriteEndElement();
            //            writer.WriteStartElement("end");
            //            writer.WriteAttributeString("position", protein.oneBasedEndPositions[i].ToString());
            //            writer.WriteEndElement();
            //            writer.WriteEndElement();
            //            writer.WriteEndElement();
            //        }

            //        IEnumerable<Tuple<int, string>> SortedMods = protein.OneBasedPossibleLocalizedModifications.SelectMany(
            //            b => b.Value.Select(c => new Tuple<int, string>(b.Key, c.NameInXML)
            //            ));
            //        IEnumerable<Tuple<int, string>> FinalSortedMods;
            //        if (Mods.ContainsKey(protein.Accession))
            //            FinalSortedMods = SortedMods.Union(Mods[protein.Accession]).OrderBy(b => b.Item1);
            //        else
            //            FinalSortedMods = SortedMods.OrderBy(b => b.Item1);
            //        foreach (var ye in FinalSortedMods)
            //        {
            //            writer.WriteStartElement("feature");
            //            writer.WriteAttributeString("type", "modified residue");
            //            writer.WriteAttributeString("description", ye.Item2);
            //            writer.WriteStartElement("location");
            //            writer.WriteStartElement("position");
            //            writer.WriteAttributeString("position", ye.Item1.ToString());
            //            writer.WriteEndElement();
            //            writer.WriteEndElement();
            //            writer.WriteEndElement();
            //        }

            //        writer.WriteStartElement("sequence");
            //        writer.WriteAttributeString("length", protein.Length.ToString());
            //        writer.WriteString(protein.BaseSequence);
            //        writer.WriteEndElement();

            //        writer.WriteEndElement();
            //    }

            //    writer.WriteEndElement();
            //    writer.WriteEndDocument();
            //}
        }

        private static bool ModFits(MorpheusModification attemptToLocalize, char v1, char prevAA, int peptideIndex, int peptideLength, int proteinIndex, int proteinLength)
        {
            if (!attemptToLocalize.AminoAcid.Equals('\0') && !attemptToLocalize.AminoAcid.Equals(v1))
                return false;
            if (!attemptToLocalize.PrevAminoAcid.Equals('\0') && !attemptToLocalize.PrevAminoAcid.Equals(prevAA))
                return false;
            if (attemptToLocalize.Type == ModificationType.ProteinNTerminus &&
                ((proteinIndex > 2) || (proteinIndex == 2 && prevAA != 'M')))
                return false;
            if (attemptToLocalize.Type == ModificationType.PeptideNTerminus && peptideIndex > 1)
                return false;
            if (attemptToLocalize.Type == ModificationType.PeptideCTerminus && peptideIndex < peptideLength)
                return false;
            if (attemptToLocalize.Type == ModificationType.ProteinCTerminus && proteinIndex < proteinLength)
                return false;
            return true;
        }

        private static IEnumerable<MorpheusModification> GetMod(double massDiff, bool isotopeErrors, IEnumerable<MorpheusModification> allMods, IEnumerable<Tuple<double, double>> combos, double tol)
        {
            foreach (var Mod in allMods)
            {
                if (Mod.MonoisotopicMassShift > massDiff - tol && Mod.MonoisotopicMassShift < massDiff + tol)
                    yield return Mod;
                if (isotopeErrors && Mod.MonoisotopicMassShift > massDiff - tol - 1.003 && Mod.MonoisotopicMassShift < massDiff + tol - 1.003)
                    yield return Mod;
                if (!double.IsNaN(Mod.AlternativeMassShift) && Mod.AlternativeMassShift > massDiff - tol && Mod.AlternativeMassShift < massDiff + tol)
                    yield return Mod;
                if (!double.IsNaN(Mod.AlternativeMassShift) && isotopeErrors && Mod.AlternativeMassShift > massDiff - tol - 1.003 && Mod.AlternativeMassShift < massDiff + tol - 1.003)
                    yield return Mod;
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (combined > massDiff - tol && combined < massDiff + tol)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                }
                if (isotopeErrors && combined > massDiff - tol - 1.003 && combined < massDiff + tol - 1.003)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                }
            }
        }
    }
}