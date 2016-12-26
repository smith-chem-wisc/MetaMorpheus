using System;
using System.Collections.Generic;
using System.Xml;

namespace MetaMorpheus
{
    public class XMLdb
    {
        public bool Use { get; set; }
        public string FileName { get; private set; }

        public XMLdb(string FileName)
        {
            Use = true;
            this.FileName = FileName;
        }

        public IEnumerable<Protein> getProteins(bool onTheFlyDecoys, IDictionary<string, List<MorpheusModification>> allModifications)
        {
            using (XmlReader xml = XmlReader.Create(FileName))
            {
                string[] nodes = new string[6];

                string dataset = null;
                string accession = null;
                string name = null;
                string full_name = null;
                bool fragment = false;
                string organism = null;
                string gene_name = null;
                string protein_existence = null;
                string sequence_version = null;
                string sequence = null;
                string feature_type = null;
                string feature_description = null;
                int oneBasedfeature_position = -1;
                int oneBasedbeginPosition = -1;
                int oneBasedendPosition = -1;
                List<int> oneBasedBeginPositions = new List<int>();
                List<int> oneBasedEndPositions = new List<int>();
                List<string> peptideTypes = new List<string>();
                Dictionary<int, List<MorpheusModification>> oneBasedModifications = new Dictionary<int, List<MorpheusModification>>();
                int offset = 0;
                while (xml.Read())
                {
                    switch (xml.NodeType)
                    {
                        case XmlNodeType.Element:
                            nodes[xml.Depth] = xml.Name;
                            switch (xml.Name)
                            {
                                case "entry":
                                    dataset = xml.GetAttribute("dataset");
                                    break;

                                case "accession":
                                    if (accession == null)
                                    {
                                        accession = xml.ReadElementString();
                                    }
                                    break;

                                case "name":
                                    if (xml.Depth == 2)
                                    {
                                        name = xml.ReadElementString();
                                    }
                                    else if (nodes[2] == "gene")
                                    {
                                        if (gene_name == null)
                                        {
                                            gene_name = xml.ReadElementString();
                                        }
                                    }
                                    else if (nodes[2] == "organism")
                                    {
                                        if (organism == null)
                                        {
                                            organism = xml.ReadElementString();
                                        }
                                    }
                                    break;

                                case "fullName":
                                    if (full_name == null)
                                    {
                                        full_name = xml.ReadElementString();
                                    }
                                    break;

                                case "proteinExistence":
                                    protein_existence = xml.GetAttribute("type");
                                    break;

                                case "feature":
                                    feature_type = xml.GetAttribute("type");
                                    feature_description = xml.GetAttribute("description");
                                    break;

                                case "position":
                                    oneBasedfeature_position = int.Parse(xml.GetAttribute("position"));
                                    break;

                                case "begin":
                                    try
                                    {
                                        oneBasedbeginPosition = int.Parse(xml.GetAttribute("position"));
                                    }
                                    catch
                                    {
                                    }
                                    break;

                                case "end":
                                    try
                                    {
                                        oneBasedendPosition = int.Parse(xml.GetAttribute("position"));
                                    }
                                    catch
                                    {
                                    }
                                    break;

                                case "sequence":
                                    fragment = xml.GetAttribute("fragment") != null;
                                    sequence_version = xml.GetAttribute("version");
                                    sequence = xml.ReadElementString().Replace("\n", null);
                                    break;
                            }
                            break;

                        case XmlNodeType.EndElement:
                            switch (xml.Name)
                            {
                                case "feature":
                                    if (feature_type == "modified residue" && allModifications != null && !feature_description.Contains("variant") && allModifications.ContainsKey(feature_description))
                                    {
                                        List<MorpheusModification> residue_modifications;
                                        if (!oneBasedModifications.TryGetValue(oneBasedfeature_position, out residue_modifications))
                                        {
                                            residue_modifications = new List<MorpheusModification>();
                                            oneBasedModifications.Add(oneBasedfeature_position, residue_modifications);
                                        }
                                        int semicolon_index = feature_description.IndexOf(';');
                                        if (semicolon_index >= 0)
                                        {
                                            feature_description = feature_description.Substring(0, semicolon_index);
                                        }
                                        residue_modifications.AddRange(allModifications[feature_description]);
                                    }
                                    else if ((feature_type == "peptide" || feature_type == "signal peptide" || feature_type == "propeptide" || feature_type == "chain") && oneBasedbeginPosition >= 0 && oneBasedendPosition >= 0)
                                    {
                                        oneBasedBeginPositions.Add(oneBasedbeginPosition);
                                        oneBasedEndPositions.Add(oneBasedendPosition);
                                        peptideTypes.Add(feature_type);
                                    }
                                    oneBasedbeginPosition = -1;
                                    oneBasedendPosition = -1;

                                    break;

                                case "entry":
                                    string dataset_abbreviation;
                                    if (dataset != null && dataset.Equals("Swiss-Prot", StringComparison.InvariantCultureIgnoreCase))
                                        dataset_abbreviation = "sp";
                                    else
                                        dataset_abbreviation = "uk";

                                    if (!SpecificProteinSelection.enabled || SpecificProteinSelection.ConsiderProtein(accession))
                                    {
                                        Protein protein = new Protein(sequence, accession, dataset_abbreviation, oneBasedModifications, oneBasedBeginPositions.ToArray(), oneBasedEndPositions.ToArray(), peptideTypes.ToArray(), name, full_name, offset, false);

                                        yield return protein;

                                        offset += protein.Length;
                                        if (onTheFlyDecoys)
                                        {
                                            char[] sequence_array = sequence.ToCharArray();
                                            Dictionary<int, List<MorpheusModification>> decoy_modifications = null;
                                            if (sequence.StartsWith("M"))
                                            {
                                                // Do not include the initiator methionine in reversal!!!
                                                Array.Reverse(sequence_array, 1, sequence.Length - 1);
                                                if (oneBasedModifications != null)
                                                {
                                                    decoy_modifications = new Dictionary<int, List<MorpheusModification>>(oneBasedModifications.Count);
                                                    foreach (KeyValuePair<int, List<MorpheusModification>> kvp in oneBasedModifications)
                                                    {
                                                        if (kvp.Key == 1)
                                                        {
                                                            decoy_modifications.Add(1, kvp.Value);
                                                        }
                                                        else if (kvp.Key > 1)
                                                        {
                                                            decoy_modifications.Add(sequence.Length - kvp.Key + 2, kvp.Value);
                                                        }
                                                    }
                                                }
                                            }
                                            else
                                            {
                                                Array.Reverse(sequence_array);
                                                if (oneBasedModifications != null)
                                                {
                                                    decoy_modifications = new Dictionary<int, List<MorpheusModification>>(oneBasedModifications.Count);
                                                    foreach (KeyValuePair<int, List<MorpheusModification>> kvp in oneBasedModifications)
                                                    {
                                                        decoy_modifications.Add(sequence.Length - kvp.Key + 1, kvp.Value);
                                                    }
                                                }
                                            }
                                            string reversed_sequence = new string(sequence_array);
                                            int[] decoybeginPositions = new int[oneBasedBeginPositions.Count];
                                            int[] decoyendPositions = new int[oneBasedEndPositions.Count];
                                            string[] decoyBigPeptideTypes = new string[oneBasedEndPositions.Count];
                                            for (int i = 0; i < decoybeginPositions.Length; i++)
                                            {
                                                decoybeginPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedEndPositions[i] + 1;
                                                decoyendPositions[oneBasedBeginPositions.Count - i - 1] = sequence.Length - oneBasedBeginPositions[i] + 1;
                                                decoyBigPeptideTypes[oneBasedBeginPositions.Count - i - 1] = peptideTypes[i];
                                            }
                                            Protein decoy_protein = new Protein(reversed_sequence, "DECOY_" + accession, dataset_abbreviation, decoy_modifications, decoybeginPositions, decoyendPositions, decoyBigPeptideTypes, name, full_name, offset, true);
                                            yield return decoy_protein;
                                            offset += protein.Length;
                                        }
                                    }
                                    dataset = null;
                                    accession = null;
                                    name = null;
                                    full_name = null;
                                    fragment = false;
                                    organism = null;
                                    gene_name = null;
                                    protein_existence = null;
                                    sequence_version = null;
                                    sequence = null;
                                    feature_type = null;
                                    feature_description = null;
                                    oneBasedfeature_position = -1;
                                    oneBasedModifications = new Dictionary<int, List<MorpheusModification>>();

                                    oneBasedBeginPositions = new List<int>();
                                    oneBasedEndPositions = new List<int>();
                                    peptideTypes = new List<string>();
                                    break;
                            }
                            break;
                    }
                }
            }
        }
    }
}