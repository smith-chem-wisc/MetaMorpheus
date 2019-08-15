using EngineLayer;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    class ProteinForTreeView
    {
        public ProteinForTreeView(string displayName, string accession, Dictionary<string, PeptideForTreeView> uniquePep, Dictionary<string, PeptideForTreeView> sharedPep, Dictionary<string, PeptideForTreeView> allPep)
        {
            Children = new ObservableCollection<PeptideForTreeView>();
            Expanded = false;
            DisplayName = displayName;
            Accession = accession;
            UniquePeptides = uniquePep;
            SharedPeptides = sharedPep;
            AllPeptides = allPep;
        }

        public Dictionary<string, PeptideForTreeView> UniquePeptides { get; }
        public Dictionary<string, PeptideForTreeView> SharedPeptides { get; }
        public Dictionary<string, PeptideForTreeView> AllPeptides { get; } 
        public string DisplayName { get; }
        public string Accession { get; }
        public ObservableCollection<PeptideForTreeView> Children { get; set; }
        public bool Expanded { get; set; }
    }

    class PeptideForTreeView
    {
        public PeptideForTreeView(string displayName, ProteinForTreeView parent)
        {
            Children = new ObservableCollection<PsmForTreeView>();
            Parent = parent;
            DisplayName = displayName;
            Expanded = false;
        }

        public ProteinForTreeView Parent { get; }
        public string DisplayName { get; }
        public ObservableCollection<PsmForTreeView> Children { get; set; }
        public bool Expanded { get; set; }


    }

    class PsmForTreeView
    {
        public PsmForTreeView(int scan, PsmFromTsv psm, PeptideForTreeView parent)
        {
            Parent = parent;
            Psm = psm;
            ScanNo = scan;
        }

        public PeptideForTreeView Parent { get; }
        public PsmFromTsv Psm { get; }
        public int ScanNo { get; }
    }

    //class Node
    //{
    //    public Node()
    //    {
    //        Children = new ObservableCollection<ProteinForTreeView>();
    //    }

    //    public string Name { get; set; }
    //    public ObservableCollection<ProteinForTreeView> Children { get; set; }
    //}
}
