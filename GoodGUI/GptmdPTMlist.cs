using MetaMorpheus;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace GoodGUI
{
    internal class GptmdPTMlist
    {
        public bool Use { get; set; }

        public string FileName
        {
            get
            {
                return Path.GetFileName(FullFileName);
            }
        }

        public int Count
        {
            get { return mods.Count(); }
        }

        public string Description { get; private set; }

        public GptmdPTMlist(string v1, bool v2)
        {
            this.FullFileName = Path.GetFullPath(v1);
            this.Use = v2;
            mods = ProteomeDatabaseReader.ReadModFile(FullFileName).ToList();
            Description = File.ReadLines(FullFileName).First();
        }

        private List<MorpheusModification> mods;

        public List<MorpheusModification> getMods()
        {
            return mods;
        }

        private string FullFileName;
    }
}