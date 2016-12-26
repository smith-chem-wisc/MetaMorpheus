using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MetaMorpheus
{
    public class ModList
    {
        private List<MorpheusModification> mods;
        public bool Localize { get; set; }
        public bool Fixed { get; set; }
        public bool Variable { get; set; }

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

        public ModList(string fileName, bool Localize, bool Variable, bool Fixed)
        {
            this.FullFileName = Path.GetFullPath(fileName);
            this.Localize = Localize;
            this.Variable = Variable;
            this.Fixed = Fixed;
            mods = ProteomeDatabaseReader.ReadModFile(FullFileName).ToList();
            Description = File.ReadLines(FullFileName).First();
        }

        public List<MorpheusModification> getMods()
        {
            return mods;
        }

        private string FullFileName;
    }
}