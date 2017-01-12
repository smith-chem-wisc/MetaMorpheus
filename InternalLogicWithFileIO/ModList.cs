using OldInternalLogic;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InternalLogicTaskLayer
{
    public class ModList
    {
        private List<MorpheusModification> mods;

        public string FileName
        {
            get
            {
                return Path.GetFileName(FullFileName);
            }
        }

        public int Count
        {
            get { return mods.Count; }
        }

        public string Description { get; private set; }

        public ModList(string fileName)
        {
            this.FullFileName = Path.GetFullPath(fileName);
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