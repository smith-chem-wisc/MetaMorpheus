using OldInternalLogic;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
    public class ModListForGPTMD
    {
        public bool Use { get; set; }
        public bool Fixed { get; set; }
        public bool Variable { get; set; }
        public bool Localize { get; set; }

        public string FileName
        {
            get
            {
                return uu.FileName;
            }
        }

        public string Description { get { return uu.Description; } }

        private ModList uu;

        public ModListForGPTMD(ModList uu)
        {
            this.uu = uu;
        }

        public IEnumerable<MorpheusModification> getMods()
        {
            return uu.getMods();
        }
    }
}