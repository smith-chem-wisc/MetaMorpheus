//Object that is placed in the DataGrid for ModSelection when writing a pruned DB
namespace MetaMorpheusGUI
{
    public class ModTypeForGrid
    {
        public ModTypeForGrid(string modName)
        {
            ModName = modName;
            Item2 = true;
        }

        //types
        public string ModName { get; set; }

        public bool Item2 { get; set; }

        public bool Item3 { get; set; }

        public bool Item4 { get; set; }

        public bool Item5 { get; set; }
    }
}