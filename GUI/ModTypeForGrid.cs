//Object that is placed in the DataGrid for ModSelection when writing a pruned DB
namespace MetaMorpheusGUI
{
    public class ModTypeForGrid
    {
        #region Public Constructors

        public ModTypeForGrid(string modName)
        {
            ModName = modName;
            Item2 = true;
        }

        #endregion Public Constructors

        #region Public Properties

        //types
        public string ModName { get; set; }

        public bool Item2 { get; set; }

        public bool Item3 { get; set; }

        public bool Item4 { get; set; }

        public bool Item5 { get; set; }

        #endregion Public Properties
    }
}