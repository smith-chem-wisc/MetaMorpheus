using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//Object that is placed in the DataGrid for ModSelection when writing a pruned DB
namespace MetaMorpheusGUI
{

    public class ModTypeForGrid
    {

        public ModTypeForGrid(string modName, bool item2, bool item3, bool item4, bool item5)
        {
            ModName = modName;
            Item2 = item2;
            Item3 = item3;
            Item4 = item4;
            Item5 = item5;

        }

        //types
        public string ModName { get; set; }

        public bool Item2 { get; set; }

        public bool Item3 { get; set; }

        public bool Item4 { get; set; }

        public bool Item5 {get; set; }





    }
}
