using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    public class DoubleTextBox : TextBox
    {
        protected override void OnPreviewTextInput(TextCompositionEventArgs e)
        {
            foreach (var character in e.Text)
            {
                if (!char.IsDigit(character) && !(character == '.'))
                {
                    e.Handled = true;
                    return;
                }

                if (((TextBox)e.Source).Text.Contains('.') && character == '.')
                {
                    e.Handled = true;
                    return;
                }
            }
            e.Handled = false;
        }
    }
}
