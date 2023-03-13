using System.Windows.Controls;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// This text box requires input text to be integer only.
    /// </summary>
    public class IntegerTexBoxControl : TextBox
    {
        protected override void OnPreviewTextInput(TextCompositionEventArgs e)
        {
            foreach (var character in e.Text)
            {
                if (!char.IsDigit(character))
                {
                    e.Handled = true;
                    return;
                }
            }
            e.Handled = false;
        }
    }
}
