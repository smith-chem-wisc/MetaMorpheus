using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// This text box requires input text to be decimal only.
    /// </summary>
    public class DoubleTextBoxControl : TextBox
    {
        public DoubleTextBoxControl()
        {
            HorizontalContentAlignment = HorizontalAlignment.Center;
            VerticalContentAlignment = VerticalAlignment.Center;
        }

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

        /// <summary>
        /// Cursor is removed from text box on pressing Return
        /// </summary>
        /// <param name="e"></param>
        protected override void OnKeyDown(KeyEventArgs e)
        {
            base.OnKeyDown(e);
            if (e.Key == Key.Return)
                Keyboard.ClearFocus();
        }
    }
}