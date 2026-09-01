using System.Windows.Controls;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MassDiffAcceptorControl.xaml
    /// </summary>
    public partial class MassDifferenceAcceptorControl : UserControl
    {
        public MassDifferenceAcceptorControl()
        {
            InitializeComponent();
        }

        private void TextBox_KeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                var textBox = sender as TextBox;
                if (textBox != null)
                {
                    // Force update of binding source
                    var binding = textBox.GetBindingExpression(TextBox.TextProperty);
                    binding?.UpdateSource();
                }
                Keyboard.ClearFocus();
                e.Handled = true;
            }
        }
    }
}
