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

        public static readonly DependencyProperty LowerBoundProperty =
            DependencyProperty.Register(
                nameof(LowerBound),
                typeof(double),
                typeof(DoubleTextBoxControl),
                new PropertyMetadata(double.MinValue));

        public static readonly DependencyProperty UpperBoundProperty =
            DependencyProperty.Register(
                nameof(UpperBound),
                typeof(double),
                typeof(DoubleTextBoxControl),
                new PropertyMetadata(double.MaxValue));

        public double LowerBound
        {
            get => (double)GetValue(LowerBoundProperty);
            set => SetValue(LowerBoundProperty, value);
        }

        public double UpperBound
        {
            get => (double)GetValue(UpperBoundProperty);
            set => SetValue(UpperBoundProperty, value);
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

        protected override void OnTextChanged(TextChangedEventArgs e)
        {
            base.OnTextChanged(e);

            if (double.TryParse(Text, out double value))
            {
                if (value < LowerBound)
                    Text = LowerBound.ToString();
                else if (value > UpperBound)
                    Text = UpperBound.ToString();
            }
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
