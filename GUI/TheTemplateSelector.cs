using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using TaskLayer;

namespace MetaMorpheusGUI
{
    class TheTemplateSelector : DataTemplateSelector
    {
        public override System.Windows.DataTemplate SelectTemplate(object item, System.Windows.DependencyObject container)
        {
            if (item is Parameter)
            {
                Parameter settings = item as Parameter;
                if (settings.ValueType == "ComboBoxProtease")
                {
                    return ComboBoxProtease;

                }
                else if (settings.ValueType == "Bool")
                {
                    return Bool;
                }
                else if (settings.ValueType == "ComboBoxInit")
                {
                    return ComboBoxInit;
                }
                else if (settings.ValueType == "ProductMassToleranceList")
                {
                    return ComboBoxTolerance;
                }
                else if (settings.ValueType == "TextBox")
                {
                    return TextBox;
                }


            }
            return null;
        }
        public DataTemplate ComboBoxProtease { get; set; }
        public DataTemplate ComboBoxInit { get; set; }
        public DataTemplate ComboBoxTolerance { get; set; }
        public DataTemplate Bool { get; set; }
        public DataTemplate TextBox { get; set; }
    }
}
