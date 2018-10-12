using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace MetaMorpheusGUI
{
    public partial class CustomMsgBox : Form
    {
        public CustomMsgBox()
        {
            InitializeComponent();
            this.button1.Click += new EventHandler(button1_Click);
            this.button2.Click += new EventHandler(button2_Click);
            this.button3.Click += new EventHandler(button3_Click);
        }

        //internal GuiGlobalParams GuiGlobalParams = new GuiGlobalParams();
        static CustomMsgBox MsgBox;
        static DialogResult result = DialogResult.No;
        public static DialogResult Show(string text, string caption, string btnYes, string btnNo, string btnStop)
        {
            MsgBox = new CustomMsgBox();
            MsgBox.label1.Text = text;
            MsgBox.button1.Text = btnYes;
            MsgBox.button2.Text = btnNo;
            MsgBox.button3.Text = btnStop;
            MsgBox.Text = caption;
            MsgBox.ShowDialog();
            return result;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            result = DialogResult.Yes;
            MsgBox.Close();
        }

        private void button2_Click(object sender, EventArgs e)
        {
            result = DialogResult.No;
            MsgBox.Close();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            result = DialogResult.OK;
            MsgBox.Close();
        }
    }
}
