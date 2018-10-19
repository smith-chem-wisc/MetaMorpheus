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
            this.button1.Click += new EventHandler(Button1_Click);
            this.button2.Click += new EventHandler(Button2_Click);
            this.button3.Click += new EventHandler(Button3_Click);
        }

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

        private void Button1_Click(object sender, EventArgs e)
        {
            result = DialogResult.Yes;
            MsgBox.Close();
        }

        private void Button2_Click(object sender, EventArgs e)
        {
            result = DialogResult.No;
            MsgBox.Close();
        }

        private void Button3_Click(object sender, EventArgs e)
        {
            result = DialogResult.OK;
            MsgBox.Close();
        }
    }
}
