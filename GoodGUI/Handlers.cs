using System;
using System.Diagnostics;

namespace GoodGUI
{
    public partial class MainWindow
    {
        private void NewoutLabelStatus(object sender, string s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutLabelStatus(sender, s)));
            }
            else
            {
                outLabelStatus.Content = s;
            }
        }

        private void NewoutTextBox(object sender, string s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutTextBox(sender, s)));
            }
            else
            {
                outTextBox.Text = s;
            }
        }

        private void NewoutProgressBar(object sender, int s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutProgressBar(sender, s)));
            }
            else
            {
                Debug.Assert(s >= 0 && s <= 100);
                outProgressBar.Value = s;
            }
        }

        private void NewoutRichTextBox(object sender, string tup)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutRichTextBox(sender, tup)));
            }
            else
            {
                RegOutput(tup);
            }
        }

        private void NewRefreshBetweenTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewRefreshBetweenTasks(sender, e)));
            }
            else
            {
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        private void NewSuccessfullyStartingAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyStartingAllTasks(sender, e)));
            }
            else
            {
                XMLdbPanel.IsEnabled = false;
                DatafilesStackPanel.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;
                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                tasksPanel.IsEnabled = false;
                dataGridDatafiles.Items.Refresh();
            }
        }

        private void NewSuccessfullyFinishedAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedAllTasks(sender, e)));
            }
            else
            {
                XMLdbPanel.IsEnabled = true;
                DatafilesStackPanel.IsEnabled = true;
                addSearchTaskButton.IsEnabled = true;
                addCalibrateTaskButton.IsEnabled = true;
                addGPTMDTaskButton.IsEnabled = true;
                tasksPanel.IsEnabled = true;
                dataGridDatafiles.Items.Refresh();
            }
        }

        private void NewSuccessfullyFinishedFile(object sender, string v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedFile(sender, v)));
            }
            else
            {
                addFile(v);
            }
        }

        private void NewUnSuccessfullyFinishedTask(object sender, string v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewUnSuccessfullyFinishedTask(sender, v)));
            }
            else
            {
                ErrorOutput(v);
            }
        }
    }
}