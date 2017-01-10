using MetaMorpheusLogic;
using System;

namespace MetaMorpheusGUI
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
                outProgressBar.IsIndeterminate = true;
                statusLabel.Content = s;
            }
        }

        private void NewoutProgressBar(object sender, ProgressEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutProgressBar(sender, s)));
            }
            else
            {
                outProgressBar.IsIndeterminate = false;
                outProgressBar.Value = s.new_progress;
                statusLabel.Content = s.v;
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
                //TODO: Check those
                XMLdbPanel.IsEnabled = false;
                DatafilesStackPanel.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;
                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                tasksPanel.IsEnabled = false;

                statusLabel.Content = "Starting all tasks...";
                outProgressBar.IsIndeterminate = true;

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
                //TODO: Check those
                XMLdbPanel.IsEnabled = true;
                DatafilesStackPanel.IsEnabled = true;
                addSearchTaskButton.IsEnabled = true;
                addCalibrateTaskButton.IsEnabled = true;
                addGPTMDTaskButton.IsEnabled = true;
                tasksPanel.IsEnabled = true;

                statusLabel.Content = "Finished all tasks!";
                outProgressBar.IsIndeterminate = false;
                outProgressBar.Value = 100;

                dataGridDatafiles.Items.Refresh();
            }
        }

        private void NewSuccessfullyFinishedWritingFile(object sender, SingleFileEventArgs v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedWritingFile(sender, v)));
            }
            else
            {
                addFinishedFile(v.writtenFile);
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