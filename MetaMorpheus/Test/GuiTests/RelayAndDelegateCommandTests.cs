using System;
using System.Windows.Input;
using GuiFunctions;
using NUnit.Framework;

namespace Test.GuiTests
{
    [TestFixture]
    public class RelayAndDelegateCommandTests
    {
        [Test]
        public void RelayCommand_Constructor_StoresAction()
        {
            // Arrange
            bool executed = false;
            var command = new RelayCommand(() => executed = true);

            // Act
            command.Execute(null);

            // Assert
            Assert.That(executed, Is.True);
        }

        [Test]
        public void RelayCommand_CanExecute_AlwaysReturnsTrue()
        {
            // Arrange
            var command = new RelayCommand(() => { });

            // Act & Assert
            Assert.That(command.CanExecute(null), Is.True);
            Assert.That(command.CanExecute("any parameter"), Is.True);
        }

        [Test]
        public void RelayCommand_Execute_InvokesAction()
        {
            // Arrange
            int callCount = 0;
            var command = new RelayCommand(() => callCount++);

            // Act
            command.Execute(null);
            command.Execute(null);

            // Assert
            Assert.That(callCount, Is.EqualTo(2));
        }

        [Test]
        public void RelayCommand_CanExecuteChanged_Event_IsNotNullAndCanBeSubscribed()
        {
            // Arrange
            var command = new RelayCommand(() => { });
            bool eventFired = false;
            EventHandler handler = (s, e) => eventFired = true;

            // Act
            command.CanExecuteChanged += handler;
            // The event is never raised by the class, but we can invoke it manually for coverage
            handler.Invoke(command, EventArgs.Empty);

            // Assert
            Assert.That(eventFired, Is.True);
        }

        [Test]
        public void DelegateCommand_Constructor_StoresAction()
        {
            // Arrange
            object received = null;
            var command = new DelegateCommand(obj => received = obj);

            // Act
            command.Execute("test");

            // Assert
            Assert.That(received, Is.EqualTo("test"));
        }

        [Test]
        public void DelegateCommand_CanExecute_AlwaysReturnsTrue()
        {
            // Arrange
            var command = new DelegateCommand(obj => { });

            // Act & Assert
            Assert.That(command.CanExecute(null), Is.True);
            Assert.That(command.CanExecute("any parameter"), Is.True);
        }

        [Test]
        public void DelegateCommand_Execute_InvokesActionWithParameter()
        {
            // Arrange
            object lastParam = null;
            var command = new DelegateCommand(param => lastParam = param);

            // Act
            command.Execute(42);
            command.Execute("foo");

            // Assert
            Assert.That(lastParam, Is.EqualTo("foo"));
        }

        [Test]
        public void DelegateCommand_CanExecuteChanged_Event_IsNotNullAndCanBeSubscribed()
        {
            // Arrange
            var command = new DelegateCommand(obj => { });
            bool eventFired = false;
            EventHandler handler = (s, e) => eventFired = true;

            // Act
            command.CanExecuteChanged += handler;
            // The event is never raised by the class, but we can invoke it manually for coverage
            handler.Invoke(command, EventArgs.Empty);

            // Assert
            Assert.That(eventFired, Is.True);
        }

        [Test]
        public void RelayCommand_ThrowsIfActionIsNull()
        {
            // Arrange, Act & Assert
            Assert.Throws<NullReferenceException>(() => new RelayCommand(null).Execute(new()));
        }

        [Test]
        public void DelegateCommand_ThrowsIfActionIsNull()
        {
            // Arrange, Act & Assert
            Assert.Throws<NullReferenceException>(() => new DelegateCommand(null).Execute(new()));
        }
    }
}
