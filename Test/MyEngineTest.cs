using InternalLogicEngineLayer;
using NUnit.Framework;
using System.Text;

namespace Test
{
    [TestFixture]
    public class MyEngineTest
    {
        #region Public Methods

        [Test]
        public void TestMyEngine()
        {
            MyEngine level0engine = new TestEngine(0, null);
            Assert.That(() => level0engine.Run(),
            Throws.TypeOf<EngineValidationException>()
                .With.Property("Message").EqualTo("param1 cannot be null"));

            level0engine = new TestEngine(0, new object());
            level0engine.Run();
        }

        #endregion Public Methods

        #region Private Classes

        private class TestEngine : MyEngine
        {
            #region Private Fields

            private object param1;

            #endregion Private Fields

            #region Public Constructors

            public TestEngine(int level, object param1) : base(level)
            {
                this.param1 = param1;
            }

            #endregion Public Constructors

            #region Protected Methods

            protected override MyResults RunSpecific()
            {
                return new TestResults(this);
            }

            protected override void ValidateParams()
            {
                if (param1 == null)
                    throw new EngineValidationException("param1 cannot be null");
            }

            #endregion Protected Methods

            #region Private Classes

            private class TestResults : MyResults
            {
                #region Public Constructors

                public TestResults(MyEngine e) : base(e)
                {
                }

                #endregion Public Constructors

                #region Protected Methods

                protected override string GetStringForOutput()
                {
                    var sb = new StringBuilder();
                    sb.Append("String for the TestResults results class");
                    return sb.ToString();
                }

                #endregion Protected Methods
            }

            #endregion Private Classes
        }

        #endregion Private Classes
    }
}