using EngineLayer;
using NUnit.Framework;
using System.Text;

namespace Test
{
    [TestFixture]
    public class MyEngineTest
    {

        #region Public Methods

        [Test]
        public static void TestMyEngine()
        {
            MetaMorpheusEngine level0engine = new TestEngine(0);

            level0engine = new TestEngine(0);
            level0engine.Run();
        }

        #endregion Public Methods

        #region Private Classes

        private class TestEngine : MetaMorpheusEngine
        {

            #region Public Constructors

            public TestEngine(int level)
            {
            }

            #endregion Public Constructors

            #region Protected Methods

            protected override MetaMorpheusEngineResults RunSpecific()
            {
                return new TestResults(this);
            }

            #endregion Protected Methods

            #region Private Classes

            private class TestResults : MetaMorpheusEngineResults
            {

                #region Public Constructors

                public TestResults(MetaMorpheusEngine e) : base(e)
                {
                }

                #endregion Public Constructors

                #region Public Methods

                public override string ToString()
                {
                    var sb = new StringBuilder();
                    sb.AppendLine(base.ToString());
                    sb.Append("String for the TestResults results class");
                    return sb.ToString();
                }

                #endregion Public Methods

            }

            #endregion Private Classes

        }

        #endregion Private Classes

    }
}