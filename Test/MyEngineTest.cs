using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.Text;

namespace Test
{
    [TestFixture]
    public static class MyEngineTest
    {

        [Test]
        public static void TestMyEngine()
        {
            MetaMorpheusEngine level0engine = new TestEngine(0);

            level0engine = new TestEngine(0);
            level0engine.Run();
        }

        private class TestEngine : MetaMorpheusEngine
        {

            public TestEngine(int level) : base(new CommonParameters(), new List<string>())
            {
            }

            protected override MetaMorpheusEngineResults RunSpecific()
            {
                return new TestResults(this);
            }

            private class TestResults : MetaMorpheusEngineResults
            {

                public TestResults(MetaMorpheusEngine e) : base(e)
                {
                }

                public override string ToString()
                {
                    var sb = new StringBuilder();
                    sb.AppendLine(base.ToString());
                    sb.Append("String for the TestResults results class");
                    return sb.ToString();
                }

            }

        }

    }
}