using InternalLogicEngineLayer;
using NUnit.Framework;
using OldInternalLogic;
using System.Collections.Generic;
using System.Linq;
using System;
using System.Text;

namespace Test
{
    [TestFixture]
    public class MyEngineTest
    {
        [Test]
        public void TestMyEngine()
        {
            MyEngine level0engine = new TestEngine(0, null);
            Assert.That(() => level0engine.Run(),
            Throws.TypeOf<EngineValidationException>()
                .With.Property("Message").EqualTo("param1 cannot be null"));

            level0engine = new TestEngine(0, new object());
            var myResults = level0engine.Run();
            Console.WriteLine(myResults.ToString());
        }

        private class TestEngine : MyEngine
        {
            private object param1;

            public TestEngine(int level, object param1) : base(level)
            {
                this.param1 = param1;
            }

            protected override MyResults RunSpecific()
            {
                return new TestResults(this);
            }

            protected override void ValidateParams()
            {
                if (this.param1 == null)
                    throw new EngineValidationException("param1 cannot be null");
            }

            private class TestResults : MyResults
            {
                public TestResults(MyEngine e) : base(e)
                {
                }

                protected override string GetStringForOutput()
                {
                    StringBuilder sb = new StringBuilder();
                    sb.Append("String for the TestResults results class");
                    return sb.ToString();
                }
            }
        }
    }
}