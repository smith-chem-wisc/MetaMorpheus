using FragmentGeneration;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    [TestFixture]
    public class Class1
    {
        [Test]
        public void FirstTest()
        {
            int a = 1;
            Assert.AreEqual(1, a);
            Bin b = new Bin(2);
            Assert.AreEqual(1, a);
        }
    }
}
