using NUnit.Framework;
using TaskLayer;
using Moq;

namespace Test
{
    [TestFixture]
    internal class DbForTaskTests
    {

        [Test]
        public void MockDbForTask_IsSpectralLibrary_ReturnsTrue()
        {
            var dbMock = new Mock<IdbForTask>();
            dbMock.Setup(x => x.FilePath).Returns("path");
            dbMock.Setup(x => x.IsContaminant).Returns(false);
            dbMock.Setup(x => x.FileName).Returns("name");
            dbMock.Setup(x => x.IsSpectralLibrary).Returns(true);

            Assert.That(dbMock.Object.IsSpectralLibrary, Is.True);
        }
    }
}