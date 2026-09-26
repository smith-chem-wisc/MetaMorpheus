using System;
using System.Net;
using System.Net.Http;
using System.Net.Sockets;
using System.Threading.Tasks;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Offline stubs for every classification branch of <see cref="ExternalServiceTestHelper"/>.
    /// The live canaries only reach these branches when a service is actually down, so without
    /// these tests a broken classifier would go unnoticed until the next outage turned CI red.
    /// No test here touches the network.
    /// </summary>
    [TestFixture]
    public static class ExternalServiceTestHelperTests
    {
        private static Task Throwing(Exception e) => Task.FromException(e);

        [Test]
        public static void RunAsync_SkipsOnUnavailableMarker() =>
            Assert.ThrowsAsync<IgnoreException>(() => ExternalServiceTestHelper.RunAsync("Stub",
                () => Throwing(new ExternalServiceUnavailableException("stub outage"))));

        [Test]
        public static void RunAsync_SkipsOnHttpRequestException() =>
            Assert.ThrowsAsync<IgnoreException>(() => ExternalServiceTestHelper.RunAsync("Stub",
                () => Throwing(new HttpRequestException("stub transport failure"))));

        [Test]
        public static void RunAsync_SkipsOnTaskCanceled() =>
            Assert.ThrowsAsync<IgnoreException>(() => ExternalServiceTestHelper.RunAsync("Stub",
                () => Throwing(new TaskCanceledException("stub timeout"))));

        [Test]
        public static void RunAsync_SkipsOnSocketException() =>
            Assert.ThrowsAsync<IgnoreException>(() => ExternalServiceTestHelper.RunAsync("Stub",
                () => Throwing(new SocketException((int)SocketError.HostNotFound))));

        // A genuine contract break must FAIL, not be swallowed as an outage.
        [Test]
        public static void RunAsync_PropagatesUnclassifiedException() =>
            Assert.ThrowsAsync<InvalidOperationException>(() => ExternalServiceTestHelper.RunAsync("Stub",
                () => Throwing(new InvalidOperationException("stub parse failure"))));

        [Test]
        public static void RunAsync_CompletesWhenBodySucceeds() =>
            Assert.DoesNotThrowAsync(() => ExternalServiceTestHelper.RunAsync("Stub", () => Task.CompletedTask));

        [Test]
        [TestCase(HttpStatusCode.RequestTimeout)]
        [TestCase(HttpStatusCode.TooManyRequests)]
        [TestCase(HttpStatusCode.InternalServerError)]
        [TestCase(HttpStatusCode.ServiceUnavailable)]
        public static void ThrowIfUnavailable_ThrowsOnAvailabilityStatus(HttpStatusCode status)
        {
            using var response = new HttpResponseMessage(status);
            Assert.Throws<ExternalServiceUnavailableException>(() => ExternalServiceTestHelper.ThrowIfUnavailable(response));
        }

        // 4xx other than 408/429 means our request is wrong: a real failure, not an outage.
        [Test]
        [TestCase(HttpStatusCode.OK)]
        [TestCase(HttpStatusCode.BadRequest)]
        [TestCase(HttpStatusCode.NotFound)]
        public static void ThrowIfUnavailable_PassesOtherStatus(HttpStatusCode status)
        {
            using var response = new HttpResponseMessage(status);
            Assert.DoesNotThrow(() => ExternalServiceTestHelper.ThrowIfUnavailable(response, "<uniprot>normal body</uniprot>"));
        }

        [Test]
        public static void ThrowIfUnavailable_ThrowsOnHttp200ErrorBody()
        {
            using var response = new HttpResponseMessage(HttpStatusCode.OK);
            Assert.Throws<ExternalServiceUnavailableException>(() =>
                ExternalServiceTestHelper.ThrowIfUnavailable(response, "Error encountered when streaming data. Please try again later."));
        }
    }
}
