using System;
using System.Net.Http;
using System.Net.Sockets;
using System.Threading.Tasks;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Support for tests that exercise a live external web service (UniProt, Koina, PRIDE, ...).
    /// Such tests carry <c>[Category("ExternalService")]</c> so CI can run them in a dedicated,
    /// non-blocking job (see .github/workflows/Test.yml) instead of the required unit-test run.
    ///
    /// The hard part is telling two failures apart:
    ///   * the service is unavailable (down, rate-limited, 5xx, timeout) - NOT our bug, so the
    ///     test should be <b>skipped</b> ("we tried, the service is down, don't worry"); versus
    ///   * the service answered but the contract is broken (our URL is wrong, the response no
    ///     longer parses, an expected value is missing) - that is a real regression and must FAIL.
    ///
    /// <see cref="RunAsync"/> wraps a test body and converts availability failures into
    /// <see cref="Assert.Ignore(string)"/> (reported as Skipped) while letting genuine assertion
    /// failures propagate. A test signals "unavailable" from inside the body either by letting a
    /// transport exception bubble up or by calling <see cref="ThrowIfUnavailable"/> on the response.
    /// </summary>
    public static class ExternalServiceTestHelper
    {
        /// <summary>
        /// Runs <paramref name="testBody"/>; if the external service proves unavailable, marks the
        /// test Skipped (via Assert.Ignore) instead of Failed. Real assertion failures propagate.
        /// </summary>
        public static async Task RunAsync(string serviceName, Func<Task> testBody)
        {
            try
            {
                await testBody();
            }
            catch (ExternalServiceUnavailableException e)
            {
                Skip(serviceName, $"unavailable ({e.Message})");
            }
            catch (HttpRequestException e)
            {
                Skip(serviceName, $"unavailable ({e.Message})");
            }
            catch (TaskCanceledException)
            {
                Skip(serviceName, "timed out");
            }
            catch (SocketException e)
            {
                Skip(serviceName, $"unreachable ({e.Message})");
            }
        }

        /// <summary>
        /// Records the skip reason on both the console/CI log (via TestContext.Progress, which is
        /// flushed immediately regardless of verbosity) and the NUnit test result, then skips the
        /// test. This is what surfaces the "we tried, the service is down, don't worry" message.
        /// </summary>
        private static void Skip(string serviceName, string reason)
        {
            string message = $"Skipping external-service test: {serviceName} {reason}. " +
                             "This is a third-party availability problem, not a code failure.";
            TestContext.Progress.WriteLine(message);
            Assert.Ignore(message);
        }

        /// <summary>
        /// Probes an external service's health/readiness URL and skips the calling test (or, from a
        /// [OneTimeSetUp], the whole fixture) if it is unreachable. Use this for tests whose service
        /// call is buried deep in production code and cannot easily route through <see cref="RunAsync"/>.
        /// </summary>
        public static void EnsureReachable(string serviceName, string healthUrl)
        {
            try
            {
                using var client = new HttpClient { Timeout = TimeSpan.FromSeconds(20) };
                var response = client.GetAsync(healthUrl).GetAwaiter().GetResult();
                ThrowIfUnavailable(response);
                if (!response.IsSuccessStatusCode)
                {
                    throw new ExternalServiceUnavailableException($"health check returned HTTP {(int)response.StatusCode}");
                }
            }
            catch (ExternalServiceUnavailableException e)
            {
                Skip(serviceName, $"unavailable ({e.Message})");
            }
            catch (HttpRequestException e)
            {
                Skip(serviceName, $"unavailable ({e.Message})");
            }
            catch (TaskCanceledException)
            {
                Skip(serviceName, "timed out");
            }
            catch (SocketException e)
            {
                Skip(serviceName, $"unreachable ({e.Message})");
            }
        }

        /// <summary>
        /// Classifies an HTTP response as a service-availability problem and, if so, throws
        /// <see cref="ExternalServiceUnavailableException"/> so <see cref="RunAsync"/> skips the test.
        /// Catches transport-level status codes (408/429/5xx) and known error bodies that some
        /// services (e.g. UniProt's streaming endpoint) return with an HTTP 200.
        /// </summary>
        public static void ThrowIfUnavailable(HttpResponseMessage response, string bodyText = null)
        {
            int status = (int)response.StatusCode;
            if (status == 408 || status == 429 || status >= 500)
            {
                throw new ExternalServiceUnavailableException($"HTTP {status} {response.ReasonPhrase}");
            }

            // Some services answer 200 but stream an error message in the body.
            if (!string.IsNullOrWhiteSpace(bodyText) &&
                bodyText.Contains("Error encountered when streaming data", StringComparison.OrdinalIgnoreCase))
            {
                throw new ExternalServiceUnavailableException("service returned an error body instead of data");
            }
        }
    }

    /// <summary>
    /// Marker exception meaning "the external service is unavailable" (as opposed to a real bug).
    /// <see cref="ExternalServiceTestHelper.RunAsync"/> turns this into a skipped test.
    /// </summary>
    public class ExternalServiceUnavailableException : Exception
    {
        public ExternalServiceUnavailableException(string message) : base(message) { }
    }
}
