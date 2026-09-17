using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Threading.Tasks;
using EngineLayer;
using MetaMorpheusAvalonia.ViewModels;
using NUnit.Framework;
using TaskLayer;

namespace Test.AvaloniaGui;

/// <summary>
/// The engines report through static events raised from inside Parallel.For, so they arrive on
/// arbitrary thread-pool threads. Getting the marshalling wrong loses them silently, which is worse
/// than not subscribing at all: the progress bar and log look merely quiet rather than broken.
/// </summary>
public class EngineEventTests
{
    /// <summary>
    /// Raises each subscribed event from a pool thread and requires it to land. This is what the old
    /// marshalling failed: with no Avalonia application, it Post()ed to a ManagedDispatcherImpl bound
    /// to whichever thread first touched it, and anything raised elsewhere was dropped.
    /// </summary>
    [Test]
    public async Task EventsRaisedFromAPoolThreadArrive()
    {
        using var viewModel = new MainWindowViewModel();
        int before = viewModel.Log.Length;

        await Task.Run(() =>
        {
            Raise(typeof(MetaMorpheusTask), "LogHandler", new StringEventArgs("log-from-pool-thread", null));
            Raise(typeof(MetaMorpheusEngine), "WarnHandler", new StringEventArgs("engine-warning", null));
            Raise(typeof(MyFileManager), "WarnHandler", new StringEventArgs("file-warning", null));
            Raise(typeof(MetaMorpheusEngine), "OutProgressHandler",
                new ProgressEventArgs(63, "engine-progress", new List<string>()));
            Raise(typeof(MetaMorpheusEngine), "OutLabelStatusHandler", new StringEventArgs("engine-status", null));
        });

        Assert.Multiple(() =>
        {
            Assert.That(viewModel.Log.Length, Is.GreaterThan(before), "nothing reached the log at all");
            Assert.That(viewModel.Log, Does.Contain("log-from-pool-thread"));
            Assert.That(viewModel.Log, Does.Contain("engine-warning"), "engine warnings must reach the log");
            Assert.That(viewModel.Log, Does.Contain("file-warning"), "MyFileManager warnings must reach the log");
            Assert.That(viewModel.Progress, Is.EqualTo(63), "inner-engine progress must move the bar");
            Assert.That(viewModel.Status, Is.EqualTo("engine-status"));
        });
    }

    /// <summary>
    /// Gate C. Every static event on the engine types is either handled or listed as a deliberate
    /// exclusion with a reason, so adding one to the engine fails here rather than being forgotten.
    /// </summary>
    [Test]
    public void EveryEngineEventIsHandledOrExplicitlyExcluded()
    {
        Type[] sources = { typeof(MetaMorpheusTask), typeof(EverythingRunnerEngine), typeof(MetaMorpheusEngine), typeof(MyFileManager) };

        var accountedFor = MainWindowViewModel.SubscribedEventNames
            .Concat(MainWindowViewModel.DeliberatelyUnsubscribedEvents.Keys)
            .ToHashSet();

        var unaccounted = sources
            .SelectMany(t => t.GetEvents(BindingFlags.Public | BindingFlags.Static)
                .Select(e => $"{t.Name}.{e.Name}"))
            .Where(name => !accountedFor.Contains(name))
            .OrderBy(name => name, StringComparer.Ordinal)
            .ToList();

        Assert.That(unaccounted, Is.Empty,
            "these engine events are neither subscribed nor listed as deliberate exclusions: "
            + string.Join(", ", unaccounted));
    }

    /// <summary>The two lists must describe reality, not intent.</summary>
    [Test]
    public void TheDeclaredEventListsNameRealEvents()
    {
        Type[] sources = { typeof(MetaMorpheusTask), typeof(EverythingRunnerEngine), typeof(MetaMorpheusEngine), typeof(MyFileManager) };
        var real = sources
            .SelectMany(t => t.GetEvents(BindingFlags.Public | BindingFlags.Static).Select(e => $"{t.Name}.{e.Name}"))
            .ToHashSet();

        var declared = MainWindowViewModel.SubscribedEventNames
            .Concat(MainWindowViewModel.DeliberatelyUnsubscribedEvents.Keys)
            .ToList();

        Assert.That(declared.Where(name => !real.Contains(name)), Is.Empty,
            "an entry names an event that does not exist");
        Assert.That(declared, Is.Unique);
    }

    /// <summary>
    /// Disposing has to detach every subscription, or a discarded view model keeps handling the next
    /// one's events - which showed up as one test's run writing into another's log.
    /// </summary>
    [Test]
    public void DisposeStopsEventsArriving()
    {
        var viewModel = new MainWindowViewModel();
        viewModel.Dispose();
        int after = viewModel.Log.Length;

        // No Assert.Fail on an empty delegate here: with nothing else subscribed, empty is the point.
        RaiseIfAnyoneIsListening(typeof(MetaMorpheusTask), "LogHandler", new StringEventArgs("after-dispose", null));

        Assert.That(viewModel.Log, Does.Not.Contain("after-dispose"),
            "a disposed view model must not still be listening");
        Assert.That(viewModel.Log.Length, Is.EqualTo(after));
    }

    /// <summary>
    /// Static events can only be raised from their declaring type, so a test has to reach the backing
    /// delegate directly. This is the price of the engines reporting through statics.
    /// </summary>
    private static void Raise(Type declaringType, string eventName, EventArgs args)
    {
        if (!RaiseIfAnyoneIsListening(declaringType, eventName, args))
        {
            Assert.Fail($"{declaringType.Name}.{eventName} has no subscribers, so nothing is listening");
        }
    }

    private static bool RaiseIfAnyoneIsListening(Type declaringType, string eventName, EventArgs args)
    {
        FieldInfo field = declaringType.GetField(eventName, BindingFlags.Static | BindingFlags.NonPublic)
            ?? throw new InvalidOperationException($"no backing field for {declaringType.Name}.{eventName}");

        if (field.GetValue(null) is not Delegate handler)
        {
            return false;
        }
        handler.DynamicInvoke(null, args);
        return true;
    }
}
