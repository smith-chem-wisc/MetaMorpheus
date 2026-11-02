using System.Runtime.CompilerServices;

// GuiGlobalParams and friends deliberately keep internal setters so only the settings layer writes
// them. Splitting this assembly out of GuiFunctions put those writers on the other side of an
// assembly boundary, so grant them access rather than widening the public API.
[assembly: InternalsVisibleTo("GuiFunctions")]
[assembly: InternalsVisibleTo("Test")]
