## Base image is the Alpine Linux distro with .NET Core runtime
FROM mcr.microsoft.com/dotnet/core/runtime:3.1-alpine AS build

## Copies contents of the build folder into the Docker image
ADD CMD/bin/Release/netcoreapp3.1/ /metamorpheus/

## Set the entrypoint of the Docker image to CMD.dll
ENTRYPOINT ["dotnet", "metamorpheus/CMD.dll"]