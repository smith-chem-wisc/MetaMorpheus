## Base image is the Alpine Linux distro with .NET Core runtime
FROM mcr.microsoft.com/dotnet/runtime:8.0-noble-chiseled

## Copies contents of the build folder into the Docker image
ADD CMD/bin/Release/net8.0/ /metamorpheus/

## Installs bash (seemingly necessary for NextFlow)
RUN apk add --no-cache bash

## Installs procps
RUN apk add --no-cache procps

## Set the entrypoint of the Docker image to CMD.dll
ENTRYPOINT ["dotnet", "/MetaMorpheus/CMD.dll"]