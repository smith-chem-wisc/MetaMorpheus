FROM microsoft/dotnet:sdk AS build-env
LABEL maintainer="Anthony Cesnik <cesnik@wisc.edu>"
WORKDIR /app
COPY . ./
RUN dotnet restore CMD/CMD.csproj
RUN dotnet publish -c Release -f netcoreapp2.0 CMD/CMD.csproj
WORKDIR /app/CMD/bin/Release/netcoreapp2.0
ENTRYPOINT ["dotnet", "CMD.dll"]