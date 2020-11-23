Write-Host "Restoring NuGet Packages" -BackgroundColor Blue
nuget restore MetaMorpheus.sln -verbosity quiet

Write-Host "Building" -BackgroundColor Blue
dotnet msbuild -v:minimal -p:Configuration=Release -p:UseSharedCompilation=false

Write-Host "Building Docker Image" -BackgroundColor Blue
docker build -t smithchemwisc/metamorpheus:dev .

docker run --rm -v C:/Users/rmillikin/Desktop:/mnt/data smithchemwisc/metamorpheus:dev --test -o ./mnt/data/DockerMicrovignetteOutput -v minimal

#Write-Host "Deploying Build to Docker Hub" -BackgroundColor Blue
#docker login docker.io -u smithchemwiscappveyor -p 09e1eb36-31af-404a-becc-ad408ef2150e
#docker push smithchemwisc/metamorpheus:lrproteogenomics