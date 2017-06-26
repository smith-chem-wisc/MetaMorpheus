node {
  
     
//def msbuild = tool 'Main'; sadfasdfasdfasdfasdfasdfasdfaasdfs

}

pipeline {
    agent any
    stages {
        stage('Build') {
            
            steps {
                checkout scm
                echo 'Now Building...'
                bat 'set'
                bat '"E:\\Stefan\\Bin\\nuget.exe" restore MetaMorpheus.sln'
               
            }
        }
        
        stage('Test'){
            
            steps {
                echo 'Now Testing...'
                
                bat '"C:\\Windows\\Microsoft.NET\\Framework\\v4.0.30319\\msbuild.exe" MetaMorpheus.sln'
                //bat "\"${tool 'MSBuild'}\" MetaMorpheus.sln /p:Configuration=Release /p:Platform=\"Any CPU\" /p:ProductVersion=1.0.0.${env.BUILD_NUMBER}"
                input "Does the staging environment look ok?"
             }
                  }
     
        stage('Deploy'){
            steps{
                echo 'Now Deploying...' 
        
                }
            }
    
    }
    post{
        always{
             echo("Pipeline is correct")
        }   
    }
        
 }

