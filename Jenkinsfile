pipeline {
    agent any
    stages {
        stage('Build') {
            def msbuild = tool 'Main';
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
                
                bat "\"${tool 'MSBuild'}\" SolutionName.sln /p:Configuration=Release /p:Platform=\"Any CPU\" /p:ProductVersion=1.0.0.${env.BUILD_NUMBER}"
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

