node {
//def msbuild = tool 'Main';
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

