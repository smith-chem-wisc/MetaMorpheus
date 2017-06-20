pipeline {
    agent any
    stages {
        stage('Build') {
            steps {
                checkout scm
                echo 'Now Building...'
                bat 'set'
                bat '"E:\\Stefan\\Bin" restore MetaMorpheus.sln'
            }
        }
        
        stage('Test'){
            steps {
                echo 'Now Testing...'
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

