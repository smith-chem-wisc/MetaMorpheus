pipeline {
    agent any
    stages {
        stage('Build') {
            steps {
                echo 'Now Building...'
                bat 'set'
                bat 'E:\Stefan\Bin\nuget.exe restore C:\Program Files (x86)\Jenkins\workspace\Test_MetaMorpheus_master-TH6RBEEZ56ABBE24CBJBDWDCE3DM64VH56NCJ5X6ZTGQ2KCFPUWQ\MetaMorpheus.sln'
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

