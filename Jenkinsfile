pipeline {
    agent any
    stages {
        stage('Build') {
            steps {
                echo 'Now Building...'
                bat 'set'
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
    
    
    post{
        always{
             echo("Pipeline is correct")
        }
    }
        
}
