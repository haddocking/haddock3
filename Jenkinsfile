pipeline {
  agent {
    docker {
      image 'python:3'
    }

  }
  stages {
    stage('Install') {
      steps {
        sh 'pip install -r requirements.txt'
        sh 'chmod +x `pwd`/haddock/src/*'
      }
    }
    stage('Test') {
      steps {
        sh 'python -m coverage run -m unittest discover'
        sh 'codecov'
      }
    }
  }
  environment {
    HADDOCK3 = "${env.WORKSPACE}"
    CNS_EXE = ""
    CODECOV_TOKEN = '2ef88c60-f7a3-46bb-bf66-fa112b7896b7'
  }
}
