pipeline {
  agent {
    docker {
      image 'continuumio/anaconda3'
    }

  }
  stages {
    stage('Install') {
      steps {
        sh '/opt/conda/bin/conda env create -f haddock3_env.yml'
        sh 'export HADDOCK3=`pwd`'
        sh 'chmod +x `pwd`/haddock/src/*'
        sh 'export CODECOV_TOKEN=2ef88c60-f7a3-46bb-bf66-fa112b7896b7'
      }
    }
    stage('Test') {
      steps {
        sh 'python -m coverage run -m unittest discover'
        sh 'codecov'
      }
    }
  }
}