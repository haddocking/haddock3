pipeline {
  agent {
    docker {
      image 'python:3'
    }

  }
  stages {
    stage('Install') {
      steps {
        sh '''pip install -r requirements.txt'''
      }
    }
    stage('Test') {
      steps {
        sh '''export HADDOCK3=`pwd`
chmod +x `pwd`/haddock/src/*
python -m coverage run -m unittest discover
export CODECOV_TOKEN=2ef88c60-f7a3-46bb-bf66-fa112b7896b7
codecov'''
      }
    }
  }
}