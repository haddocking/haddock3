pipeline {
  agent {
    docker {
      image 'ubuntu:latest'
    }

  }
  stages {
    stage('Install') {
      steps {
        sh '''pip install codecov
apt-get update
apt-get install gawk'''
      }
    }
    stage('Test') {
      steps {
        sh '''export HADDOCK3=`pwd`
chmod +x `pwd`/haddock/src/*
python -m coverage run -m unittest discover'''
      }
    }
  }
}