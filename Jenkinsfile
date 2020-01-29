#!/usr/bin/env groovy

@Library('lco-shared-libs@0.0.10') _

pipeline {
	agent any
	environment {
		dockerImage = null
		PROJ_NAME = projName()
		GIT_DESCRIPTION = gitDescribe()
		DOCKER_IMG = dockerImageName("${LCO_DOCK_REG}", "${PROJ_NAME}", "${GIT_DESCRIPTION}")
	}
	options {
		timeout(time: 1, unit: 'HOURS')
	}
	stages {
		stage('Build image') {
			steps {
				script {
					dockerImage = docker.build("${DOCKER_IMG}")
				}
			}
		}
		stage('Push image') {
			steps {
				script {
					dockerImage.push()
				}
			}
		}
	    stage('Deploy') {
	        when {
                buildingTag();
	            }
	        steps {
	            withKubeConfig([credentialsId: 'prod-kube-config']) {
	                sh('helm upgrade --install floyds-ogg-en06 helm/ -f helm/ogg-values.yaml --force --set image.tag="${GIT_DESCRIPTION}"')
	                sh('helm upgrade --install floyds-coj-en12 helm/ -f helm/coj-values.yaml --force --set image.tag="${GIT_DESCRIPTION}"')
	            }
	        }
	    }
	}
}

