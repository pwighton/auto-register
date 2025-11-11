.PHONY: all docker

all: docker

docker:
	docker build -f ./docker/Dockerfile -t pwighton/areg .
